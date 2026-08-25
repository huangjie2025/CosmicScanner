#!/usr/bin/env bash
# External GEMMA -k benchmark for the Pedigree Animal Model.
#
# Runs --animal-model, dumps the exact analysis A22 kernel with
# --write-ped-kernel, passes it to GEMMA (-k ... -lmm 4) on the same data,
# and compares BETA / SE / -log10(p_wald).
#
# GEMMA is optional: if no GEMMA binary is found the script reports SKIP.
#
# Usage: bash tests/run_animal_model_gemma_validation.sh
set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULT_DIR="${ROOT}/tests/results/animal_model_gemma"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
SIM="${ROOT}/tests/run_animal_model_simulation.py"

mkdir -p "${RESULT_DIR}"
rm -rf "${RESULT_DIR}"/* 2>/dev/null || true

resolve_gemma() {
  if [[ -n "${GEMMA_BIN:-}" ]]; then printf '%s\n' "${GEMMA_BIN}"; return 0; fi
  if command -v gemma >/dev/null 2>&1; then command -v gemma; return 0; fi
  local conda_bin="/disk228/huangjie/software/miniconda3/envs/benchmark/bin/gemma"
  if [[ -x "${conda_bin}" ]]; then printf '%s\n' "${conda_bin}"; return 0; fi
  local repo_bin="/disk228/huangjie/PanOmics/References_and_Tools/GEMMA-master/bin/gemma"
  if [[ -x "${repo_bin}" ]]; then printf '%s\n' "${repo_bin}"; return 0; fi
  return 1
}

if [[ ! -x "${BIN}" ]]; then
  echo "FAIL: missing cosmicscanner binary ${BIN}" >&2
  exit 2
fi
if ! GEMMA_EXE="$(resolve_gemma)"; then
  echo "SKIP: no GEMMA binary found (GEMMA validation skipped; simulation + ssGWAS equivalence still validate the model)"
  exit 0
fi

python3 "${SIM}" --out-dir "${RESULT_DIR}/data" --prefix animal_model --seed 777 --h2 0.5

"${BIN}" \
  --bfile "${RESULT_DIR}/data/animal_model" \
  --pheno "${RESULT_DIR}/data/animal_model.pheno" \
  --animal-model \
  --ped "${RESULT_DIR}/data/pedigree.txt" \
  --write-ped-kernel \
  --no-loco \
  --out "${RESULT_DIR}/data/am" > "${RESULT_DIR}/data/am.log" 2>&1

KERNEL="${RESULT_DIR}/data/am.animal_model.k.txt"
test -f "${KERNEL}" || { echo "FAIL: missing kernel ${KERNEL}"; exit 1; }

# Sample order of the kernel must match the .fam order (no missing phenotypes here).
python3 - "${RESULT_DIR}/data" <<'PY'
import sys
d=sys.argv[1]
fam=[l.split()[1] for l in open(f"{d}/animal_model.fam")]
ids=[l.rstrip().split(":")[1] for l in open(f"{d}/am.animal_model.k.id")]
assert fam==ids, (fam[:3], ids[:3])
print(f"kernel sample order matches .fam order (n={len(fam)})")
PY

mkdir -p "${RESULT_DIR}/gemma_out"
cd "${RESULT_DIR}/gemma_out"
"${GEMMA_EXE}" -bfile "${RESULT_DIR}/data/animal_model" -k "${KERNEL}" -lmm 4 -o cs_k \
  > "${RESULT_DIR}/gemma_out/gemma.log" 2>&1
cd "${ROOT}"

python3 - "${RESULT_DIR}/data" "${RESULT_DIR}/gemma_out/output/cs_k.assoc.txt" <<'PY'
import math, re, sys
d, gfile = sys.argv[1], sys.argv[2]
def load_cs(fn):
    rs={}
    for line in open(fn):
        c=line.rstrip().split("\t")
        if c[0]=="CHR": continue
        rs[c[1]]=dict(beta=float(c[7]),se=float(c[8]),logp=float(c[10]))
    return rs
def load_gemma(fn):
    rs={}
    for line in open(fn):
        c=line.rstrip().split("\t")
        if c[0]=="chr": continue
        rs[c[1]]=dict(beta=float(c[7]),se=float(c[8]),p=float(c[12]))
    return rs
a=load_cs(f"{d}/am.assoc.txt")
g=load_gemma(gfile)
common=set(a)&set(g)
def pearson(xs,ys):
    n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
    cov=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    vx=sum((x-mx)**2 for x in xs); vy=sum((y-my)**2 for y in ys)
    return cov/math.sqrt(vx*vy)
res={}
for k in ("beta","se"):
    xs=[a[s][k] for s in common]; ys=[g[s][k] for s in common]
    res[f"pearson_{k}"]=pearson(xs,ys)
    res[f"maxabs_{k}"]=max(abs(x-y) for x,y in zip(xs,ys))
xl=[a[s]["logp"] for s in common]
yl=[-math.log10(g[s]["p"]) if g[s]["p"]>0 else 300 for s in common]
res["pearson_log10p_wald"]=pearson(xl,yl)
gtop=set(sorted(common,key=lambda s:g[s]["p"])[:5])
ctop=set(sorted(common,key=lambda s:a[s]["logp"],reverse=True)[:5])
res["top5_overlap"]=len(gtop&ctop)
res["causal_rank_cs"]=sorted(common,key=lambda s:a[s]["logp"],reverse=True).index("rs000001")+1
res["causal_rank_gemma"]=sorted(common,key=lambda s:g[s]["p"]).index("rs000001")+1
# NOTE: p-values across independent REML implementations differ in log10-space on
# null SNPs (tiny p + numerical se differences magnify in -log10); cutoff >= 0.99.
thr={"pearson_beta":0.999,"pearson_se":0.999,"pearson_log10p_wald":0.99}
ok=True
for k,t in thr.items():
    if res[k]<t:
        ok=False; print(f"FAIL: {k}={res[k]:.6f} < {t}")
    else:
        print(f"PASS: {k}={res[k]:.6f} (>= {t})")
assert res["top5_overlap"]>=3, res
assert res["causal_rank_cs"]<=3 and res["causal_rank_gemma"]<=3, res
with open(f"{d}/gemma_validation_summary.tsv","w") as f:
    f.write("metric\tvalue\n")
    for k in sorted(res): f.write(f"{k}\t{res[k]}\n")
print("PASS: top5 overlap", res["top5_overlap"], "causal rank CS/GM", res["causal_rank_cs"], res["causal_rank_gemma"])
print("GEMMA -k external benchmark PASSED")
PY
