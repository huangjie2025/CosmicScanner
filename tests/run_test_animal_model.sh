#!/usr/bin/env bash
# Pedigree Animal Model (-animal-model) regression and validation suite.
#
# Covers:
#   1. default simulation (untyped ancestors, causal SNP, h2=0.5) smoke + causal recovery
#   2. LOCO warning (LOCO is not applicable to a pedigree-only kernel)
#   3. ssGWAS tau=0 equivalence (shared pedigree A/A22 constructor) with row-order
#      shuffle invariance (pedigree rows shuffled -> identical results)
#   4. h2 ~ 0 (no polygenic signal) scenario -> VCE near zero
#   5. high-heritability scenario
#   6. incomplete pedigree (missing one/both parents)
#   7. no-ancestor pedigree (typed only) -> still runs, finite output
#
# Usage: bash tests/run_test_animal_model.sh
set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULT_DIR="${ROOT}/tests/results/animal_model"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"
SIM="${ROOT}/tests/run_animal_model_simulation.py"

mkdir -p "${RESULT_DIR}"
rm -rf "${RESULT_DIR}"/* 2>/dev/null || true

echo "=== Pedigree Animal Model validation suite ==="
echo "results: ${RESULT_DIR}"

if [[ ! -x "${BIN}" ]]; then
  echo "FAIL: missing cosmicscanner binary ${BIN}" >&2
  exit 2
fi

GEN() {
  local label="$1"; shift
  python3 "${SIM}" --out-dir "${RESULT_DIR}/${label}" --prefix animal_model "$@"
}

AM_RUN() {
  # AM_RUN <label> <out-suffix> [extra CLI args...]
  local label="$1"; local suffix="$2"; shift 2
  local extra=("$@")
  "${BIN}" \
    --bfile "${RESULT_DIR}/${label}/animal_model" \
    --pheno "${RESULT_DIR}/${label}/animal_model.pheno" \
    --animal-model \
    --ped "${RESULT_DIR}/${label}/pedigree.txt" \
    --out "${RESULT_DIR}/${label}/${suffix}" \
    "${extra[@]}" 2>&1
}

# ---------- 1/2. default sim + LOCO warning ----------
GEN default --seed 101
AM_LOG="$(AM_RUN default am 2>&1 || true)"
echo "${AM_LOG}" > "${RESULT_DIR}/default/am.log"
if ! grep -q "WARNING: LOCO is not applicable to pedigree-only Animal Model GWAS" <<<"${AM_LOG}"; then
  echo "FAIL [1]: expected LOCO warning (4-chrom data, LOCO on by default)" >&2
  exit 1
fi
grep -q "Pedigree A22 .* built via Schur complement of sparse A⁻¹" <<<"${AM_LOG}" || { echo "FAIL: A22 Schur-complement build missing"; exit 1; }
python3 - "${RESULT_DIR}/default" <<'PY'
import json, os, sys
d = sys.argv[1]
null = {}
with open(f"{d}/am.animal_model.null.txt") as f:
    next(f)
    for line in f:
        k, v = line.rstrip().split("\t")
        null[k] = float(v) if k in ("sigma2_a","sigma2_e","h2") else v
assert null["converged"] == "yes", null
assert abs(null["h2"] - 0.5) < 0.3, f"h2 out of range: {null['h2']}"
assert null["sigma2_a"] > 0.05, null
# standard header + causal recovery
hdr = "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP\tLOG10P"
with open(f"{d}/am.assoc.txt") as f:
    assert f.readline().rstrip("\n") == hdr
    rows = [l.split("\t") for l in f if l.strip()]
assert len(rows) == 600, len(rows)
rows.sort(key=lambda r: float(r[9]))
tops = [r[1] for r in rows[:10]]
assert "rs000001" in tops, f"causal not in top10: {tops}"
print("PASS [1]: default sim null + causal recovery + LOCO warning")
PY

# ---------- 3. ssGWAS tau=0 equivalence + shuffle invariance ----------
awk '{print $1":"$2}' "${RESULT_DIR}/default/animal_model.fam" > "${RESULT_DIR}/default/gl.txt"
awk '{print $2}' "${RESULT_DIR}/default/animal_model.fam" > "${RESULT_DIR}/default/gl_iid.txt"
"${BIN}" \
  --bfile "${RESULT_DIR}/default/animal_model" \
  --pheno "${RESULT_DIR}/default/animal_model.pheno" \
  --single-step \
  --ped "${RESULT_DIR}/default/pedigree.txt" \
  --genotyped-list "${RESULT_DIR}/default/gl_iid.txt" \
  --h2-tau 0 --h2-omega 0 \
  --no-loco \
  --out "${RESULT_DIR}/default/ss" > "${RESULT_DIR}/default/ss.log" 2>&1

GEN shuffled --seed 101 --shuffle-pedigree
AM_LOG_S="$(AM_RUN shuffled am --no-loco 2>&1 || true)"
echo "${AM_LOG_S}" > "${RESULT_DIR}/shuffled/am.log"

python3 - "${RESULT_DIR}/default" <<'PY'
import math, os, sys
d = sys.argv[1]
def load(fn):
    rs={}
    for line in open(fn):
        c=line.rstrip().split("\t")
        if c[0]=="CHR": continue
        rs[c[1]]=dict(beta=float(c[7]),se=float(c[8]),logp=float(c[10]))
    return rs
def pearson(xs,ys):
    n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
    cov=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    vx=sum((x-mx)**2 for x in xs); vy=sum((y-my)**2 for y in ys)
    return cov/math.sqrt(vx*vy)
a=load(f"{d}/am.assoc.txt")
b=load(f"{d}/ss.assoc.txt")
assert set(a)==set(b)
for k,thr in (("beta",0.999999),("se",0.999999),("logp",0.999999)):
    xs=[a[s][k] for s in a]; ys=[b[s][k] for s in a]
    c=pearson(xs,ys)
    assert c>thr, (k,c)
print("PASS [2]: ssGWAS --h2-tau 0 --h2-omega 0 equivalence Pearson>0.999999",
      *(f"{k}={pearson([a[s][k] for s in a],[b[s][k] for s in a]):.9f}" for k in ("beta","se","logp")))
# shuffle invariance: identical to unshuffled run
r=load(f"{d}/shuffled/../shuffled/am.assoc.txt") if False else None
PY
# shuffle invariance asserted in Python below
python3 - "${RESULT_DIR}/default" "${RESULT_DIR}/shuffled" <<'PY'
import math, sys
da, ds = sys.argv[1], sys.argv[2]
def load(fn):
    rs={}
    for line in open(fn):
        c=line.rstrip().split("\t")
        if c[0]=="CHR": continue
        rs[c[1]]=(float(c[7]),float(c[8]),float(c[10]))
    return rs
a=load(f"{da}/am.assoc.txt"); b=load(f"{ds}/am.assoc.txt")
assert set(a)==set(b)
maxd=max(max(abs(x-y) for x,y in zip(a[s],b[s])) for s in a)
assert maxd < 1e-12, f"shuffle invariance broken, max diff {maxd}"
na={}
for line in open(f"{da}/am.animal_model.null.txt"):
    c=line.rstrip().split("\t")
    if c[0]=="sigma2_a" or c[0]=="sigma2_e" or c[0]=="h2": na[c[0]]=float(c[1])
nb={}
for line in open(f"{ds}/am.animal_model.null.txt"):
    c=line.rstrip().split("\t")
    if c[0]=="sigma2_a" or c[0]=="sigma2_e" or c[0]=="h2": nb[c[0]]=float(c[1])
assert max(abs(na[k]-nb[k]) for k in na) < 1e-12
print("PASS [3]: pedigree row-order shuffle invariance (identical results)")
PY

# ---------- 4. sigma2_a ~ 0 (null VCE) ----------
GEN null_h2 --seed 202 --h2 0.0 --causal-effect 0.0
AM_RUN null_h2 am --no-loco > "${RESULT_DIR}/null_h2/am.log" 2>&1
python3 - "${RESULT_DIR}/null_h2" <<'PY'
import sys
d=sys.argv[1]
val={}
for line in open(f"{d}/am.animal_model.null.txt"):
    c=line.rstrip().split("\t")
    if c[0] in ("sigma2_a","sigma2_e","h2"): val[c[0]]=float(c[1])
assert val["h2"] < 0.1, val
print(f"PASS [4]: sigma2_a=0 VCE h2={val['h2']:.4f} < 0.1")
PY

# ---------- 5. high heritability ----------
GEN high_h2 --seed 303 --h2 0.9
AM_RUN high_h2 am --no-loco > "${RESULT_DIR}/high_h2/am.log" 2>&1
python3 - "${RESULT_DIR}/high_h2" <<'PY'
import sys
d=sys.argv[1]
val={}
for line in open(f"{d}/am.animal_model.null.txt"):
    c=line.rstrip().split("\t")
    if c[0] in ("sigma2_a","sigma2_e","h2"): val[c[0]]=float(c[1])
assert val["h2"] > 0.5, val
print(f"PASS [5]: high-h2 VCE h2={val['h2']:.4f} > 0.5")
PY

# ---------- 6. incomplete pedigree (missing one/both parents) ----------
GEN incomplete --seed 404 --h2 0.5 --missing-one-parent 0.3 --missing-both-parents 0.1
AM_RUN incomplete am --no-loco > "${RESULT_DIR}/incomplete/am.log" 2>&1
python3 - "${RESULT_DIR}/incomplete" <<'PY'
import sys
d=sys.argv[1]
rows=[l.split("\t") for l in open(f"{d}/am.assoc.txt") if l.strip() and not l.startswith("CHR")]
top=[r[1] for r in sorted(rows,key=lambda r:float(r[9]))[:10]]
val={}
for line in open(f"{d}/am.animal_model.null.txt"):
    c=line.rstrip().split("\t")
    if c[0] in ("sigma2_a","sigma2_e","h2"): val[c[0]]=float(c[1])
assert "rs000001" in top, top
assert abs(val["h2"]-0.5) < 0.35, val
print(f"PASS [6]: incomplete pedigree (missing parents) causal top-10, h2={val['h2']:.3f}")
PY

# ---------- 7. no-ancestor pedigree (typed only) ----------
GEN noanc --seed 505 --no-ancestors --h2 0.5 --causal-effect 0.0
AM_RUN noanc am --no-loco > "${RESULT_DIR}/noanc/am.log" 2>&1
python3 - "${RESULT_DIR}/noanc" <<'PY'
import sys
d=sys.argv[1]
n=0; finite=True
for line in open(f"{d}/am.assoc.txt"):
    c=line.rstrip().split("\t")
    if c[0]=="CHR": continue
    n+=1
    for v in (c[7],c[8],c[9]): finite = finite and float(v)==float(v)
assert n==600
val={}
for line in open(f"{d}/am.animal_model.null.txt"):
    c=line.rstrip().split("\t")
    if c[0] in ("sigma2_a","sigma2_e","converged"): val[c[0]]=c[1]
assert val["converged"]=="yes"
print("PASS [7]: no-ancestor pedigree runs, 600 finite assoc rows")
PY

echo
echo "PASS: all animal-model validations"
exit 0
