#!/bin/bash
# Common-environment GWAS test: simulated quantitative trait with
# direct additive genetic (GRM), common-environment (group) and residual
# variance components, plus a causal SNP.
set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="${ROOT}/build/cosmicscanner"
DATA_DIR="${ROOT}/tests/data/common_env_simulation"
RESULT_DIR="${ROOT}/tests/results/common_env_simulation"
mkdir -p "${DATA_DIR}" "${RESULT_DIR}"
rm -f "${DATA_DIR}"/* "${RESULT_DIR}"/*

echo "=== Step 1: simulate genotypes + records + groups ==="
python3 - <<'PY'
import numpy as np
from pathlib import Path

rng = np.random.default_rng(20240815)
data_dir = Path("/disk228/huangjie/PanOmics/CosmicScanner/tests/data/common_env_simulation")

# Parameters
N = 120          # subjects
M1, M2 = 100, 60 # SNPs on chr1 / chr2
M = M1 + M2
G_GROUPS = 12    # common-environment groups (10 subjects each)
SIG_A2, SIG_C2, SIG_E2 = 0.45, 0.25, 0.30
BETA_CAUSAL = 0.7
CAUSAL = 7       # SNP index (chr1)

# --- genotypes: 0/1/2 dosage, MAF ~ U[0.08, 0.45]
maf = rng.uniform(0.08, 0.45, M)
geno = (rng.random((N, M)) < maf).astype(np.float64)          # one allele
geno += (rng.random((N, M)) < maf).astype(np.float64)         # two alleles
p = geno.mean(axis=0) / 2.0
# center/scale (VanRaden) -> G = W W' / m, diag ~ 1
W = (geno - 2.0 * p) / np.sqrt(2.0 * p * (1.0 - p) + 1e-12)
G = (W @ W.T) / M + 1e-6 * np.eye(N)

# --- causal SNP on chr1
x_causal = geno[:, CAUSAL] - 2.0 * p[CAUSAL]

# --- random effects
u_a = rng.multivariate_normal(np.zeros(N), SIG_A2 * G, method="cholesky")
u_c = rng.normal(0.0, np.sqrt(SIG_C2), G_GROUPS)

# --- records: subject i has 2 + (i % 3) records, group = i // 10
records = []
for i in range(N):
    n_rec = 2 + (i % 3)
    grp = i // (N // G_GROUPS)
    for r in range(n_rec):
        y = 10.0 + u_a[i] + u_c[grp] + BETA_CAUSAL * x_causal[i] \
            + rng.normal(0.0, np.sqrt(SIG_E2))
        records.append((i, grp, y))
n_rec = len(records)
print(f"N={N} M={M} records={n_rec} groups={G_GROUPS}")

# --- write PLINK files (SNP-major bed)
bim_rows = []
for c, s, base in ((1, M1, 1), (2, M2, M1 + 1)):
    for k in range(s):
        idx = base - 1 + k
        a1 = "A" if idx % 2 == 0 else "C"
        a2 = "C" if idx % 2 == 0 else "A"
        bim_rows.append((str(c), f"rs{idx + 1}", "0", str((k + 1) * 10000), a1, a2))
with open(data_dir / "sim.bim", "w") as f:
    for row in bim_rows:
        f.write("\t".join(row) + "\n")
with open(data_dir / "sim.fam", "w") as f:
    for i in range(N):
        f.write(f"F{i} I{i} 0 0 0 -9\n")

bed = bytearray(b"\x6c\x1b\x01")
n_bps = (N + 3) // 4
for v in range(M):
    block = bytearray(n_bps)
    for i in range(N):
        code = {0: 0b00, 1: 0b10, 2: 0b11}[int(geno[i, v])]
        block[i // 4] |= code << ((i % 4) * 2)
    bed.extend(block)
with open(data_dir / "sim.bed", "wb") as f:
    f.write(bed)

# --- long pheno: FID IID VALUE (mpheno 1)
with open(data_dir / "sim.pheno.long", "w") as f:
    for (i, grp, y) in records:
        f.write(f"F{i} I{i} {y:.6f}\n")

# --- subject-level common-env file: FID IID GROUP
with open(data_dir / "sim.common_env.subject.txt", "w") as f:
    for i in range(N):
        grp = i // (N // G_GROUPS)
        f.write(f"F{i} I{i} G{grp}\n")

# --- record-level common-env file: FID IID GROUP (one line per pheno record,
#     matching pheno file order)
with open(data_dir / "sim.common_env.record.txt", "w") as f:
    for (i, grp, y) in records:
        f.write(f"F{i} I{i} G{grp}\n")

# --- ground truth summary
np.savez(data_dir / "truth.npz",
         sigma_a2=SIG_A2, sigma_c2=SIG_C2, sigma_e2=SIG_E2,
         causal_index=CAUSAL, causal_maf=p[CAUSAL], n_rec=n_rec)
PY

echo "=== Step 2: run common-environment GWAS (no-LOCO) ==="
"${BIN}" \
  --bfile "${DATA_DIR}/sim" \
  --pheno "${DATA_DIR}/sim.pheno.long" \
  --common-env-model \
  --common-env-file "${DATA_DIR}/sim.common_env.subject.txt" \
  --out "${RESULT_DIR}/ce_global" \
  --no-loco \
  --threads 1 > "${RESULT_DIR}/ce_global.log" 2>&1
tail -25 "${RESULT_DIR}/ce_global.log"

echo "=== Step 3: run common-environment GWAS (LOCO) ==="
"${BIN}" \
  --bfile "${DATA_DIR}/sim" \
  --pheno "${DATA_DIR}/sim.pheno.long" \
  --common-env-model \
  --common-env-file "${DATA_DIR}/sim.common_env.subject.txt" \
  --out "${RESULT_DIR}/ce_loco" \
  --threads 1 > "${RESULT_DIR}/ce_loco.log" 2>&1
tail -12 "${RESULT_DIR}/ce_loco.log"

echo "=== Step 4: run record-level common-env mapping ==="
"${BIN}" \
  --bfile "${DATA_DIR}/sim" \
  --pheno "${DATA_DIR}/sim.pheno.long" \
  --common-env-model \
  --common-env-file "${DATA_DIR}/sim.common_env.record.txt" \
  --common-env-record-level \
  --out "${RESULT_DIR}/ce_record_level" \
  --no-loco \
  --threads 1 > "${RESULT_DIR}/ce_record_level.log" 2>&1
tail -12 "${RESULT_DIR}/ce_record_level.log"

echo "=== Step 5: run with fixed variance components ==="
"${BIN}" \
  --bfile "${DATA_DIR}/sim" \
  --pheno "${DATA_DIR}/sim.pheno.long" \
  --common-env-model \
  --common-env-file "${DATA_DIR}/sim.common_env.subject.txt" \
  --out "${RESULT_DIR}/ce_fixedvc" \
  --no-loco \
  --reml-vg 0.45 --reml-vc 0.25 --reml-ve 0.30 \
  --threads 1 > "${RESULT_DIR}/ce_fixedvc.log" 2>&1
tail -8 "${RESULT_DIR}/ce_fixedvc.log"

echo "=== Step 6: validate ==="
python3 - <<'PY'
import numpy as np
from pathlib import Path

res_dir = Path("/disk228/huangjie/PanOmics/CosmicScanner/tests/results/common_env_simulation")
data_dir = Path("/disk228/huangjie/PanOmics/CosmicScanner/tests/data/common_env_simulation")
truth = dict(np.load(data_dir / "truth.npz"))
fails = []

# --- VCE recovery (global run)
null = {}
for line in (res_dir / "ce_global.common_env.null.txt").read_text().splitlines():
    parts = line.split("\t")
    if len(parts) == 2:
        null[parts[0]] = parts[1]
sa = float(null["sigma2_a"]); sc = float(null["sigma2_common"]); se = float(null["sigma2_e"])
print(f"VCE: sigma2_a={sa:.4f} (true {truth['sigma_a2']}), "
      f"sigma2_c={sc:.4f} (true {truth['sigma_c2']}), "
      f"sigma2_e={se:.4f} (true {truth['sigma_e2']})")
if abs(sa - truth["sigma_a2"]) > 0.20: fails.append(f"sigma2_a off: {sa}")
if abs(sc - truth["sigma_c2"]) > 0.15: fails.append(f"sigma2_c off: {sc}")
if abs(se - truth["sigma_e2"]) > 0.15: fails.append(f"sigma2_e off: {se}")
n_rec = int(null["n_records"]); n_sub = int(null["n_subjects"])
if n_rec != int(truth["n_rec"]): fails.append(f"n_records {n_rec} != {int(truth['n_rec'])}")
if n_sub != 120: fails.append(f"n_subjects {n_sub} != 120")

# --- causal signal (global run)
rows = [l.split() for l in (res_dir / "ce_global.assoc.txt").read_text().splitlines()]
data = [r for r in rows[1:] if len(r) >= 10 and r[0] not in ("#CHR",)]
ps = [(i, float(r[9])) for i, r in enumerate(data) if r[9] not in ("NA", "nan")]
if len(data) != 160: fails.append(f"assoc rows {len(data)} != 160")
causal_rows = [r for r in data if r[1] == "rs8"]
if not causal_rows:
    fails.append("causal SNP rs8 missing from output")
else:
    cr = causal_rows[0]
    p_causal = float(cr[9])
    print(f"causal SNP: rs8 beta={cr[7]} se={cr[8]} p={p_causal:.3e}")
    if p_causal > 1e-3: fails.append(f"causal SNP p={p_causal} not significant")
    best = min(ps, key=lambda x: x[1])
    if best[0] != 7: fails.append(f"causal SNP not the top hit (top={best})")
    else: print("causal SNP is the top hit")
n_finite = sum(1 for _, p in ps if np.isfinite(p))
if n_finite < 150: fails.append(f"only {n_finite}/160 finite p-values")

# --- fixed VC run used the supplied components
null2 = {}
for line in (res_dir / "ce_fixedvc.common_env.null.txt").read_text().splitlines():
    parts = line.split("\t")
    if len(parts) == 2:
        null2[parts[0]] = parts[1]
if abs(float(null2["sigma2_a"]) - 0.45) > 1e-9: fails.append("fixed sigma2_a not honored")
if abs(float(null2["sigma2_common"]) - 0.25) > 1e-9: fails.append("fixed sigma2_c not honored")
if abs(float(null2["sigma2_e"]) - 0.30) > 1e-9: fails.append("fixed sigma2_e not honored")

# --- LOCO and record-level runs share the same data surface
for tag, base in (("LOCO", "ce_loco"), ("record-level", "ce_record_level")):
    rows2 = [l.split() for l in (res_dir / f"{base}.assoc.txt").read_text().splitlines()]
    d2 = [r for r in rows2[1:] if len(r) >= 10 and r[0] not in ("#CHR",)]
    if len(d2) != 160: fails.append(f"{tag}: assoc rows {len(d2)} != 160")
    cr2 = [r for r in d2 if r[1] == "rs8"]
    if not cr2: fails.append(f"{tag}: causal SNP missing")
    else:
        p2 = float(cr2[0][9])
        print(f"{tag}: causal p={p2:.3e}")
        if p2 > 1e-2: fails.append(f"{tag}: causal p={p2} degraded")

if fails:
    print("FAILURES:")
    for f_ in fails:
        print("  -", f_)
    raise SystemExit(1)
print("ALL COMMON-ENVIRONMENT CHECKS PASSED")
PY

echo "=== common-environment simulation test PASSED ==="