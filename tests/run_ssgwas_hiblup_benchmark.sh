#!/usr/bin/env bash
# Phase 6.6 external ssGBLUP benchmark vs HIBLUP (reference implementation).
#   bash tests/run_ssgwas_hiblup_benchmark.sh [--build]
#
# Model (identical on both sides, intercept-only, fixed VCs):
#   y = mu + a + e,  a ~ N(0, H*Va),
#   H^-1 = A^-1 + [tau*Gw^-1 - omega*A22^-1], tau=omega=1,
#   Gw = 0.95*G* + 0.05*A22,  G* = b*G + a (tuneG, mean diag/off-diag of A22),
#   G = M M'/S, S = sum_j 2p_j(1-p_j), M = dosage - mean (full genotyped set).
# VCs: Va=0.616532, Ve=1.85018 (HIBLUP AI-REML estimates on the smoke data).
# HIBLUP applies an internal factor to --vc-prior (effective lambda = 0.98 x
# provided value at N=2000); --vc-prior is therefore given as Ve/0.98 so that
# both tools solve the same effective lambda (see benchmark summary).
set -euo pipefail
cd "$(dirname "$0")/.."

BIN=${BIN:-build/cosmicscanner}
HIBLUP=${HIBLUP:-/disk228/huangjie/software/miniconda3/envs/benchmark/bin/hiblup}
SMP=${SMP:-/disk228/huangjie/PanOmics/CosmicSolver/tests/results/ssgblup_smoke}
OUT=${OUT:-tests/results/ssgwas_hiblup_benchmark}

VA=0.616532
VE=1.85018
LAM=$(python3 -c "print('%.15f' % ($VE/$VA))")
VE_HIB=$(python3 -c "print('%.8f' % ($VE/0.98))")

if [[ "${1:-}" == "--build" ]]; then
  make -C build cosmicscanner -j4
fi
mkdir -p "$OUT"

# ---- input data (from the CosmicSolver smoke simulator) ----
cp "$SMP/sim.ped" "$SMP/sim.phe" "$SMP/sim.bed" "$SMP/sim.bim" "$SMP/sim.fam" "$OUT/"
tail -n +2 "$OUT/sim.ped" > "$OUT/sim.cs.ped"

# 4-class phenotype masks: typed+pheno, typed no-pheno, untyped+pheno, untyped no-pheno
python3 - "$OUT" <<'EOF'
import sys
D = sys.argv[1]
lines = [l.split() for l in open(f"{D}/sim.phe")][1:]
keep = set()
for r in lines:
    i = int(r[0])
    if 1 <= i <= 200: keep.add(i)        # founders untyped, phenotyped
    if 401 <= i <= 1200: keep.add(i)     # typed + phenotyped
    if 1401 <= i <= 1700: keep.add(i)    # untyped gen2 + phenotyped
rows = [r for r in lines if int(r[0]) in keep]
with open(f"{D}/sim_m4.phe", "w") as f:
    f.write("FID IID SEX BATCH Pheno_1\n")
    for r in rows: f.write(f"{r[0]} {r[0]} {r[3]} {r[4]} {r[5]}\n")
with open(f"{D}/sim_m4_y.phe", "w") as f:
    f.write("FID IID Pheno_1\n")
    for r in rows: f.write(f"{r[0]} {r[0]} {r[5]}\n")
print(f"4-class pheno: {len(rows)} records")
EOF

# ---- HIBLUP reference: MME with fixed VCs (aligned effective lambda) ----
/usr/bin/time -f "hiblup mme wall: %e s" "$HIBLUP" --pedigree "$OUT/sim.ped" --bfile "$OUT/sim" \
  --pheno "$OUT/sim_m4.phe" --pheno-pos 5 --mme --vc-prior "$VA,$VE_HIB" \
  --threads 8 --out "$OUT/hi_m4_adj" 2>&1 | tail -1

# ---- CosmicScanner reference-compatible scenario (tuneG + A22 blend) ----
mkdir -p "$OUT/cs_m4"
/usr/bin/time -f "cosmicscanner wall: %e s" "$BIN" --bfile "$OUT/sim" --pheno "$OUT/sim_m4_y.phe" \
  --gwas --single-step --single-step-method backsolve --ped "$OUT/sim.cs.ped" \
  --reml-vg "$VA" --reml-ve "$VE" --no-loco --ssgwas-window-snps 10 --ssgwas-solver ldlt \
  --ssgwas-tune-g --ssgwas-blend-a22 --out "$OUT/cs_m4/m2" 2>&1 | tail -1

# ---- pedigree row-order shuffle invariance (CosmicScanner) ----
python3 - "$OUT" <<'EOF'
import random, sys
D = sys.argv[1]
lines = open(f"{D}/sim.cs.ped").read().splitlines()
random.seed(7); random.shuffle(lines)
open(f"{D}/sim_shuf.ped", "w").write("\n".join(lines) + "\n")
EOF
mkdir -p "$OUT/cs_shuf"
"$BIN" --bfile "$OUT/sim" --pheno "$OUT/sim_m4_y.phe" --gwas \
  --single-step --single-step-method backsolve --ped "$OUT/sim_shuf.ped" \
  --reml-vg "$VA" --reml-ve "$VE" --no-loco --ssgwas-window-snps 10 --ssgwas-solver ldlt \
  --ssgwas-tune-g --ssgwas-blend-a22 --out "$OUT/cs_shuf/m2" >/dev/null 2>&1

# ---- comparison / summary ----
python3 - "$OUT" "$VA" "$VE" <<'EOF'
import sys, numpy as np
D, VA, VE = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
def rand(p): return {l.split()[0]: float(l.split()[1]) for l in open(p) if not l.startswith("ID")}
def ebv(p):
    return {l.split()[0]: (float(l.split()[1]), int(l.split()[2]), int(l.split()[3]))
            for l in open(p) if not l.startswith("ID")}
hib = rand(f"{D}/hi_m4_adj.rand")
cs = ebv(f"{D}/cs_m4/m2.ss_backsolve.ebv.txt")
ids = sorted(hib, key=int)
h = np.array([hib[i] for i in ids]); a = np.array([cs[i][0] for i in ids])
G1 = np.array([cs[i][1] for i in ids]); P1 = np.array([cs[i][2] for i in ids])
def rep(x, y):
    return (np.corrcoef(x, y)[0, 1],
            np.corrcoef(np.argsort(np.argsort(x)), np.argsort(np.argsort(y)))[0, 1],
            np.polyfit(x, y, 1)[0], np.polyfit(x, y, 1)[1],
            np.sqrt(np.mean((x - y) ** 2)), np.abs(x - y).max())
inh = np.array([rand(f"{D}/cs_shuf/m2.ss_backsolve.ebv.txt")[i] for i in ids])
lines = [f"# ssGBLUP external benchmark vs HIBLUP (fixed VCs Va={VA} Ve={VE})"]
lines.append(f"LAMBDA_SOLVED = {VE/VA:.8f} (HIBLUP --vc-prior given as {VA},{VE/0.98:.8f}:"
             " HIBLUP applies an internal ~0.98 factor to --vc-prior, calibrated so both solve the same lambda)")
lines.append("EBV overall: pearson=%.8f spearman=%.8f slope=%.6f intercept=%.6f rmse=%.3e maxabs=%.3e" % rep(h, a))
for nm, msk in [("genotyped+phenotyped", (G1 == 1) & (P1 == 1)), ("genotyped,no pheno", (G1 == 1) & (P1 == 0)),
                ("ungenotyped+phenotyped", (G1 == 0) & (P1 == 1)), ("ungenotyped,no pheno", (G1 == 0) & (P1 == 0))]:
    x, y = h[msk], a[msk]
    if len(x):
        lines.append("EBV %-22s n=%4d pearson=%.8f slope=%.6f intercept=%.6f rmse=%.3e maxabs=%.3e"
                     % ((nm, len(x)) + rep(x, y)[0:1] + rep(x, y)[2:]))
# shuffle invariance
p2, s2, sl2, ic2, rm2, ma2 = rep(a, inh)
lines.append("Shuffle invariance (CS, reordered pedigree): corr=%.12f maxabs=%.2e rmse=%.2e"
             % (np.corrcoef(a, inh)[0, 1], np.abs(a - inh).max(), np.sqrt(np.mean((a - inh) ** 2))))
# SNP effects (HIBLUP needs a separate --snp-effect run; run it here if missing)
import os
if not os.path.exists(f"{D}/hi_snp_m4.snpeff"):
    os.system("/disk228/huangjie/software/miniconda3/envs/benchmark/bin/hiblup --single-trait"
              f" --pedigree {D}/sim.ped --bfile {D}/sim --pheno {D}/sim_m4.phe --pheno-pos 5"
              f" --snp-effect --vc-prior {VA},{VE/0.98:.8f} --threads 8 --out {D}/hi_snp_m4")
snpe = [l.split() for l in open(f"{D}/hi_snp_m4.snpeff")][1:]
add = np.array([float(x[4]) for x in snpe]); sn = np.array([x[0] for x in snpe])
ea = {}
for l in open(f"{D}/cs_m4/m2.ss_backsolve.snp.txt"):
    t = l.split()
    if t[0] != "CHR": ea[t[1]] = float(t[6])
ea = np.array([ea[s] for s in sn])
coef = np.polyfit(ea, add, 1)
res = add - (coef[0] * ea + coef[1])
lines.append("SNP effect add_a1 vs EFFECT_ALLELE: N=%d corr=%.8f OLS slope=%.6f intercept=%.1e rmse=%.3e"
             % (len(sn), np.corrcoef(ea, add)[0, 1], coef[0], coef[1], np.sqrt(np.mean(res ** 2))))
lines.append("NOTE: HIBLUP applies a data-dependent internal factor to its SNP-effect scale;"
             " corr/definition alignment is reported, exact scale parity is not claimable.")
open(f"{D}/benchmark_summary.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
EOF
echo
echo "Summary written to $OUT/benchmark_summary.txt"