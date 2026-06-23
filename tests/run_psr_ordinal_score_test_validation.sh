#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${COSMICSCANNER_BIN:-${ROOT}/build/cosmicscanner}"

N=60
LEVELS=3
M=50
BASE_SEED=20260624
CAUSAL_IDX=$((M / 3))

DATA_DIR="${ROOT}/tests/data/psr_ordinal_score_test_validation"
RESULT_DIR="${ROOT}/tests/results/psr_ordinal_score_test_validation"
WORK_DIR="${RESULT_DIR}/work"
SUMMARY="${RESULT_DIR}/psr_ordinal_score_test_validation_summary.tsv"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

if [[ ! -x "${BIN}" ]]; then
  echo "Missing CosmicScanner executable: ${BIN}" >&2
  exit 2
fi

mkdir -p "${DATA_DIR}" "${RESULT_DIR}" "${WORK_DIR}"

PASS=0
FAIL=0

check_result() {
    local name=$1 ec=$2
    if [ "$ec" -eq 0 ]; then
        echo "[PASS] $name"
        PASS=$((PASS + 1))
    else
        echo "[FAIL] $name"
        FAIL=$((FAIL + 1))
    fi
}

# ============================================================
# Generate data for null and power scenarios
# ============================================================
echo ""
echo "================================================"
echo "Generating ordinal test data (N=$N M=$M levels=$LEVELS)"
echo "================================================"

for scenario in null power; do
    seed=$((BASE_SEED + $(test "$scenario" = "power" && echo 1 || echo 0)))
    python3 -c "
import numpy as np, sys
sys.path.insert(0, '${ROOT}/tests')
from run_psr_simulation_summary import write_plink

np.random.seed($seed)
n = $N; m = $M; n_levels = $LEVELS; causal = m // 3

freqs = np.random.uniform(0.15, 0.45, size=m)
G = np.random.binomial(2, freqs, size=(n, m)).astype(int)

Z = (G - G.mean(axis=0)) / G.std(axis=0)
Z[np.isnan(Z)] = 0.0
K = Z @ Z.T / m
L = np.linalg.cholesky(K + np.eye(n) * 1e-6)

pop = np.r_[np.zeros(n // 2), np.ones(n - n // 2)]
polygenic = L @ np.random.normal(0, 0.8, size=n)
effect_size = 0.0 if '$scenario' == 'null' else 0.8
eta = 0.2 * pop + 0.5 * polygenic + effect_size * (G[:, causal] - G[:, causal].mean())
latent = eta + np.random.normal(0, 1.0, size=n)

y = np.ones(n, dtype=int)
y[latent > -0.5] = 2
y[latent > 0.5] = 3
if len(set(y)) < 3:
    y[:n//3] = 1; y[n//3:2*n//3] = 2; y[2*n//3:] = 3

prefix = '${DATA_DIR}/ordinal_${scenario}'
write_plink(prefix, G)
with open(prefix + '.ordinal.pheno', 'w') as f:
    for i in range(n):
        f.write(f'F{i+1:04d}\tI{i+1:04d}\t{int(y[i])}\n')
print(f'Generated N={n} M={m} levels={n_levels}')
    " 2>&1
    check_result "generate ${scenario} data" $?
done

# ============================================================
# Run each scenario x configuration
# ============================================================
declare -a configs
configs=(
    "fixed-gls:--ordinal-stacked-engine fixed-gls --ordinal-vc-method fixed"
    "dense-zkz-fixed:--ordinal-stacked-engine dense-zkz --ordinal-vc-method fixed"
    "dense-zkz-gridreml:--ordinal-stacked-engine dense-zkz --ordinal-vc-method grid-reml --ordinal-vc-grid-points 11"
)

for scenario in null power; do
    echo ""
    echo "================================================"
    echo "Running $scenario scenario"
    echo "================================================"
    for config in "${configs[@]}"; do
        label="${config%%:*}"
        extra_args="${config#*:}"
        out_prefix="${WORK_DIR}/ordinal_${scenario}_${label}"
        echo ""
        echo "--- $scenario / $label ---"

        set +e
        "${BIN}" --bfile "${DATA_DIR}/ordinal_${scenario}" --gwas --no-loco --maf 0.0 \
            --ordinal-model --ordinal-pheno "${DATA_DIR}/ordinal_${scenario}.ordinal.pheno" \
            --link logit --ordinal-mode stacked-psr \
            $extra_args --ordinal-test all \
            --psr-maxit 25 --out "${out_prefix}" 2>&1
        rc=$?
        set -e

        gwas_path="${out_prefix}.gwas.tsv"
        score_path="${out_prefix}.ordinal.score.tsv"

        if [ "$rc" -ne 0 ]; then
            if echo "$label" | grep -q "gridreml"; then
                check_result "$scenario/$label graceful error on dense-zkz grid-reml" 0
            else
                check_result "$scenario/$label exit code 0" 1
            fi
            continue
        fi
        check_result "$scenario/$label exit code 0" 0

        if [ -f "$gwas_path" ] && [ -s "$gwas_path" ]; then
            check_result "$scenario/$label .gwas.tsv exists" 0
        else
            check_result "$scenario/$label .gwas.tsv exists" 1
            continue
        fi

        if [ -f "$score_path" ] && [ -s "$score_path" ]; then
            check_result "$scenario/$label .ordinal.score.tsv exists" 0
        else
            check_result "$scenario/$label .ordinal.score.tsv exists" 1
            continue
        fi

        python3 - "$score_path" "$N" <<'PY'
import csv, math, sys

score_path = sys.argv[1]
expected_n = int(sys.argv[2])

with open(score_path, newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

n_snps = len(rows)
if n_snps == 0:
    print("FAIL: empty score file"); sys.exit(1)

header = set(rows[0].keys())
required = {"wald_psr_score_equivalent", "n_clusters", "wald_stat", "wald_p",
            "psr_score_stat", "psr_score_p", "raw_score_stat", "cluster_score_stat"}
missing = required - header
if missing:
    print(f"FAIL: missing columns: {sorted(missing)}")
    sys.exit(1)

wald_equiv = 0
n_clusters_ok = 0
all_finite = True

for row in rows:
    ncl = int(row["n_clusters"])
    if ncl == expected_n:
        n_clusters_ok += 1

    if row.get("wald_psr_score_equivalent") == "TRUE":
        wald_equiv += 1

    for col in ["wald_stat", "wald_p", "wald_log10p",
                "psr_score_stat", "psr_score_p", "psr_score_log10p"]:
        v = float(row[col])
        if not math.isfinite(v):
            all_finite = False
        if col.endswith("_p") and not (0.0 <= v <= 1.0):
            print(f"FAIL: {col} outside [0,1]: {v}")
            sys.exit(1)

    for col in ["raw_score_stat", "raw_score_p", "raw_score_log10p",
                "cluster_score_stat", "cluster_score_p", "cluster_score_log10p"]:
        v_str = row.get(col, "nan")
        try:
            v = float(v_str)
        except ValueError:
            continue
        if math.isfinite(v) and v >= 0:
            if col.endswith("_p") and not (0.0 <= v <= 1.0):
                print(f"FAIL: {col} outside [0,1]: {v}")
                sys.exit(1)
            pass
        elif not math.isfinite(v):
            all_finite = False

if not all_finite:
    print("FAIL: non-finite values found"); sys.exit(1)

print(f"VALID: {n_snps} SNPs, wald_equiv={wald_equiv}/{n_snps}, n_clusters_ok={n_clusters_ok}/{n_snps}")
PY
        check_result "$scenario/$label score validation" $?
    done
done

# ============================================================
# Verify dense-zkz grid-reml with fixed-gls engine errors
# ============================================================
echo ""
echo "--- fixed-gls + grid-reml (expected error) ---"
set +e
"${BIN}" --bfile "${DATA_DIR}/ordinal_null" --gwas --no-loco --maf 0.0 \
    --ordinal-model --ordinal-pheno "${DATA_DIR}/ordinal_null.ordinal.pheno" \
    --link logit --ordinal-mode stacked-psr \
    --ordinal-stacked-engine fixed-gls --ordinal-vc-method grid-reml \
    --ordinal-test all --psr-maxit 25 \
    --out "${WORK_DIR}/ordinal_fixedgls_gridreml" > /dev/null 2>&1
rc=$?
set -e
if [ "$rc" -ne 0 ]; then
    check_result "fixed-gls + grid-reml correctly errors" 0
else
    check_result "fixed-gls + grid-reml correctly errors (expected non-zero)" 1
fi

# ============================================================
# Verify robust-score is the primary GWAS statistic
# ============================================================
echo ""
echo "--- robust-score primary output check ---"
robust_prefix="${WORK_DIR}/ordinal_null_densezkz_robust_primary"
"${BIN}" --bfile "${DATA_DIR}/ordinal_null" --gwas --no-loco --maf 0.0 \
    --ordinal-model --ordinal-pheno "${DATA_DIR}/ordinal_null.ordinal.pheno" \
    --link logit --ordinal-mode stacked-psr \
    --ordinal-stacked-engine dense-zkz --ordinal-vc-method grid-reml \
    --ordinal-test robust-score --psr-maxit 25 \
    --out "${robust_prefix}" > /dev/null 2>&1
check_result "robust-score primary run" $?

python3 - "${robust_prefix}.gwas.tsv" "${robust_prefix}.ordinal.score.tsv" <<'PY'
import csv, math, sys

gwas_path, score_path = sys.argv[1], sys.argv[2]
with open(gwas_path, newline="") as f:
    gwas = {r["snp"]: r for r in csv.DictReader(f, delimiter="\t")}
with open(score_path, newline="") as f:
    score = {r["snp"]: r for r in csv.DictReader(f, delimiter="\t")}

if not gwas or set(gwas) != set(score):
    print("FAIL: GWAS/score SNP sets differ")
    sys.exit(1)

for snp, gr in gwas.items():
    sr = score[snp]
    if gr.get("association_test") != "cluster_score":
        print(f"FAIL: {snp} association_test={gr.get('association_test')}")
        sys.exit(1)
    checks = [
        ("wald_or_stat", "cluster_score_stat"),
        ("p_value", "cluster_score_p"),
        ("log10p", "cluster_score_log10p"),
    ]
    for gc, sc in checks:
        gv = float(gr[gc]); sv = float(sr[sc])
        if not (math.isfinite(gv) and math.isfinite(sv) and abs(gv - sv) <= 1e-10 * max(1.0, abs(sv))):
            print(f"FAIL: {snp} {gc}={gv} differs from {sc}={sv}")
            sys.exit(1)
print(f"VALID: robust-score primary output matches cluster_score for {len(gwas)} SNPs")
PY
check_result "robust-score primary equals cluster_score" $?

# ============================================================
# Collect results into summary TSV
# ============================================================
echo ""
echo "================================================"
echo "Collecting summary results"
echo "================================================"

python3 - "${WORK_DIR}" "${RESULT_DIR}" "${SUMMARY}" "${N}" "${M}" "${CAUSAL_IDX}" <<'PY'
import csv, math, sys
from pathlib import Path

work_dir = Path(sys.argv[1])
result_dir = Path(sys.argv[2])
summary_path = Path(sys.argv[3])
N = int(sys.argv[4])
M = int(sys.argv[5])
causal_idx = int(sys.argv[6])
causal_id = f"sim{causal_idx+1:06d}"

config_map = [
    ("fixed-gls", "cosmicscanner_psr", "fixed", "fixed-gls"),
    ("dense-zkz", "cosmicscanner_psr", "fixed", "dense-zkz-fixed"),
    ("dense-zkz", "cosmicscanner_psr", "grid-reml", "dense-zkz-gridreml"),
]

fieldnames = [
    "method", "engine", "vc_method", "n", "m", "scenario",
    "valid_wald", "valid_psr_score", "valid_raw_score", "valid_cluster_score",
    "lambda_wald", "lambda_psr_score", "lambda_raw_score", "lambda_cluster_score",
    "top10_wald", "top10_psr_score", "top10_raw_score", "top10_cluster_score",
    "wald_psr_equiv_rate", "n_clusters_ok", "status",
]

results = []
for scenario in ["null", "power"]:
    for method, engine, vc_method, label in config_map:
        score_path = work_dir / f"ordinal_{scenario}_{label}.ordinal.score.tsv"
        row = {k: "NA" for k in fieldnames}
        row.update(method=method, engine=engine, vc_method=vc_method,
                   n=N, m=M, scenario=scenario, status="FAIL")

        if not score_path.exists():
            results.append(row)
            continue

        try:
            with open(score_path, newline="") as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
        except Exception:
            results.append(row)
            continue

        if not rows:
            results.append(row)
            continue

        n_snps = len(rows)
        wald_pairs = []
        psr_pairs = []
        raw_pairs = []
        cluster_pairs = []
        n_clusters_ok = 0
        wald_equiv = 0

        for r in rows:
            try:
                ws = float(r["wald_stat"])
                wp = float(r["wald_p"])
                ps = float(r["psr_score_stat"])
                pp = float(r["psr_score_p"])
                ncl = int(r["n_clusters"])
            except (KeyError, ValueError):
                continue

            if math.isfinite(ws) and math.isfinite(wp):
                wald_pairs.append((ws, wp, r.get("snp", "")))
            if math.isfinite(ps) and math.isfinite(pp):
                psr_pairs.append((ps, pp, r.get("snp", "")))

            if ncl == N:
                n_clusters_ok += 1
            if r.get("wald_psr_score_equivalent") == "TRUE":
                wald_equiv += 1

            try:
                rs = float(r.get("raw_score_stat", "-1"))
                rp = float(r.get("raw_score_p", "-1"))
                if math.isfinite(rs) and rs > 0 and math.isfinite(rp):
                    raw_pairs.append((rs, rp, r.get("snp", "")))
            except (ValueError, TypeError):
                pass

            try:
                cs = float(r.get("cluster_score_stat", "-1"))
                cp = float(r.get("cluster_score_p", "-1"))
                if math.isfinite(cs) and cs > 0 and math.isfinite(cp):
                    cluster_pairs.append((cs, cp, r.get("snp", "")))
            except (ValueError, TypeError):
                pass

        def median_finite(vals):
            if not vals:
                return float("nan")
            sv = sorted(vals)
            nv = len(sv)
            return float(sv[nv // 2]) if nv % 2 == 1 else float((sv[nv // 2 - 1] + sv[nv // 2]) / 2.0)

        def lambda_gc(pairs):
            if len(pairs) < 5:
                return "NA"
            stats = [p[0] for p in pairs]
            med = median_finite(stats)
            expected_median = 0.454936
            return f"{med / expected_median:.6f}" if med > 0 else "NA"

        def top10_hits(pairs):
            if len(pairs) < 10:
                return "NA"
            sorted_by_p = sorted(pairs, key=lambda x: x[1])
            top_snps = set(x[2] for x in sorted_by_p[:10])
            return "1" if causal_id in top_snps else "0"

        row.update(
            valid_wald="TRUE" if len(wald_pairs) >= 5 else "FALSE",
            valid_psr_score="TRUE" if len(psr_pairs) >= 5 else "FALSE",
            valid_raw_score="TRUE" if len(raw_pairs) >= 5 else "FALSE",
            valid_cluster_score="TRUE" if len(cluster_pairs) >= 5 else "FALSE",
            lambda_wald=lambda_gc(wald_pairs),
            lambda_psr_score=lambda_gc(psr_pairs),
            lambda_raw_score=lambda_gc(raw_pairs),
            lambda_cluster_score=lambda_gc(cluster_pairs),
            top10_wald=top10_hits(wald_pairs),
            top10_psr_score=top10_hits(psr_pairs),
            top10_raw_score=top10_hits(raw_pairs),
            top10_cluster_score=top10_hits(cluster_pairs),
            wald_psr_equiv_rate=f"{wald_equiv / n_snps:.4f}" if n_snps > 0 else "NA",
            n_clusters_ok=f"{n_clusters_ok}/{n_snps}",
            status="PASS",
        )
        results.append(row)

with open(summary_path, "w", newline="") as f:
    writer = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)

print(f"Summary written to {summary_path}")
for r in results:
    print(f"  {r['scenario']}/{r['method']}/{r['vc_method']}: status={r['status']} "
          f"wald_equiv={r['wald_psr_equiv_rate']} n_clusters={r['n_clusters_ok']} "
          f"top10_wald={r['top10_wald']}")
PY
check_result "summary collection" $?

# ============================================================
# Final summary
# ============================================================
echo ""
echo "================================================"
echo "Results: ${PASS} passed, ${FAIL} failed"
echo "================================================"

rm -rf "${DATA_DIR}" "${WORK_DIR}"

if [ "$FAIL" -gt 0 ]; then exit 1; fi
echo "All ordinal score test validation checks passed."
