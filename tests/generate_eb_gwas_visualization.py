import csv
import math
import sys
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
except ImportError:
    print("Error: matplotlib not installed", file=sys.stderr)
    sys.exit(1)

ROOT = Path(__file__).parent.parent
RESULT_DIR = ROOT / "tests" / "results"
OUTPUT_DIR = RESULT_DIR / "eb_gwas_visualization"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def read_tsv(path):
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def read_kv(path):
    out = {}
    if not path.exists():
        return out
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t", 1)
        if len(parts) == 2:
            out[parts[0]] = parts[1]
    return out

def fval(row, *names):
    lower = {k.lower(): k for k in row}
    for name in names:
        key = lower.get(name.lower())
        if key and row.get(key, "") not in ("", "NA", "nan", "-nan"):
            try:
                return float(row[key])
            except ValueError:
                pass
    return math.nan

print("[EB GWAS Visualization] Loading performance data...")
perf_rows = read_tsv(RESULT_DIR / "eb_performance_benchmark" / "eb_performance.summary.tsv")
print(f"  Loaded {len(perf_rows)} performance records")

print("[EB GWAS Visualization] Loading validation data...")
val_rows = read_tsv(RESULT_DIR / "eb_gwas_validation" / "eb_gwas_validation.summary.tsv")
print(f"  Loaded {len(val_rows)} validation records")

print("[EB GWAS Visualization] Loading calibration data...")
calib_rows = read_tsv(RESULT_DIR / "eb_calibration_benchmark" / "eb_calibration.aggregate.tsv")
print(f"  Loaded {len(calib_rows)} calibration records")

print("[EB GWAS Visualization] Loading comparison data...")
comp_data = read_kv(RESULT_DIR / "eb_comparison_report" / "standard_vs_eb.summary.tsv")
print(f"  Loaded comparison data: {list(comp_data.keys())}")

print("[EB GWAS Visualization] Loading multivariate benchmark data...")
mv_rows = read_tsv(RESULT_DIR / "multivariate_system_benchmark_n240_m1200_t2" / "method_metrics.tsv")
print(f"  Loaded {len(mv_rows)} multivariate benchmark records")

plt.style.use("seaborn-v0_8-whitegrid")
colors = ["#3b82f6", "#10b981", "#f59e0b", "#ef4444", "#8b5cf6", "#ec4899"]

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

standard_row = next((r for r in perf_rows if r["case"] == "standard"), None)
eb_rows = [r for r in perf_rows if r["case"].startswith("eb_t")]

if standard_row and eb_rows:
    std_time = fval(standard_row, "elapsed_sec")
    std_rss = fval(standard_row, "max_rss_kb") / 1024
    
    eb_times = [fval(r, "elapsed_sec") for r in eb_rows]
    eb_rss = [fval(r, "max_rss_kb") / 1024 for r in eb_rows]
    eb_threads = [int(r["threads"]) for r in eb_rows]
    
    labels = ["Standard LMM"] + [f"EB-LMM (t={t})" for t in eb_threads]
    times = [std_time] + eb_times
    rss_vals = [std_rss] + eb_rss
    
    x = range(len(labels))
    bars1 = axes[0].bar([i - 0.18 for i in x], times, width=0.35, color="#3b82f6", label="Runtime (sec)")
    bars2 = axes[1].bar([i - 0.18 for i in x], rss_vals, width=0.35, color="#10b981", label="Memory (MB)")
    
    for bar, val in zip(bars1, times):
        axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.02, f"{val:.3f}", ha="center", va="bottom", fontsize=8)
    
    for bar, val in zip(bars2, rss_vals):
        axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.02, f"{val:.1f}", ha="center", va="bottom", fontsize=8)
    
    axes[0].set_title("EB-LMM vs Standard LMM Runtime", fontsize=12)
    axes[0].set_ylabel("Elapsed Time (seconds)", fontsize=10)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=20, ha="right", fontsize=9)
    axes[0].legend()
    
    axes[1].set_title("EB-LMM vs Standard LMM Memory Usage", fontsize=12)
    axes[1].set_ylabel("Max RSS (MB)", fontsize=10)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=20, ha="right", fontsize=9)
    axes[1].legend()

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "eb_vs_standard_performance.png", dpi=180, bbox_inches="tight")
plt.savefig(OUTPUT_DIR / "eb_vs_standard_performance.pdf", bbox_inches="tight")
plt.close()
print("  Saved: eb_vs_standard_performance.png/pdf")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

null_rows = [r for r in val_rows if r["scenario"] == "null"]
power_rows = [r for r in val_rows if r["scenario"] == "power"]

models = [r["model"] for r in null_rows]
null_effective = [fval(r, "effective_tests") for r in null_rows]
null_bonf = [int(r["bonf_hits"]) for r in null_rows]
null_eb = [int(r["eb_mbonf_hits"]) for r in null_rows]

power_effective = [fval(r, "effective_tests") for r in power_rows]
power_causal = [int(r["causal_rank"]) for r in power_rows]

x = range(len(models))

bars_ef_null = axes[0].bar([i - 0.3 for i in x], null_effective, width=0.28, color="#8b5cf6", label="Effective tests")
bars_bonf = axes[0].bar([i for i in x], null_bonf, width=0.28, color="#ef4444", label="Bonf hits")
bars_eb = axes[0].bar([i + 0.3 for i in x], null_eb, width=0.28, color="#10b981", label="EB-mBonf hits")

for bar, val in zip(bars_ef_null, null_effective):
    axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.02, f"{val:.1f}", ha="center", va="bottom", fontsize=8)

axes[0].set_title("EB GWAS Validation - Null Scenarios", fontsize=12)
axes[0].set_ylabel("Count", fontsize=10)
axes[0].set_xticks(x)
axes[0].set_xticklabels(models, fontsize=10)
axes[0].legend()

bars_ef_power = axes[1].bar([i - 0.15 for i in x], power_effective, width=0.3, color="#8b5cf6", label="Effective tests")
bars_causal = axes[1].bar([i + 0.15 for i in x], power_causal, width=0.3, color="#f59e0b", label="Causal rank")

for bar, val in zip(bars_ef_power, power_effective):
    axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.02, f"{val:.1f}", ha="center", va="bottom", fontsize=8)

for bar, val in zip(bars_causal, power_causal):
    axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.02, str(val), ha="center", va="bottom", fontsize=8)

axes[1].set_title("EB GWAS Validation - Power Scenarios", fontsize=12)
axes[1].set_ylabel("Value", fontsize=10)
axes[1].set_xticks(x)
axes[1].set_xticklabels(models, fontsize=10)
axes[1].legend()

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "eb_validation_matrix.png", dpi=180, bbox_inches="tight")
plt.savefig(OUTPUT_DIR / "eb_validation_matrix.pdf", bbox_inches="tight")
plt.close()
print("  Saved: eb_validation_matrix.png/pdf")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

if calib_rows:
    scenarios = [r["scenario"] for r in calib_rows]
    reps = [int(r["reps"]) for r in calib_rows]
    bonf_hits = [int(r["standard_bonf_reps_with_hits"]) for r in calib_rows]
    eb_bonf_hits = [int(r["eb_bonf_reps_with_hits"]) for r in calib_rows]
    eb_mbonf_hits = [int(r["eb_mbonf_reps_with_hits"]) for r in calib_rows]
    eff_tests = [fval(r, "mean_effective_tests") for r in calib_rows]
    causal_top10 = [int(r["power_causal_top10"]) for r in calib_rows]
    
    x = range(len(scenarios))
    
    bars_bonf = axes[0].bar([i - 0.3 for i in x], bonf_hits, width=0.28, color="#ef4444", label="Standard Bonf")
    bars_eb_bonf = axes[0].bar([i for i in x], eb_bonf_hits, width=0.28, color="#3b82f6", label="EB Bonf")
    bars_eb_mbonf = axes[0].bar([i + 0.3 for i in x], eb_mbonf_hits, width=0.28, color="#10b981", label="EB-mBonf")
    
    for bar, val, max_val in zip(bars_bonf, bonf_hits, reps):
        axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.02, f"{val}/{max_val}", ha="center", va="bottom", fontsize=8)
    for bar, val, max_val in zip(bars_eb_bonf, eb_bonf_hits, reps):
        axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.02, f"{val}/{max_val}", ha="center", va="bottom", fontsize=8)
    for bar, val, max_val in zip(bars_eb_mbonf, eb_mbonf_hits, reps):
        axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.02, f"{val}/{max_val}", ha="center", va="bottom", fontsize=8)
    
    axes[0].set_title("EB Calibration - Type I Error Control", fontsize=12)
    axes[0].set_ylabel("Hit reps / Total reps", fontsize=10)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(scenarios, fontsize=10)
    axes[0].legend()
    
    bars_eff = axes[1].bar([i - 0.2 for i in x], eff_tests, width=0.4, color="#8b5cf6", label="Mean effective tests")
    bars_causal = axes[1].bar([i + 0.2 for i in x], causal_top10, width=0.4, color="#f59e0b", label="Causal top10 reps")
    
    for bar, val in zip(bars_eff, eff_tests):
        axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.02, f"{val:.1f}", ha="center", va="bottom", fontsize=8)
    
    for bar, val, max_val in zip(bars_causal, causal_top10, reps):
        axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.02, f"{val}/{max_val}", ha="center", va="bottom", fontsize=8)
    
    axes[1].set_title("EB Calibration - Power & Effective Tests", fontsize=12)
    axes[1].set_ylabel("Value", fontsize=10)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(scenarios, fontsize=10)
    axes[1].legend()

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "eb_calibration_benchmark.png", dpi=180, bbox_inches="tight")
plt.savefig(OUTPUT_DIR / "eb_calibration_benchmark.pdf", bbox_inches="tight")
plt.close()
print("  Saved: eb_calibration_benchmark.png/pdf")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

if mv_rows:
    methods = [r["method"] for r in mv_rows]
    times = [fval(r, "seconds") for r in mv_rows]
    rss = [fval(r, "max_rss_mb") for r in mv_rows]
    gemma_corr = [fval(r, "log10p_corr_vs_gemma") for r in mv_rows]
    emmv_corr = [fval(r, "log10p_corr_vs_emmvgwas") for r in mv_rows]
    
    x = range(len(methods))
    
    bars_time = axes[0].bar([i - 0.18 for i in x], times, width=0.35, color="#3b82f6", label="Runtime (sec)")
    bars_rss = axes[0].bar([i + 0.18 for i in x], rss, width=0.35, color="#10b981", label="Memory (MB)")
    
    for bar, val in zip(bars_time, times):
        axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.02, f"{val:.1f}", ha="center", va="bottom", fontsize=8)
    
    for bar, val in zip(bars_rss, rss):
        axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.02, f"{val:.1f}", ha="center", va="bottom", fontsize=8)
    
    axes[0].set_yscale("log")
    axes[0].set_title("Multivariate GWAS Runtime & Memory", fontsize=12)
    axes[0].set_ylabel("Seconds / MB (log scale)", fontsize=10)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(methods, rotation=25, ha="right", fontsize=8)
    axes[0].legend()
    
    bars_gemma = axes[1].bar([i - 0.18 for i in x], gemma_corr, width=0.35, color="#8b5cf6", label="Corr vs GEMMA")
    bars_emmvgwas = axes[1].bar([i + 0.18 for i in x], emmv_corr, width=0.35, color="#ec4899", label="Corr vs EMmvGWAS")
    
    for bar, val in zip(bars_gemma, gemma_corr):
        if math.isfinite(val):
            axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.02, f"{val:.4f}", ha="center", va="bottom", fontsize=8)
    
    for bar, val in zip(bars_emmvgwas, emmv_corr):
        if math.isfinite(val):
            axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.02, f"{val:.4f}", ha="center", va="bottom", fontsize=8)
    
    axes[1].set_ylim(0.95, 1.01)
    axes[1].set_title("Multivariate GWAS Accuracy (vs GEMMA/EMmvGWAS)", fontsize=12)
    axes[1].set_ylabel("-log10(P) Correlation", fontsize=10)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(methods, rotation=25, ha="right", fontsize=8)
    axes[1].legend()

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "mv_system_benchmark.png", dpi=180, bbox_inches="tight")
plt.savefig(OUTPUT_DIR / "mv_system_benchmark.pdf", bbox_inches="tight")
plt.close()
print("  Saved: mv_system_benchmark.png/pdf")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

if comp_data:
    markers = int(comp_data.get("markers", 0))
    top10_overlap = int(comp_data.get("top10_overlap", 0))
    mean_eb_confidence = float(comp_data.get("mean_eb_confidence", 0))
    median_lambda = float(comp_data.get("median_marker_lambda", 0))
    
    labels = ["Standard LMM", "EB-LMM"]
    top10_data = [10, top10_overlap]
    mean_conf_data = [1.0, mean_eb_confidence]
    
    axes[0].bar(labels, top10_data, color=["#3b82f6", "#10b981"])
    axes[0].set_title("Standard vs EB-LMM Top10 Overlap", fontsize=12)
    axes[0].set_ylabel("Top10 Overlap", fontsize=10)
    axes[0].set_ylim(0, 12)
    for i, v in enumerate(top10_data):
        axes[0].text(i, v + 0.3, str(v), ha="center", va="bottom", fontsize=10)
    
    axes[1].bar(labels, mean_conf_data, color=["#3b82f6", "#10b981"])
    axes[1].set_title("Mean EB Confidence", fontsize=12)
    axes[1].set_ylabel("Confidence Score", fontsize=10)
    axes[1].set_ylim(0, 1.1)
    for i, v in enumerate(mean_conf_data):
        axes[1].text(i, v + 0.03, f"{v:.4f}", ha="center", va="bottom", fontsize=10)
    
    info_text = f"Markers: {markers}\nMedian Lambda: {median_lambda:.4f}"
    fig.text(0.01, 0.01, info_text, fontsize=9, bbox=dict(facecolor="white", alpha=0.8))

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "standard_vs_eb_consistency.png", dpi=180, bbox_inches="tight")
plt.savefig(OUTPUT_DIR / "standard_vs_eb_consistency.pdf", bbox_inches="tight")
plt.close()
print("  Saved: standard_vs_eb_consistency.png/pdf")

print(f"\n[EB GWAS Visualization] All figures saved to: {OUTPUT_DIR}")
print("\nGenerated files:")
for f in sorted(OUTPUT_DIR.glob("*.png")):
    print(f"  {f.name}")