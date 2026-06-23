#!/usr/bin/env python3
"""Compare two runs of the PSR-GWAS ordinal benchmark.

Reads v05_ordinal_method_summary.tsv and v05_ordinal_benchmark_runs.tsv
from two directories and produces a markdown comparison report with
per-method metric deltas, per-scenario breakdown, and summary verdict.
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

DEFAULT_BASE = (
    Path(__file__).resolve().parents[1]
    / "tests" / "results" / "psr_v05_ordinal_benchmark"
)


def find_run_dirs(base):
    if not base.is_dir():
        return []
    return sorted(
        [d for d in base.iterdir() if d.is_dir()],
        key=lambda d: d.stat().st_mtime,
        reverse=True,
    )


def read_tsv(path):
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def safe_float(v):
    if v in (None, "", "NA", "nan", "NaN"):
        return float("nan")
    try:
        return float(v)
    except (ValueError, TypeError):
        return float("nan")


def safe_int(v):
    if v in (None, "", "NA", "nan", "NaN"):
        return None
    try:
        return int(float(v))
    except (ValueError, TypeError):
        return None


def fmt_num(v, decimals=3):
    if isinstance(v, float) and np.isnan(v):
        return "NA"
    if decimals == 0:
        return f"{v:.0f}"
    return f"{v:.{decimals}f}"


def fmt_pct(v):
    if isinstance(v, float) and np.isnan(v):
        return "NA"
    return f"{v * 100:.1f}%"


def compute_delta(old, new):
    if np.isnan(old) or np.isnan(new):
        return float("nan"), float("nan")
    delta = new - old
    pct = delta / abs(old) * 100 if abs(old) > 1e-12 else float("nan")
    return delta, pct


def summarize_verdict(improved, regressed, unchanged):
    if improved > regressed and improved > 0:
        return "IMPROVED"
    if regressed > improved:
        return "REGRESSED"
    return "STABLE"


def main():
    parser = argparse.ArgumentParser(
        description="Compare two PSR-GWAS ordinal benchmark runs"
    )
    parser.add_argument("--old", default=None, help="Old benchmark directory")
    parser.add_argument("--new", default=None, help="New benchmark directory")
    parser.add_argument("--output", default=None,
                        help="Output markdown file (default: comparison_report.md in --new dir, also prints to stdout)")
    args = parser.parse_args()

    # Resolve default directories
    if args.old is None or args.new is None:
        dirs = find_run_dirs(DEFAULT_BASE)
        if args.new is None:
            args.new = str(dirs[0]) if dirs else str(DEFAULT_BASE)
        if args.old is None:
            args.old = str(dirs[1]) if len(dirs) > 1 else str(DEFAULT_BASE)

    old_dir = Path(args.old)
    new_dir = Path(args.new)

    old_summary = read_tsv(old_dir / "v05_ordinal_method_summary.tsv")
    new_summary = read_tsv(new_dir / "v05_ordinal_method_summary.tsv")
    old_runs = read_tsv(old_dir / "v05_ordinal_benchmark_runs.tsv")
    new_runs = read_tsv(new_dir / "v05_ordinal_benchmark_runs.tsv")

    if not old_summary and not new_summary:
        print("No v05_ordinal_method_summary.tsv found in either directory.",
              file=sys.stderr)
        sys.exit(1)

    old_map = {r["ordinal_method"]: r for r in old_summary}
    new_map = {r["ordinal_method"]: r for r in new_summary}
    all_methods = sorted(set(list(old_map.keys()) + list(new_map.keys())))

    lines = [
        "# PSR-GWAS v0.5 Ordinal Benchmark Comparison Report",
        "",
        f"- **Old directory**: `{old_dir}`",
        f"- **New directory**: `{new_dir}`",
        "",
        "## Per-Method Metric Comparison",
        "",
    ]

    # ---------- Per-method summary comparison ----------
    for method in all_methods:
        old = old_map.get(method, {})
        new = new_map.get(method, {})

        lines.append(f"### {method}")
        lines.append("")
        lines.append("| Metric | Old | New | Delta | %Change |")
        lines.append("| --- | --- | --- | --- | --- |")

        # N runs
        old_n = safe_int(old.get("n_runs"))
        new_n = safe_int(new.get("n_runs"))
        if old_n is not None and new_n is not None:
            d = new_n - old_n
            pct = d / old_n * 100 if old_n else float("nan")
            lines.append(f"| N runs | {old_n} | {new_n} | {d:+d} | {fmt_num(pct, 1)}% |")
        elif old_n is not None:
            lines.append(f"| N runs | {old_n} | NA | NA | NA |")
        elif new_n is not None:
            lines.append(f"| N runs | NA | {new_n} | NA | NA |")
        else:
            lines.append("| N runs | NA | NA | NA | NA |")

        # Runtime (mean ± sd)
        old_rt = safe_float(old.get("runtime_mean"))
        old_rtsd = safe_float(old.get("runtime_sd"))
        new_rt = safe_float(new.get("runtime_mean"))
        new_rtsd = safe_float(new.get("runtime_sd"))

        def rt_str(m, s):
            if np.isnan(m):
                return "NA"
            return f"{m:.3f}\u00b1{s:.3f}" if not np.isnan(s) else f"{m:.3f}"

        o_rt_str = rt_str(old_rt, old_rtsd)
        n_rt_str = rt_str(new_rt, new_rtsd)
        delta_rt, pct_rt = compute_delta(old_rt, new_rt)
        lines.append(
            f"| Runtime (mean\u00b1sd) | {o_rt_str} | {n_rt_str} | "
            f"{fmt_num(delta_rt, 3)}s | {fmt_num(pct_rt, 1)}% |"
        )

        # Scalar metrics
        for label, key, formatter in [
            ("Convergence Rate", "convergence_rate", fmt_pct),
            ("Lambda GC (mean)", "lambda_gc_mean", fmt_num),
            ("Type I Error @0.05", "type1_rate", fmt_pct),
            ("Top10 Hit Rate", "top10_hit_rate", fmt_pct),
            ("Sigma2_g (mean)", "sigma2_g_mean", fmt_num),
            ("Sigma2_g Boundary Rate", "sigma2_g_boundary_rate", fmt_pct),
        ]:
            old_v = safe_float(old.get(key))
            new_v = safe_float(new.get(key))
            if np.isnan(old_v) and np.isnan(new_v):
                continue
            delta_v, pct_v = compute_delta(old_v, new_v)
            delta_s = f"{delta_v:+.4f}" if not np.isnan(delta_v) else "NA"
            pct_s = f"{pct_v:+.2f}%" if not np.isnan(pct_v) else "NA"
            lines.append(
                f"| {label} | {formatter(old_v)} | {formatter(new_v)} | "
                f"{delta_s} | {pct_s} |"
            )
        lines.append("")

    # ---------- Per-scenario breakdown ----------
    lines.append("## Per-Scenario Breakdown")
    lines.append("")

    if old_runs or new_runs:
        scenarios = sorted(
            set(r.get("scenario", "") for r in old_runs) |
            set(r.get("scenario", "") for r in new_runs)
        )
        for sc in scenarios:
            old_sc = [r for r in old_runs if r.get("scenario") == sc]
            new_sc = [r for r in new_runs if r.get("scenario") == sc]
            if not old_sc and not new_sc:
                continue
            lines.append(f"### Scenario: `{sc}`")
            lines.append("")
            for method in all_methods:
                old_m = [r for r in old_sc if r.get("ordinal_method") == method]
                new_m = [r for r in new_sc if r.get("ordinal_method") == method]
                if not old_m and not new_m:
                    continue
                lines.append(f"**{method}**  ")

                # Null: lambda GC
                if sc.startswith("null"):
                    old_lams = [safe_float(r.get("lambda_gc")) for r in old_m
                                if not np.isnan(safe_float(r.get("lambda_gc")))]
                    new_lams = [safe_float(r.get("lambda_gc")) for r in new_m
                                if not np.isnan(safe_float(r.get("lambda_gc")))]
                    if old_lams and new_lams:
                        lines.append(
                            f"  - Lambda GC: old={np.mean(old_lams):.4f} "
                            f"(n={len(old_lams)}), new={np.mean(new_lams):.4f} "
                            f"(n={len(new_lams)})"
                        )
                    elif old_lams:
                        lines.append(f"  - Lambda GC: old={np.mean(old_lams):.4f} (n={len(old_lams)}), new=NA")
                    elif new_lams:
                        lines.append(f"  - Lambda GC: old=NA, new={np.mean(new_lams):.4f} (n={len(new_lams)})")

                # Power: top10 hit, runtime, sigma2_g
                if sc.startswith("power"):
                    old_t10 = [safe_int(r.get("top10_hit")) for r in old_m]
                    new_t10 = [safe_int(r.get("top10_hit")) for r in new_m]
                    old_t10 = [x for x in old_t10 if x is not None]
                    new_t10 = [x for x in new_t10 if x is not None]
                    if old_t10 and new_t10:
                        lines.append(
                            f"  - Top10 hit rate: old={sum(old_t10)/len(old_t10):.0%} "
                            f"({sum(old_t10)}/{len(old_t10)}), "
                            f"new={sum(new_t10)/len(new_t10):.0%} "
                            f"({sum(new_t10)}/{len(new_t10)})"
                        )
                    elif old_t10:
                        lines.append(
                            f"  - Top10 hit rate: old={sum(old_t10)/len(old_t10):.0%} "
                            f"({sum(old_t10)}/{len(old_t10)}), new=NA"
                        )
                    elif new_t10:
                        lines.append(
                            f"  - Top10 hit rate: old=NA, "
                            f"new={sum(new_t10)/len(new_t10):.0%} "
                            f"({sum(new_t10)}/{len(new_t10)})"
                        )

                    old_rt_vals = [safe_float(r.get("runtime_sec")) for r in old_m
                                   if not np.isnan(safe_float(r.get("runtime_sec")))]
                    new_rt_vals = [safe_float(r.get("runtime_sec")) for r in new_m
                                   if not np.isnan(safe_float(r.get("runtime_sec")))]
                    if old_rt_vals and new_rt_vals:
                        lines.append(
                            f"  - Runtime: old={np.mean(old_rt_vals):.3f}s, "
                            f"new={np.mean(new_rt_vals):.3f}s"
                        )

                old_sig = [safe_float(r.get("sigma2_g")) for r in old_m
                           if not np.isnan(safe_float(r.get("sigma2_g")))]
                new_sig = [safe_float(r.get("sigma2_g")) for r in new_m
                           if not np.isnan(safe_float(r.get("sigma2_g")))]
                if old_sig and new_sig:
                    lines.append(
                        f"  - Sigma2_g: old={np.mean(old_sig):.4f}, "
                        f"new={np.mean(new_sig):.4f}"
                    )

                old_bnd = [safe_float(r.get("sigma2_g_boundary")) for r in old_m
                           if not np.isnan(safe_float(r.get("sigma2_g_boundary")))]
                new_bnd = [safe_float(r.get("sigma2_g_boundary")) for r in new_m
                           if not np.isnan(safe_float(r.get("sigma2_g_boundary")))]
                if old_bnd and new_bnd:
                    lines.append(
                        f"  - Boundary rate: old={np.mean(old_bnd):.0%}, "
                        f"new={np.mean(new_bnd):.0%}"
                    )
                lines.append("")
    else:
        lines.append("No per-run data available for scenario breakdown.\n")

    # ---------- Summary verdict ----------
    lines.append("## Summary Verdict")
    lines.append("")
    lines.append("| Method | Improved | Regressed | Unchanged | Verdict |")
    lines.append("| --- | --- | --- | --- | --- |")

    for method in all_methods:
        old = old_map.get(method, {})
        new = new_map.get(method, {})

        metrics = []

        # convergence rate – higher is better, threshold 0.02
        old_cr = safe_float(old.get("convergence_rate"))
        new_cr = safe_float(new.get("convergence_rate"))
        metrics.append(("conv_rate", old_cr, new_cr, True, 0.02))

        # runtime – lower is better, threshold 0.5 s
        old_rt = safe_float(old.get("runtime_mean"))
        new_rt = safe_float(new.get("runtime_mean"))
        metrics.append(("runtime", old_rt, new_rt, False, 0.5))

        # lambda GC distance from 1 – lower distance is better, threshold 0.05
        old_lg = safe_float(old.get("lambda_gc_mean"))
        new_lg = safe_float(new.get("lambda_gc_mean"))
        if not np.isnan(old_lg) and not np.isnan(new_lg):
            metrics.append(("lambda_gc", abs(old_lg - 1.0), abs(new_lg - 1.0), False, 0.05))

        # top10 hit rate – higher is better, threshold 0.02
        old_t10 = safe_float(old.get("top10_hit_rate"))
        new_t10 = safe_float(new.get("top10_hit_rate"))
        metrics.append(("top10", old_t10, new_t10, True, 0.02))

        # boundary rate – lower is better, threshold 0.02
        old_br = safe_float(old.get("sigma2_g_boundary_rate"))
        new_br = safe_float(new.get("sigma2_g_boundary_rate"))
        metrics.append(("boundary", old_br, new_br, False, 0.02))

        improved = regressed = unchanged = 0
        for _, old_v, new_v, higher_better, threshold in metrics:
            if np.isnan(old_v) or np.isnan(new_v):
                unchanged += 1
                continue
            delta = new_v - old_v
            if abs(delta) < threshold:
                unchanged += 1
            elif (delta > 0 and higher_better) or (delta < 0 and not higher_better):
                improved += 1
            else:
                regressed += 1

        verdict = summarize_verdict(improved, regressed, unchanged)
        lines.append(f"| {method} | {improved} | {regressed} | {unchanged} | {verdict} |")

    lines.append("")
    lines.append("*Verdict based on: convergence rate, runtime, lambda GC distance from 1.0, top10 hit rate, sigma2_g boundary rate. Threshold-based: |delta| < threshold → unchanged.*")
    lines.append("")
    lines.append("---")
    lines.append(f"*Report auto-generated by `compare_psr_benchmarks.py`*")

    report = "\n".join(lines)

    # Output
    print(report)

    output_path = None
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = new_dir / "comparison_report.md"
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(report + "\n")
        print(f"\n[Compare] Report written to {output_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
