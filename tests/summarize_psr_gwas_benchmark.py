#!/usr/bin/env python3
import csv
import math
import sys
from pathlib import Path


def read_table(path):
    with open(path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fval(row, *names):
    lower = {k.lower(): k for k in row}
    for name in names:
        key = name if name in row else lower.get(name.lower())
        if key and row.get(key, "") not in ("", "NA", "nan", "-nan"):
            try:
                return float(row[key])
            except ValueError:
                return math.nan
    return math.nan


def sval(row, *names):
    lower = {k.lower(): k for k in row}
    for name in names:
        key = name if name in row else lower.get(name.lower())
        if key:
            return row[key]
    return ""


def corr(xs, ys):
    pairs = [(x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if len(pairs) < 3:
        return math.nan
    mx = sum(x for x, _ in pairs) / len(pairs)
    my = sum(y for _, y in pairs) / len(pairs)
    vx = sum((x - mx) ** 2 for x, _ in pairs)
    vy = sum((y - my) ** 2 for _, y in pairs)
    if vx <= 0.0 or vy <= 0.0:
        return math.nan
    return sum((x - mx) * (y - my) for x, y in pairs) / math.sqrt(vx * vy)


def read_timing(path):
    out = {}
    if not path.exists():
        return out
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                out[parts[0]] = parts[1]
    return out


def read_convergence(result_dir, trait):
    path = result_dir / f"cosmic_{trait}.psr.convergence.tsv"
    if not path.exists():
        return None
    rows = read_table(path)
    if not rows:
        return None
    final = rows[-1]
    total_clipped = 0
    warning_flags = []
    prev_quasi = None
    for row in rows:
        clipped = int(float(row.get("n_clipped", 0) or 0))
        total_clipped += clipped
        if clipped > 0:
            warning_flags.append("WARNING_PSR_CLIPPED")
        quasi = fval(row, "quasiLik")
        if prev_quasi is not None and math.isfinite(quasi) and math.isfinite(prev_quasi) and quasi < prev_quasi - 1e-8:
            warning_flags.append("WARNING_PSR_QUASILIK_DECREASED")
        prev_quasi = quasi
    converged = int(float(final.get("converged", 0) or 0))
    if not converged:
        warning_flags.append("WARNING_PSR_MAXIT_REACHED")
    sg = fval(final, "sigma_g2")
    if math.isfinite(sg) and sg <= 1e-5:
        warning_flags.append("WARNING_PSR_SIGMA_BOUNDARY")
    warning_flag = ",".join(sorted(set(warning_flags))) if warning_flags else "OK"
    status_class = "OK"
    if not converged:
        status_class = "NOT_CONVERGED"
    elif warning_flag != "OK":
        status_class = "CONVERGED_WITH_WARNING"
    return {
        "trait": trait,
        "n_iter": len(rows),
        "converged": converged,
        "final_sigma_g2": sg,
        "final_delta_sigma_g2": fval(final, "delta_sigma_g2"),
        "final_quasiLik": fval(final, "quasiLik"),
        "min_mu": min(fval(r, "min_mu") for r in rows),
        "max_mu": max(fval(r, "max_mu") for r in rows),
        "min_weight": min(fval(r, "min_weight") for r in rows),
        "max_weight": max(fval(r, "max_weight") for r in rows),
        "total_clipped": total_clipped,
        "warning_flag": warning_flag,
        "status_class": status_class,
    }


def read_weighted_equivalence(result_dir):
    path = result_dir / "weighted_lmm_equivalence.tsv"
    if not path.exists():
        return {"weighted_max_abs_diff": "NA", "weighted_all_pass": "NA", "weighted_tests": "NA"}
    rows = read_table(path)
    max_diff = 0.0
    all_pass = True
    tests = []
    for row in rows:
        tests.append(sval(row, "test"))
        all_pass = all_pass and sval(row, "pass") in ("1", "TRUE", "True", "true")
        for key in ("beta_abs_diff", "se_abs_diff", "wald_abs_diff", "p_abs_diff"):
            val = fval(row, key)
            if math.isfinite(val):
                max_diff = max(max_diff, abs(val))
    return {"weighted_max_abs_diff": max_diff, "weighted_all_pass": int(all_pass), "weighted_tests": ",".join(tests)}


def summarize_refinement(result_dir, trait):
    paths = list(result_dir.glob(f"*{trait}*.refined*.tsv")) + list(result_dir.glob(f"candidate_suite.{trait}.*.tsv"))
    if trait == "poisson":
        paths.append(result_dir / "cosmic_poisson.refined.glm.tsv")
    rows = []
    for path in sorted(set(paths)):
        if path.exists():
            rows.extend(read_table(path))
    if not rows:
        return {"refinement_methods": "NA", "refinement_rows": "NA", "refinement_converged": "NA", "refinement_errors": "NA"}
    methods = sorted(set(sval(r, "method") or "unknown" for r in rows))
    converged = sum(1 for r in rows if sval(r, "converged") in ("TRUE", "True", "true", "1"))
    errors = sum(1 for r in rows if sval(r, "error_message") not in ("", "NA"))
    return {"refinement_methods": ",".join(methods), "refinement_rows": len(rows), "refinement_converged": converged, "refinement_errors": errors}


def summarize_simulation(result_dir, trait):
    type1_path = result_dir / "psr_type1_summary.tsv"
    power_path = result_dir / "psr_power_summary.tsv"
    metrics_path = result_dir / "method_metrics.tsv"
    out = {
        "simulation_engine": "NA",
        "simulation_type1_1e_2_mean": "NA",
        "simulation_type1_1e_3_mean": "NA",
        "simulation_power_p_lt_0_05_mean": "NA",
        "simulation_runtime_mean": "NA",
        "simulation_status_classes": "NA",
    }
    if type1_path.exists():
        rows = [r for r in read_table(type1_path) if sval(r, "trait_type") == trait]
        if rows:
            out["simulation_engine"] = sval(rows[0], "engine")
            vals_1e2 = [fval(r, "type1_1e_2") for r in rows]
            vals_1e3 = [fval(r, "type1_1e_3") for r in rows]
            out["simulation_type1_1e_2_mean"] = sum(vals_1e2) / len(vals_1e2)
            out["simulation_type1_1e_3_mean"] = sum(vals_1e3) / len(vals_1e3)
    if power_path.exists():
        rows = [r for r in read_table(power_path) if sval(r, "trait_type") == trait]
        vals = [fval(r, "power_p_lt_0_05") for r in rows]
        if vals:
            out["simulation_power_p_lt_0_05_mean"] = sum(vals) / len(vals)
    if metrics_path.exists():
        rows = [r for r in read_table(metrics_path) if sval(r, "trait_type") == trait]
        runtimes = [fval(r, "runtime") for r in rows]
        runtimes = [x for x in runtimes if math.isfinite(x)]
        if runtimes:
            out["simulation_runtime_mean"] = sum(runtimes) / len(runtimes)
        classes = sorted(set(sval(r, "status_class") for r in rows if sval(r, "status_class")))
        if classes:
            out["simulation_status_classes"] = ",".join(classes)
    return out


def write_report(result_dir, rows, weighted):
    out = result_dir / "psr_validation_report.md"
    with out.open("w", encoding="utf-8") as f:
        f.write("# PSR-GWAS Validation Report\n\n")
        f.write("## Supported Trait Types\n\n")
        f.write("- Binary threshold PSR-GWAS (`--threshold-model --link probit|logit`)\n")
        f.write("- Poisson PSR-GWAS (`--count-model --family poisson`)\n")
        f.write("- Binomial event/trial PSR-GWAS (`--count-model --family binomial`)\n\n")
        f.write("## PSR Parameters\n\n")
        f.write("`--psr-maxit`, `--psr-tol`, `--psr-mu-eps`, `--psr-min-weight`, and `--psr-max-weight` control the PSR/PQL iteration and clipping behavior.\n\n")
        f.write("## Convergence Status\n\n")
        f.write("| Trait | Converged | Status | Warning | Iterations |\n")
        f.write("|---|---:|---|---|---:|\n")
        for row in rows:
            f.write(f"| {row['trait']} | {row['converged']} | {row['status_class']} | {row['warning_flag']} | {row['n_iter']} |\n")
        f.write("\n## Weighted Equivalence\n\n")
        f.write(f"- Max absolute difference: `{weighted['weighted_max_abs_diff']}`\n")
        f.write(f"- All tests passed: `{weighted['weighted_all_pass']}`\n")
        f.write(f"- Tests: `{weighted['weighted_tests']}`\n\n")
        f.write("## Simulation Smoke Overview\n\n")
        f.write("| Trait | Engine | Type I 1e-2 mean | Power mean | Runtime mean | Status classes |\n")
        f.write("|---|---|---:|---:|---:|---|\n")
        for row in rows:
            f.write(f"| {row['trait']} | {row['simulation_engine']} | {row['simulation_type1_1e_2_mean']} | {row['simulation_power_p_lt_0_05_mean']} | {row['simulation_runtime_mean']} | {row['simulation_status_classes']} |\n")
        f.write("\n## Candidate Refinement Status\n\n")
        f.write("| Trait | Methods | Rows | Converged | Errors |\n")
        f.write("|---|---|---:|---:|---:|\n")
        for row in rows:
            f.write(f"| {row['trait']} | {row['refinement_methods']} | {row['refinement_rows']} | {row['refinement_converged']} | {row['refinement_errors']} |\n")
        f.write("\n## Known Warnings And Limits\n\n")
        f.write("- PSR-GWAS is a PSR/PQL approximation, not an exact GLMM.\n")
        f.write("- The backend is pseudo response + working weights + weighted LMM scan.\n")
        f.write("- `CONVERGED_WITH_WARNING` is diagnostic; `NOT_CONVERGED` requires formal-analysis caution.\n")
        f.write("- Candidate refinement is script-based. `firth`, `glmer`, and `glmmTMB` require optional R packages.\n")
        f.write("- Negative binomial remains experimental and is not recommended for production claims.\n")
        f.write("- Ordinal, repeated-measure, and random-regression PSR-GWAS remain roadmap items.\n")
    print(f"[PSR-GWAS] Validation report written to {out}")


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: summarize_psr_gwas_benchmark.py <result_dir>")
    result_dir = Path(sys.argv[1])
    traits = ["binary", "poisson", "binomial"]
    rows = []
    conv_rows = []
    weighted = read_weighted_equivalence(result_dir)
    for trait in traits:
        cosmic_path = result_dir / f"cosmic_{trait}.assoc.txt"
        r_path = result_dir / f"r_baseline_{trait}.tsv"
        if not cosmic_path.exists() or not r_path.exists():
            continue
        cosmic = {sval(r, "SNP", "snp"): r for r in read_table(cosmic_path)}
        rbase = {sval(r, "SNP", "snp"): r for r in read_table(r_path)}
        common = sorted(set(cosmic) & set(rbase))
        cbeta, rbeta, cse, rse, cp, rp = [], [], [], [], [], []
        for snp in common:
            c = cosmic[snp]
            r = rbase[snp]
            cbeta.append(fval(c, "BETA", "beta"))
            rbeta.append(fval(r, "beta", "BETA"))
            cse.append(fval(c, "SE", "se"))
            rse.append(fval(r, "se", "SE"))
            pc = fval(c, "P", "P_VALUE", "p")
            pr = fval(r, "p", "P", "p_value", "P_VALUE")
            cp.append(-math.log10(pc) if pc > 0 and math.isfinite(pc) else math.nan)
            rp.append(-math.log10(pr) if pr > 0 and math.isfinite(pr) else math.nan)
        top_c = [sval(r, "SNP", "snp") for r in sorted(cosmic.values(), key=lambda x: fval(x, "P", "P_VALUE", "p"))[:20]]
        top_r = [sval(r, "SNP", "snp") for r in sorted(rbase.values(), key=lambda x: fval(x, "p", "P", "p_value", "P_VALUE"))[:20]]
        timing = read_timing(result_dir / f"timing.{trait}.tsv")
        cosmic_seconds = timing.get(f"cosmic_{trait}_seconds", "NA")
        conv = read_convergence(result_dir, trait) or {}
        refinement = summarize_refinement(result_dir, trait)
        simulation = summarize_simulation(result_dir, trait)
        if conv:
            conv_rows.append(conv)
        row = {
            "trait": trait,
            "overlap": len(common),
            "beta_corr": corr(cbeta, rbeta),
            "se_corr": corr(cse, rse),
            "log10p_corr": corr(cp, rp),
            "top20_overlap": len(set(top_c) & set(top_r)),
            "cosmic_seconds": cosmic_seconds,
            "n_iter": conv.get("n_iter", "NA"),
            "converged": conv.get("converged", "NA"),
            "final_sigma_g2": conv.get("final_sigma_g2", "NA"),
            "final_delta_sigma_g2": conv.get("final_delta_sigma_g2", "NA"),
            "min_mu": conv.get("min_mu", "NA"),
            "max_mu": conv.get("max_mu", "NA"),
            "min_weight": conv.get("min_weight", "NA"),
            "max_weight": conv.get("max_weight", "NA"),
            "total_clipped": conv.get("total_clipped", "NA"),
            "warning_flag": conv.get("warning_flag", "NA"),
            "status_class": conv.get("status_class", "NA"),
        }
        row.update(weighted)
        row.update(refinement)
        row.update(simulation)
        rows.append(row)

    out = result_dir / "psr_gwas_benchmark.summary.tsv"
    fields = ["trait", "overlap", "beta_corr", "se_corr", "log10p_corr", "top20_overlap", "cosmic_seconds",
              "n_iter", "converged", "final_sigma_g2", "final_delta_sigma_g2", "min_mu", "max_mu",
              "min_weight", "max_weight", "total_clipped", "warning_flag", "status_class",
              "weighted_max_abs_diff", "weighted_all_pass", "weighted_tests",
              "refinement_methods", "refinement_rows", "refinement_converged", "refinement_errors",
              "simulation_engine", "simulation_type1_1e_2_mean", "simulation_type1_1e_3_mean",
              "simulation_power_p_lt_0_05_mean", "simulation_runtime_mean", "simulation_status_classes"]
    with out.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    conv_out = result_dir / "psr_convergence_summary.tsv"
    if conv_rows:
        conv_fields = ["trait", "n_iter", "converged", "final_sigma_g2", "final_delta_sigma_g2", "final_quasiLik",
                       "min_mu", "max_mu", "min_weight", "max_weight", "total_clipped", "warning_flag", "status_class"]
        with conv_out.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, delimiter="\t", fieldnames=conv_fields)
            writer.writeheader()
            writer.writerows(conv_rows)
        print(f"[PSR-GWAS] Convergence summary written to {conv_out}")
    write_report(result_dir, rows, weighted)
    print(f"[PSR-GWAS] Benchmark summary written to {out}")


if __name__ == "__main__":
    main()
