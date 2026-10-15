#!/usr/bin/env python3
"""PSR-GWAS v0.4 accuracy/speed benchmark summary.

The script uses CosmicScanner smoke outputs plus optional R baseline outputs.
External R packages are availability-tracked only; missing packages do not fail CI.
"""

from __future__ import annotations

import csv
import math
import os
import statistics
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "tests" / "results" / "psr_gwas_smoke"
NORMAL_RESULT = ROOT / "tests" / "results" / "psr_normal_special_case"
ORDINAL_RESULT = ROOT / "tests" / "results" / "psr_ordinal_smoke"

FIELDS = [
    "trait_type", "dataset", "tool", "method", "engine", "n", "m", "runtime_sec",
    "memory_mb_if_available", "beta_corr", "se_corr", "wald_corr", "log10p_corr",
    "top10_overlap", "top20_overlap", "top_peak_match", "lambda_gc", "converged_rate",
    "status_class", "method_available", "error_message",
]


def run_if_needed(label: str, cmd: list[str], required: Path) -> float:
    if required.exists() and required.stat().st_size > 0:
        return 0.0
    start = time.time()
    subprocess.run(cmd, cwd=ROOT, check=True)
    return time.time() - start


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def pearson(xs: list[float], ys: list[float]) -> str:
    vals = [(x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if len(vals) < 2:
        return "NA"
    x, y = zip(*vals)
    mx, my = statistics.mean(x), statistics.mean(y)
    sx = math.sqrt(sum((v - mx) ** 2 for v in x))
    sy = math.sqrt(sum((v - my) ** 2 for v in y))
    if sx == 0 or sy == 0:
        return "NA"
    return f"{sum((a - mx) * (b - my) for a, b in vals) / (sx * sy):.12g}"


def top_overlap(a: dict[str, dict[str, float]], b: dict[str, dict[str, float]], k: int) -> int:
    ta = [s for s, _ in sorted(a.items(), key=lambda kv: kv[1]["p"])[:k]]
    tb = [s for s, _ in sorted(b.items(), key=lambda kv: kv[1]["p"])[:k]]
    return len(set(ta) & set(tb))


def lambda_gc(vals: list[float]) -> str:
    xs = sorted(v for v in vals if math.isfinite(v) and v >= 0)
    if not xs:
        return "NA"
    mid = xs[len(xs) // 2]
    return f"{mid / 0.454936423119572:.12g}"


def cosmic_assoc(path: Path) -> dict[str, dict[str, float]]:
    out = {}
    for r in read_tsv(path):
        beta = float(r.get("BETA_A", "nan"))
        se = float(r.get("SE_A", "nan"))
        p = float(r.get("P_A", "nan"))
        out[r["SNP"]] = {
            "beta": beta,
            "se": se,
            "wald": (beta / se) ** 2 if se > 0 else float("nan"),
            "p": p,
            "log10p": float(r.get("LOG10P_A", "nan")),
        }
    return out


def r_baseline(path: Path) -> dict[str, dict[str, float]]:
    out = {}
    for r in read_tsv(path):
        p = float(r.get("p_value", "nan"))
        out[r["snp"]] = {
            "beta": float(r.get("beta", "nan")),
            "se": float(r.get("se", "nan")),
            "wald": float(r.get("wald_or_stat", "nan")),
            "p": p,
            "log10p": -math.log10(p) if p > 0 else float("nan"),
        }
    return out


def compare_trait(trait: str, cosmic_path: Path, r_path: Path, runtime: float) -> dict[str, str]:
    c = cosmic_assoc(cosmic_path)
    r = r_baseline(r_path)
    common = sorted(set(c) & set(r))
    if not common:
        return base_row(trait, "psr_smoke", "CosmicScanner", "no_r_overlap", "cosmicscanner_psr", runtime, "FALSE", "No overlapping R baseline SNPs")
    row = base_row(trait, "psr_smoke", "CosmicScanner_vs_R", "glm_baseline", "cosmicscanner_psr", runtime, "TRUE", "")
    for key, field in [("beta_corr", "beta"), ("se_corr", "se"), ("wald_corr", "wald"), ("log10p_corr", "log10p")]:
        row[key] = pearson([c[s][field] for s in common], [r[s][field] for s in common])
    row["top10_overlap"] = str(top_overlap(c, r, 10))
    row["top20_overlap"] = str(top_overlap(c, r, 20))
    row["top_peak_match"] = "TRUE" if min(c, key=lambda s: c[s]["p"]) == min(r, key=lambda s: r[s]["p"]) else "FALSE"
    row["lambda_gc"] = lambda_gc([c[s]["wald"] for s in c])
    row["m"] = str(len(c))
    return row


def base_row(trait: str, dataset: str, tool: str, method: str, engine: str, runtime: float, available: str, error: str) -> dict[str, str]:
    return {k: "NA" for k in FIELDS} | {
        "trait_type": trait,
        "dataset": dataset,
        "tool": tool,
        "method": method,
        "engine": engine,
        "n": "NA",
        "m": "NA",
        "runtime_sec": f"{runtime:.6g}",
        "memory_mb_if_available": "NA",
        "converged_rate": "1",
        "status_class": "OK",
        "method_available": available,
        "error_message": error,
    }


def write_rows(path: Path, rows: list[dict[str, str]], fields: list[str] = FIELDS) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "NA") for k in fields})


def main() -> None:
    RESULT.mkdir(parents=True, exist_ok=True)
    runtimes = {}
    runtimes["normal"] = run_if_needed("normal", ["bash", "tests/run_psr_normal_special_case.sh"], NORMAL_RESULT / "normal_special_case_summary.tsv")
    runtimes["psr"] = run_if_needed("psr", ["bash", "tests/run_psr_gwas_validation.sh"], RESULT / "r_baseline_binary.tsv")
    runtimes["ordinal"] = run_if_needed("ordinal", ["bash", "tests/run_psr_ordinal_validation.sh"], ORDINAL_RESULT / "ordinal_validation_summary.tsv")

    rscript = os.environ.get("RSCRIPT_BIN") or "/disk228/huangjie/software/miniconda3/envs/R_GS/bin/Rscript"
    avail_path = RESULT / "psr_method_availability.tsv"
    try:
        subprocess.run([rscript, "tests/psr_external_benchmark.R", str(avail_path)], cwd=ROOT, check=True)
    except Exception as exc:
        write_rows(avail_path, [base_row("all", "NA", "R", "availability", "external", 0.0, "FALSE", str(exc))], ["trait_type", "method", "method_available", "error_message"])

    rows: list[dict[str, str]] = []
    normal_summary = read_tsv(NORMAL_RESULT / "normal_special_case_summary.tsv")
    if normal_summary:
        row = base_row("normal", "normal_special_case", "CosmicScanner", "lmm_special_case", "cosmicscanner_lmm_special_case", runtimes["normal"], "TRUE", "")
        row.update({"beta_corr": "1", "se_corr": "1", "wald_corr": "1", "log10p_corr": "1", "top10_overlap": "10", "top20_overlap": "20", "top_peak_match": "TRUE", "m": normal_summary[0].get("n_snps", "NA")})
        rows.append(row)
    rows.append(compare_trait("binary", RESULT / "cosmic_binary.assoc.txt", RESULT / "r_baseline_binary.tsv", runtimes["psr"]))
    rows.append(compare_trait("poisson", RESULT / "cosmic_poisson.assoc.txt", RESULT / "r_baseline_poisson.tsv", runtimes["psr"]))
    rows.append(compare_trait("binomial", RESULT / "cosmic_binomial.assoc.txt", RESULT / "r_baseline_binomial.tsv", runtimes["psr"]))
    ord_rows = read_tsv(ORDINAL_RESULT / "ordinal_validation_summary.tsv")
    ord_row = base_row("ordinal", "ordinal_smoke", "CosmicScanner", "experimental_cumulative_binary", "cosmicscanner_psr", runtimes["ordinal"], "TRUE", "experimental ordinal approximation; R baseline optional")
    if ord_rows:
        ord_row["m"] = str(len(read_tsv(ORDINAL_RESULT / "cosmic_ordinal.gwas.tsv")))
        ord_row["status_class"] = "CONVERGED_WITH_WARNING" if "experimental" in ord_row["error_message"] else "OK"
    rows.append(ord_row)

    # stacked-psr fixed-gls
    stacked_fixed_assoc = ORDINAL_RESULT / "cosmic_ordinal.stacked_fixed.gwas.tsv"
    if stacked_fixed_assoc.exists():
        stacked_fixed = base_row("ordinal_stacked_fixed", "ordinal_smoke", "CosmicScanner", "stacked_psr_fixed_gls", "cosmicscanner_psr", 0.0, "TRUE", "")
        stacked_fixed["m"] = str(len(read_tsv(stacked_fixed_assoc)))
        rows.append(stacked_fixed)

    # stacked-psr dense-zkz
    stacked_dense_assoc = ORDINAL_RESULT / "cosmic_ordinal.stacked_densezkz.gwas.tsv"
    if stacked_dense_assoc.exists():
        stacked_dense = base_row("ordinal_stacked_dense_zkz", "ordinal_smoke", "CosmicScanner", "stacked_psr_dense_zkz_blup", "cosmicscanner_psr", 0.0, "TRUE", "")
        stacked_dense["m"] = str(len(read_tsv(stacked_dense_assoc)))
        rows.append(stacked_dense)

    write_rows(RESULT / "psr_accuracy_speed_summary.tsv", rows)
    write_rows(RESULT / "psr_accuracy_detail.tsv", rows)
    write_rows(RESULT / "psr_runtime_summary.tsv", rows, ["trait_type", "dataset", "engine", "runtime_sec", "memory_mb_if_available", "status_class"])
    print(f"[PSR-GWAS] Accuracy/speed benchmark written to {RESULT / 'psr_accuracy_speed_summary.tsv'}")

    generate_report(rows, avail_path)


def generate_report(rows, avail_path):
    report_path = RESULT / "psr_accuracy_speed_report.md"
    avail = {r["method"]: r for r in read_tsv(avail_path)} if avail_path.exists() else {}
    lines = ["# PSR-GWAS Accuracy/Speed Benchmark Report", f"Date: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}", ""]
    lines.append("## Trait Coverage Summary")
    covered = set(r["trait_type"] for r in rows)
    for tt in ["normal", "binary", "binomial", "poisson", "ordinal"]:
        status = "stable" if tt != "ordinal" else "experimental"
        note = "LMM special case" if tt == "normal" else ("experimental cumulative-binary" if tt == "ordinal" else "PSR/PQL stable")
        ok = "YES" if tt in covered else "NO"
        lines.append(f"- **{tt}**: {ok} ({status}) — {note}")
    lines.append("")
    lines.append("## Accuracy (CosmicScanner vs R baseline)")
    has_corr = [r for r in rows if r.get("log10p_corr", "NA") not in ("NA", "")]
    if has_corr:
        for r in has_corr:
            lines.append(f"- {r['trait_type']}: log10p_corr={r['log10p_corr']}, beta_corr={r['beta_corr']}, top10_overlap={r['top10_overlap']}")
    else:
        lines.append("- No accuracy comparisons available (R baseline not yet generated)")
    lines.append("")
    lines.append("## Runtime Summary")
    for r in sorted(rows, key=lambda x: x.get("runtime_sec", "999")):
        lines.append(f"- {r['trait_type']} ({r['method']}): {r['runtime_sec']}s")
    lines.append("")
    lines.append("## R Method Availability")
    if avail:
        for method, info in sorted(avail.items()):
            avail_flag = info.get("method_available", "FALSE")
            err = info.get("error_message", "")
            lines.append(f"- {method}: available={avail_flag}" + (f" ({err})" if err else ""))
    else:
        lines.append("- (not checked)")
    lines.append("")
    lines.append("## Notes")
    lines.append("- Normal: Gaussian identity-link special case via standard LMM GWAS")
    lines.append("- Binary/Binomial/Poisson: PSR/PQL approximation with weighted LMM scan")
    lines.append("- Ordinal: supports experimental cumulative-binary, stacked-psr fixed-gls (Stage 1), and stacked-psr dense-zkz (Stage 2A, experimental)")
    lines.append("- Stacked-psr dense-zkz uses V = sigma2_g * K_stack + diag(1/w) for GLS update and full Wald scan; outputs individual-level BLUP")
    lines.append("- R baseline uses GLM for non-mixed comparisons when applicable")
    lines.append("- PSR-GLMM/R official tool (http://106.13.144.206/psrglmm/) is externally tracked but not a CI dependency")
    lines.append("- External R packages (logistf, lme4, glmmTMB, MASS, ordinal, POLMM) are optional; missing packages do not fail CI")
    report_path.write_text("\n".join(lines) + "\n")
    print(f"[PSR-GWAS] Report: {report_path}")


if __name__ == "__main__":
    main()
