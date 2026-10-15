#!/usr/bin/env python3
import csv
import math
import os
import subprocess
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get("COSMICSCANNER_BIN", ROOT / "build" / "cosmicscanner"))
OUT_DIR = ROOT / "tests" / "results" / "psr_gwas_smoke"
WORK_DIR = OUT_DIR / "simulation_work"
ENGINE = "cosmicscanner_psr"
BASE_SEED = 20260623


def add_ext(prefix, ext):
    return Path(str(prefix) + ext)


def logistic(x):
    return 1.0 / (1.0 + np.exp(-x))


def plink_code(dosage):
    return {2: 0b00, 1: 0b10, 0: 0b11}.get(int(dosage), 0b01)


def write_plink(prefix, G):
    n, m = G.shape
    with add_ext(prefix, ".fam").open("w", encoding="utf-8") as f:
        for i in range(n):
            f.write(f"F{i+1:04d}\tI{i+1:04d}\t0\t0\t1\t-9\n")
    with add_ext(prefix, ".bim").open("w", encoding="utf-8") as f:
        for j in range(m):
            f.write(f"{(j % 4) + 1}\tsim{j+1:06d}\t0\t{1000 + j * 10}\tA\tC\n")
    bytes_per_snp = (n + 3) // 4
    bed = bytearray(b"\x6c\x1b\x01")
    for j in range(m):
        block = bytearray(bytes_per_snp)
        for i, g in enumerate(G[:, j].astype(int)):
            block[i // 4] |= plink_code(g) << ((i % 4) * 2)
        bed.extend(block)
    add_ext(prefix, ".bed").write_bytes(bed)


def write_one_col_pheno(path, y):
    with path.open("w", encoding="utf-8") as f:
        for i, val in enumerate(y):
            f.write(f"F{i+1:04d}\tI{i+1:04d}\t{val}\n")


def write_binomial_pheno(path, events, trials):
    with path.open("w", encoding="utf-8") as f:
        for i, (event, trial) in enumerate(zip(events, trials)):
            f.write(f"F{i+1:04d}\tI{i+1:04d}\t{int(event)}\t{int(trial)}\n")


def read_table(path):
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_pvals(path):
    pvals = []
    with path.open("r", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            try:
                p = float(row.get("P_A", "nan"))
            except ValueError:
                p = math.nan
            pvals.append(p if math.isfinite(p) else 1.0)
    return np.array(pvals)


def convergence_status(out_prefix):
    path = add_ext(out_prefix, ".psr.convergence.tsv")
    if not path.exists():
        return "NA", "NA", "NA"
    rows = read_table(path)
    if not rows:
        return "NA", "NA", "NA"
    converged = int(float(rows[-1].get("converged", 0) or 0))
    warnings = []
    prev = None
    total_clipped = 0
    for row in rows:
        clipped = int(float(row.get("n_clipped", 0) or 0))
        total_clipped += clipped
        if clipped > 0:
            warnings.append("WARNING_PSR_CLIPPED")
        try:
            q = float(row.get("quasiLik", "nan"))
        except ValueError:
            q = math.nan
        if prev is not None and math.isfinite(q) and math.isfinite(prev) and q < prev - 1e-8:
            warnings.append("WARNING_PSR_QUASILIK_DECREASED")
        prev = q
    try:
        sigma_g2 = float(rows[-1].get("sigma_g2", "nan"))
    except ValueError:
        sigma_g2 = math.nan
    if math.isfinite(sigma_g2) and sigma_g2 <= 1e-5:
        warnings.append("WARNING_PSR_SIGMA_BOUNDARY")
    if not converged:
        warnings.append("WARNING_PSR_MAXIT_REACHED")
    warning_flag = ",".join(sorted(set(warnings))) if warnings else "OK"
    status_class = "NOT_CONVERGED" if not converged else ("CONVERGED_WITH_WARNING" if warning_flag != "OK" else "OK")
    return converged, warning_flag, status_class


def simulate_data(rng, n, m, trait_type, effect):
    freqs = rng.uniform(0.08, 0.45, size=m)
    G = rng.binomial(2, freqs, size=(n, m)).astype(int)
    causal = m // 3
    pop = np.r_[np.zeros(n // 2), np.ones(n - n // 2)]
    eta = 0.35 * pop + 0.25 * rng.normal(size=n) + effect * (G[:, causal] - G[:, causal].mean())
    if trait_type == "binary":
        y = rng.binomial(1, logistic(eta - 0.1))
        if y.min() == y.max():
            y[0] = 0
            y[-1] = 1
        return G, causal, {"pheno": y}
    if trait_type == "poisson":
        return G, causal, {"pheno": rng.poisson(np.exp(0.25 + eta))}
    trials = rng.integers(8, 16, size=n)
    return G, causal, {"events": rng.binomial(trials, logistic(eta - 0.2)), "trials": trials}


def run_cosmic(prefix, trait_type, payload, out_prefix, link="probit"):
    cmd = [str(BIN), "--bfile", str(prefix), "--gwas", "--no-loco", "--maf", "0.0", "--psr-maxit", "25", "--out", str(out_prefix)]
    if trait_type == "binary":
        pheno = add_ext(prefix, ".pheno")
        write_one_col_pheno(pheno, payload["pheno"])
        cmd += ["--pheno", str(pheno), "--mpheno", "1", "--threshold-model", "--link", link]
    elif trait_type == "poisson":
        pheno = add_ext(prefix, ".pheno")
        write_one_col_pheno(pheno, payload["pheno"])
        cmd += ["--pheno", str(pheno), "--mpheno", "1", "--count-model", "--family", "poisson"]
    else:
        pheno = add_ext(prefix, ".binomial.pheno")
        write_binomial_pheno(pheno, payload["events"], payload["trials"])
        cmd += ["--count-model", "--family", "binomial", "--binomial-pheno", str(pheno)]
    env = os.environ.copy()
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("OMP_NUM_THREADS", "1")
    start = time.time()
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
    runtime = time.time() - start
    return read_pvals(add_ext(out_prefix, ".assoc.txt")), runtime, convergence_status(out_prefix)


def base_row(seed, trait_type, scenario, n, m, rep, runtime, conv):
    converged, warning_flag, status_class = conv
    return {
        "engine": ENGINE,
        "seed": seed,
        "trait_type": trait_type,
        "scenario": scenario,
        "n": n,
        "m": m,
        "rep": rep,
        "runtime": runtime,
        "converged": converged,
        "warning_flag": warning_flag,
        "status_class": status_class,
    }


def main():
    if not BIN.exists():
        raise SystemExit(f"Missing CosmicScanner executable: {BIN}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    traits = ["binary", "poisson", "binomial"]
    n, m, reps = 72, 80, 3
    type1_rows, power_rows, qq_rows, metrics_rows = [], [], [], []
    for trait_type in traits:
        for scenario, effect in [("null", 0.0), ("power_low", 0.35), ("power_high", 0.7)]:
            for rep in range(reps):
                seed = BASE_SEED + rep + 100 * traits.index(trait_type) + 1000 * (0 if scenario == "null" else 1 if scenario == "power_low" else 2)
                rng = np.random.default_rng(seed)
                prefix = WORK_DIR / f"{trait_type}.{scenario}.rep{rep}"
                out = WORK_DIR / f"{trait_type}.{scenario}.rep{rep}.cosmic"
                G, causal, payload = simulate_data(rng, n, m, trait_type, effect)
                write_plink(prefix, G)
                pvals, runtime, conv = run_cosmic(prefix, trait_type, payload, out)
                row = base_row(seed, trait_type, scenario, n, m, rep, runtime, conv)
                row.update({"min_p": float(np.min(pvals)), "median_neglog10p": float(np.median(-np.log10(np.clip(pvals, 1e-300, 1.0))))})
                metrics_rows.append(dict(row, n_p_lt_0_05=int(np.sum(pvals < 0.05))))
                if scenario == "null":
                    type1_rows.append(dict(row, type1_1e_2=float(np.mean(pvals < 1e-2)), type1_1e_3=float(np.mean(pvals < 1e-3))))
                    obs = np.sort(-np.log10(np.clip(pvals, 1e-300, 1.0)))
                    exp = -np.log10((np.arange(1, len(obs) + 1) - 0.5) / len(obs))[::-1]
                    for k in np.linspace(0, len(obs) - 1, 25, dtype=int):
                        qq_rows.append(dict(row, qq_index=int(k), expected_neglog10p=float(exp[k]), observed_neglog10p=float(obs[k])))
                else:
                    top = np.argsort(pvals)[:20]
                    power_rows.append(dict(row, effect=effect, causal_snp=f"sim{causal+1:06d}", causal_p=float(pvals[causal]), power_p_lt_0_05=int(pvals[causal] < 0.05), causal_top20=int(causal in top)))

    def write(path, rows):
        with path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, delimiter="\t", fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    write(OUT_DIR / "psr_type1_summary.tsv", type1_rows)
    write(OUT_DIR / "psr_power_summary.tsv", power_rows)
    write(OUT_DIR / "psr_qq_points.tsv", qq_rows)
    write(OUT_DIR / "method_metrics.tsv", metrics_rows)
    print(f"[PSR-GWAS] CosmicScanner simulation summaries written to {OUT_DIR}")


if __name__ == "__main__":
    main()
