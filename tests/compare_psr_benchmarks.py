#!/usr/bin/env python3
from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[1]
summary = ROOT / "tests" / "results" / "psr_gwas_smoke" / "psr_accuracy_speed_summary.tsv"
if not summary.exists():
    raise SystemExit(f"Missing benchmark summary: {summary}")
with summary.open(newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
print("trait_type\tengine\tmethod\tlog10p_corr\truntime_sec\tstatus_class")
for r in rows:
    print("\t".join([r.get("trait_type", "NA"), r.get("engine", "NA"), r.get("method", "NA"), r.get("log10p_corr", "NA"), r.get("runtime_sec", "NA"), r.get("status_class", "NA")]))
