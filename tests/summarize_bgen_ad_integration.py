#!/usr/bin/env python3
"""Postprocess completed E3.2b outputs; never repeat null fits or scans."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from run_bgen_dominance_production_gate import require


def read(path):
    with path.open() as handle: return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--run",required=True); args=parser.parse_args()
    run=Path(args.run); summary=json.loads((run/"summary.json").read_text())
    certificates=[]; log_error=0.0; count_rows=0
    for path in run.glob("*/rep*.ad.null.tsv"):
        certificates.extend(read(path))
    certified=[r for r in certificates if r["CONVERGED"]=="true"]
    # Failed maxit=0 intentionally appears in its own guard directory.
    for r in certified:
        require(float(r["GRADIENT_NORM"])<=1e-7 and float(r["KKT_RESIDUAL"])<=1e-7,"certificate gradient/KKT")
        require(float(r["LIKELIHOOD_DELTA"])<=1e-9*max(1,abs(float(r["LOGL"]))),"certificate likelihood")
    for path in run.glob("matrix_*/rep*.assoc.ad.txt"):
        for r in read(path):
            count_rows+=1
            for p,l in (("P_A","LOG10P_A"),("P_D","LOG10P_D"),("P_JOINT","LOG10P_JOINT")):
                if float(r[p])>0: log_error=max(log_error,abs(float(r[l])+np.log10(float(r[p]))))
    require(log_error<3e-6,"stable writer LOG10P consistency")
    nulls=read(run/"null_calibration.tsv"); stress=[]
    for quality in sorted({r["quality"] for r in nulls}):
        selected=[r for r in nulls if r["quality"]==quality]
        stress.append({"quality":quality,"replicates":len(selected),
            "mean_pooled_type1":float(np.mean([float(r["type1"]) for r in selected])),
            "median_replicate_lambda_gc":float(np.median([float(r["lambda_gc"]) for r in selected]))})
    strict=read(run/"strict_calibration.tsv"); strict_pooled=[]
    for r in strict:
        strict_pooled.append({"quality":r["quality"],"case":r["case"],"replicates":int(r["replicates"]),
            "mean_pooled_type1":float(r["mean_pooled_type1"]),
            "median_replicate_lambda_gc":float(r["median_replicate_lambda_gc"])})
    native_error=0.0; direct_error=0.0
    for label in ("parity_global","parity_loco"):
        for fmt in ("bed","pgen"):
            for c in ("-","1","2"):
                if label=="parity_global" and c!="-": continue
                for suffix in (".ka.bin",".kd.bin"):
                    a=np.fromfile(run/(f"{label}_{fmt}_k{c}"+suffix),dtype="<f8")
                    name="rep0"+suffix if c=="-" else "rep0.loco."+("2" if c=="1" else "1")+suffix
                    b=np.fromfile(run/(label+"_common_null")/name,dtype="<f8")
                    native_error=max(native_error,float(np.max(np.abs(a-b))))
    for c in ("1","2"):
        for suffix in (".ka.bin",".kd.bin"):
            a=np.fromfile(run/("direct"+c)/("rep0"+suffix),dtype="<f8")
            b=np.fromfile(run/"parity_loco_common_null"/("rep0.loco."+("2" if c=="1" else "1")+suffix),dtype="<f8")
            direct_error=max(direct_error,float(np.max(np.abs(a-b))))
    diagnostics={"successful_null_certificates":len(certified),"intentional_failed_certificates":len(certificates)-len(certified),
        "max_gradient_norm":max(float(r["GRADIENT_NORM"]) for r in certified),
        "max_kkt_residual":max(float(r["KKT_RESIDUAL"]) for r in certified),
        "native_global_loco_kernel_max_error":native_error,"loco_direct_rebuild_max_error":direct_error,
        "matrix_unified_rows":count_rows,"max_log10p_consistency_error":float(log_error),
        "strict_calibration":strict_pooled,"polygenic_stress":stress,
        "polygenic_stress_not_strict_null":True}
    strict_ok=all(.025<=r["mean_pooled_type1"]<=.075 and .85<=r["median_replicate_lambda_gc"]<=1.15
                  for r in strict_pooled) if summary.get("smoke") is not True else True
    summary.setdefault("technical_harness_verdict_before_global_null_audit",summary["verdict"])
    summary["strict_calibration"]=strict_pooled
    summary["strict_calibration_gate"]=strict_ok
    summary["polygenic_stress"]=stress
    summary["polygenic_stress_not_strict_null"]=True
    summary["null_generation_scope"]=("strict_global=GLOBAL_KA_KD;strict_matched_loco=LOCO1_KA_KD;"
                                      "polygenic_stress=FULL_GENOME_POLYGENIC_COVARIANCE;ANALYSIS_DUAL_LOCO")
    summary["verdict"]="INTEGRATION_SMOKE_PASS" if summary.get("smoke") is True else (
        "READY FOR BGEN AD INTEGRATION" if strict_ok and summary["partial_null_no_clear_inflation"] else "BLOCKED")
    (run/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    diagnostics["verdict"]=summary["verdict"]
    (run/"diagnostic_summary.json").write_text(json.dumps(diagnostics,indent=2)+"\n")
    print(json.dumps(diagnostics,indent=2))
    require(strict_ok,"strict BGEN AD null/LOCO calibration gate blocked")


if __name__=="__main__": main()
