#!/usr/bin/env python3
"""E3.2b full internal probability -> scanner -> unified-output validation."""
import argparse
import csv
import json
import math
import os
import struct
import subprocess
import tempfile
import zlib
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy import stats
from run_bgen_ad_null_validation import phenotype_matrix, write_hardcalls, probability_oracle
from run_bgen_dominance_production_gate import posterior_probabilities, read_meta, require, write_bgen, write_tsv
from run_bgen_ad_scan_validation import columns

HEADER = ["CHR", "SNP", "BP", "A1", "A2", "FREQ", "N", "BETA_A", "SE_A", "P_A", "LOG10P_A",
          "BETA_D", "SE_D", "P_D", "LOG10P_D", "P_JOINT", "LOG10P_JOINT"]


def associations(path):
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames == HEADER, "unified 17-column schema changed")
        result = list(reader)
    for row in result:
        require(all(math.isfinite(float(row[k])) for k in HEADER[5:]), "nonfinite unified output")
        require(all(0 <= float(row[k]) <= 1 for k in ("P_A", "P_D", "P_JOINT")), "invalid unified P")
    return {r["SNP"]: r for r in result}


def mutated(source, destination, kind):
    data = bytearray(source.read_bytes())
    offset = 4 + struct.unpack_from("<I", data, 0)[0]
    for _ in range(3):
        length = struct.unpack_from("<H", data, offset)[0]; offset += 2 + length
    offset += 4
    allele_offset = offset; alleles = struct.unpack_from("<H", data, offset)[0]; offset += 2
    for _ in range(alleles):
        length = struct.unpack_from("<I", data, offset)[0]; offset += 4 + length
    block_start = offset; c, d = struct.unpack_from("<II", data, offset)
    body = bytearray(zlib.decompress(data[offset+8:offset+4+c]))
    n = struct.unpack_from("<I", body, 0)[0]
    if kind == "phased": body[8+n] = 1
    elif kind == "ploidy": body[6] = body[7] = 1; body[8] = 1
    elif kind == "bits": body[9+n] = 12
    elif kind == "multiallelic":
        struct.pack_into("<H", data, allele_offset, 3)
        data[block_start:block_start] = struct.pack("<I", 1) + b"G"
        block_start += 5; struct.pack_into("<H", body, 4, 3)
    elif kind == "uncompressed":
        flags = struct.unpack_from("<I", data, 20)[0]; struct.pack_into("<I", data, 20, flags & ~3)
        destination.write_bytes(data); return
    else: raise ValueError(kind)
    compressed = zlib.compress(body)
    old_end = block_start + 4 + c
    data[block_start:old_end] = struct.pack("<II", len(compressed)+4, len(body)) + compressed
    destination.write_bytes(data)


def flip_onehot(source, destination):
    data = source.read_bytes(); offset = 4 + struct.unpack_from("<I", data, 0)[0]
    output = bytearray(data[:offset]); m = struct.unpack_from("<I", data, 8)[0]
    for _ in range(m):
        start = offset
        for _ in range(3): length = struct.unpack_from("<H", data, offset)[0]; offset += 2+length
        offset += 6
        allele_start = offset
        length1 = struct.unpack_from("<I", data, offset)[0]; allele1 = data[offset:offset+4+length1]; offset += 4+length1
        length2 = struct.unpack_from("<I", data, offset)[0]; allele2 = data[offset:offset+4+length2]; offset += 4+length2
        output.extend(data[start:allele_start]); output.extend(allele2); output.extend(allele1)
        c, d = struct.unpack_from("<II", data, offset)
        body = bytearray(zlib.decompress(data[offset+8:offset+4+c]))
        n = struct.unpack_from("<I", body, 0)[0]; bits = body[9+n]; begin = 10+n
        denominator = (1 << bits)-1
        for i in range(n):
            if bits == 8:
                q2, q1 = body[begin+2*i:begin+2*i+2]; body[begin+2*i] = denominator-q2-q1
            else:
                q2, q1 = struct.unpack_from("<HH", body, begin+4*i)
                struct.pack_into("<H", body, begin+4*i, denominator-q2-q1)
        compressed = zlib.compress(body)
        output.extend(struct.pack("<II", len(compressed)+4, len(body))); output.extend(compressed)
        offset += 4+c
    destination.write_bytes(output)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--harness", required=True); parser.add_argument("--scanner", required=True)
    parser.add_argument("--null-driver", required=True); parser.add_argument("--workdir", required=True)
    parser.add_argument("--reps", type=int, default=100); parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--plink2", default=os.environ.get("PLINK2_BIN", "plink2"))
    args = parser.parse_args(); work = Path(args.workdir); work.mkdir(parents=True, exist_ok=True)
    if args.smoke:
        # Keep the CTest smoke deterministic and small; R31 retains the full
        # production-scale fitted-null calibration evidence.
        args.reps = min(args.reps, 4)
    run = Path(tempfile.mkdtemp(prefix="e32b-", dir=work))
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    def execute(command, expect_success=True):
        p = subprocess.run([str(v) for v in command], env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if expect_success: require(p.returncode == 0, p.stdout[-6000:])
        return p
    n, m, causal = 160, 120, 15; rng = np.random.default_rng(3203)
    samples = [f"I{i:04d}" for i in range(n)]; freq = rng.uniform(.15, .45, m)
    latent = rng.binomial(2, freq[:, None], (m, n)); missing = rng.random((m, n)) < .02
    latent[0] = 0; latent[1] = 1; missing[1] = False; missing[2] = True
    missing[3, :14] = True; missing[causal] = False
    prefix = run / "hardcall"; write_hardcalls(prefix, np.where(missing, -1, latent), samples)
    execute([args.plink2, "--bfile", prefix, "--make-pgen", "--out", str(prefix)+"_pgen"])
    execute([args.plink2, "--bfile", prefix, "--export", "bgen-1.2", "--out", str(prefix)+"_bgen"])
    paths = {"one_hot": Path(str(prefix)+"_bgen.bgen")}; probabilities = {"one_hot": np.eye(3)[latent]}
    noise = rng.normal(size=(m, n))
    for quality, tau in (("high_confidence", .2), ("medium_confidence", .65), ("low_confidence", 1.2)):
        pi = posterior_probabilities(latent, freq, tau, noise, np.ones((m,n)))
        pi[0] = (1,0,0); pi[1] = (0,1,0)
        pi[7] = ((1-freq[7])**2, 2*freq[7]*(1-freq[7]), freq[7]**2)
        paths[quality] = run / (quality+".bgen"); probabilities[quality] = write_bgen(paths[quality], samples, pi, missing)
    cov = np.column_stack([np.ones(n), rng.normal(size=n)]).astype(np.float32).astype(float)
    covpath = run / "cov.txt"; phenotype_matrix(covpath, cov)
    def harness(quality, values, label, config=None, covariates=cov, expected=True, source=None):
        directory = run / label; directory.mkdir()
        pheno = run / (label+".pheno"); cp = run / (label+".cov"); cfg = run / (label+".options")
        phenotype_matrix(pheno, values); phenotype_matrix(cp, covariates)
        settings = {"maxit": 200, **(config or {})}
        cfg.write_text("".join(f"{k}\t{v}\n" for k,v in settings.items()))
        p = execute([args.harness, source or paths[quality], pheno, cp, directory, cfg], expected)
        (run/(label+".log")).write_text(p.stdout)
        if expected:
            for rep in range(values.shape[1]):
                meta = read_meta(directory/f"rep{rep}.scan.meta.txt")
                require(meta["RUN_STATUS"] == "SUCCESS" and meta["NULL_CONVERGED"] == "true", "integration status")
                require(meta["PER_SNP_NULL_REFITS"] == "0", "per-SNP refit")
                require(meta["GENOTYPE_BACKEND"] == "BGEN" and meta["PROBABILITY_MODE"] == "true", "backend metadata")
                require(meta["QC_POLICY"] == "BGEN_DOM_V1_SHARED_AD_RESEARCH", "QC policy")
                associations(directory/f"rep{rep}.assoc.ad.txt")
        else:
            require(p.returncode != 0 and not (directory/"rep0.assoc.ad.txt").exists(), "failure emitted normal GWAS")
            meta = read_meta(directory/"rep0.scan.meta.txt"); require(meta["RUN_STATUS"] == "FAILED", "failure diagnostics")
        return directory
    preliminary = rng.normal(size=(n,1))
    seed_dir = harness("one_hot", preliminary, "seed", {"loco":"false"})
    ka, kd = [np.fromfile(seed_dir/("rep0"+s),dtype="<f8").reshape(n,n) for s in (".ka.bin",".kd.bin")]
    # Interior-rich phenotype for real native-CLI one-hot covariance parity.
    healthy = np.linalg.cholesky(.8*ka+.8*kd+.2*np.eye(n)) @ rng.normal(size=(n,1)) + .4*cov[:,1,None]
    parity_results = []
    parity_dirs = {}
    for loco in (False,True):
        label = "parity_loco" if loco else "parity_global"
        bgen = harness("one_hot", healthy, label, {"loco":str(loco).lower()}); parity_dirs[loco] = bgen
        accepted = bgen/"rep0.ad.kernel.markers.txt"
        native = run/(label+"_native")
        execute([args.plink2,"--bfile",prefix,"--extract",accepted,"--make-bed","--out",native])
        execute([args.plink2,"--bfile",native,"--make-pgen","--out",str(native)+"_pgen"])
        bgen = harness("one_hot", healthy, label+"_common_null",
            {"loco":str(loco).lower(),"native_bed":native,"native_pgen":str(native)+"_pgen.pgen"})
        parity_dirs[loco] = bgen
        pheno = run/(label+"_native.pheno"); qc = run/(label+"_native.qcov")
        pheno.write_text("".join(f"{s} {s} {healthy[i,0]:.17g}\n" for i,s in enumerate(samples)))
        qc.write_text("".join(f"{s} {s} {cov[i,1]:.17g}\n" for i,s in enumerate(samples)))
        b_rows = associations(bgen/"rep0.assoc.ad.txt")
        for fmt, flag, source in (("bed","--bfile",native),("pgen","--pfile",str(native)+"_pgen")):
            out = run/(label+"_"+fmt)
            execute([args.scanner,flag,source,"--gwas-model","additive_dominance","--pheno",pheno,"--qcovar",qc,
                     "--maf",".01","--reml-maxit","200","--reml-tol","1e-5","--out",out]+([] if loco else ["--no-loco"]))
            native_rows = associations(str(out)+".assoc.ad.txt"); require(set(native_rows)==set(b_rows),"cross-format marker set")
            differences = {}
            for key in ("BETA_A","BETA_D","SE_A","SE_D","P_A","P_D","P_JOINT"):
                differences[key] = max(abs(float(b_rows[s][key])-float(native_rows[s][key])) for s in b_rows)
            controlled = associations(bgen/("rep0.reference."+fmt+".ad.txt"))
            controlled_diff = {key:max(abs(float(b_rows[s][key])-float(controlled[s][key])) for s in b_rows)
                               for key in differences}
            require(set(controlled)==set(b_rows) and max(controlled_diff.values()) < 2e-6,
                    "one-hot controlled-null beta/SE/P parity")
            require(all(b_rows[s][k]==native_rows[s][k] for s in b_rows for k in ("CHR","BP","A1","A2","N")),"signed allele metadata parity")
            native_meta = read_meta(str(out)+".scan.meta.txt"); b_meta = read_meta(bgen/"rep0.scan.meta.txt")
            vc_diff = max(abs(float(native_meta[k])-float(b_meta[k])) for k in native_meta if
                          k.startswith(("SIGMA_A2","SIGMA_D2","SIGMA_E2")) and native_meta[k]!="PER_NULL_MODEL")
            require(vc_diff<5e-4,"one-hot VCE parity")
            parity_results.append({"loco":loco,"format":fmt,"vc_max_difference":vc_diff,
                "native_cli_statistic_differences":differences,"controlled_null_differences":controlled_diff})
            # Actual native builder, for both global and outside-chromosome pairs.
            for chr_filter in ("-","1","2"):
                matrix_out = run/f"{label}_{fmt}_k{chr_filter}"
                execute([args.null_driver,fmt,str(native)+"_pgen.pgen" if fmt=="pgen" else native,
                         matrix_out,"-",chr_filter,"0","-","0"])
                for suffix in (".ka.bin",".kd.bin"):
                    actual = np.fromfile(str(matrix_out)+suffix,dtype="<f8").reshape(n,n)
                    expected_file = bgen/("rep0"+suffix) if chr_filter=="-" else bgen/("rep0.loco."+("2" if chr_filter=="1" else "1")+suffix)
                    if not loco and chr_filter!="-": continue
                    expected_matrix=np.fromfile(expected_file,dtype="<f8").reshape(n,n)
                    require(np.max(np.abs(actual-expected_matrix))<2e-6,"native global/LOCO kernel parity")
    # Direct BGEN outside rebuild must equal stored LOCO subtraction <1e-10.
    for chromosome in ("1","2"):
        ids = [line.strip() for line in (parity_dirs[True]/"rep0.ad.kernel.markers.txt").read_text().splitlines()
               if (int(line[2:])-1 < m//2)==(chromosome=="1")]
        extract = run/("chr"+chromosome+".extract"); extract.write_text("\n".join(ids)+"\n")
        direct = harness("one_hot",healthy,"direct"+chromosome,{"loco":"false","extract":extract})
        for suffix in (".ka.bin",".kd.bin"):
            a=np.fromfile(direct/("rep0"+suffix),dtype="<f8"); b=np.fromfile(parity_dirs[True]/("rep0.loco."+("2" if chromosome=="1" else "1")+suffix),dtype="<f8")
            require(np.max(np.abs(a-b))<1e-10,"BGEN LOCO subtraction vs direct rebuild")
    flipped=run/"flipped.bgen"; flip_onehot(paths["one_hot"],flipped)
    flipdir=harness("one_hot",healthy,"orientation_flip",source=flipped)
    original_rows=associations(parity_dirs[True]/"rep0.assoc.ad.txt"); flip_rows=associations(flipdir/"rep0.assoc.ad.txt")
    for snp,r in original_rows.items():
        f=flip_rows[snp]; require(f["A1"]==r["A2"] and f["A2"]==r["A1"],"allele flip labels")
        require(abs(float(f["BETA_A"])+float(r["BETA_A"]))<2e-5 and abs(float(f["BETA_D"])-float(r["BETA_D"]))<2e-5,"allele flip beta signs")
        require(abs(float(f["FREQ"])+float(r["FREQ"])-1)<1e-8,"allele flip frequency")
    # Missing phenotype/covariates and list filters, with physical-subset oracle.
    ys=healthy.copy(); cs=cov.copy(); ys[::7]=np.nan; cs[3::13,1]=np.nan
    kept=np.isfinite(ys[:,0]) & np.isfinite(cs).all(axis=1)
    extract=run/"filter.extract"; exclude=run/"filter.exclude"
    extract.write_text("\n".join(f"rs{j+1:06d}" for j in range(m) if j%3!=0)+"\n")
    exclude.write_text("rs000017\nrs000021\n")
    filtered=harness("one_hot",ys,"filtered_subset",{"loco":"false","extract":extract,"exclude":exclude},covariates=cs)
    physical=run/"physical_subset.bgen"; write_bgen(physical,[s for s,k in zip(samples,kept) if k],
        probabilities["one_hot"][:,kept],missing[:,kept])
    # Use PLINK2 physical subset for exact one-hot quantization/order equivalence.
    keepfile=run/"keep.txt"; keepfile.write_text("".join(f"{s} {s}\n" for s,k in zip(samples,kept) if k))
    execute([args.plink2,"--bfile",prefix,"--keep",keepfile,"--export","bgen-1.2","--out",run/"physical"])
    physicaldir=harness("one_hot",ys[kept],"physical_subset",{"loco":"false","extract":extract,"exclude":exclude},covariates=cs[kept],source=run/"physical.bgen")
    fr=associations(filtered/"rep0.assoc.ad.txt"); pr=associations(physicaldir/"rep0.assoc.ad.txt")
    require(fr==pr,"analysis-sample subset must equal materialized subset")
    require(all(int(r["N"])==int(kept.sum()) for r in fr.values()),"effective analysis N")
    # Unsupported/failure guards must preserve diagnostics and produce no GWAS.
    guard_results=[]
    harness("one_hot",healthy,"guard_empty_loco",{"extract":run/"chr1.extract"},expected=False)
    guard_results.append("empty_loco_outside")
    rankdir=harness("one_hot",healthy,"rank_deficient_validation",{"loco":"false","min_info_d":0})
    rank_meta=read_meta(rankdir/"rep0.scan.meta.txt")
    require(int(rank_meta["RANK_DEFICIENT_MARKERS"])>=1 and "rs000002" not in associations(rankdir/"rep0.assoc.ad.txt"),
            "rank-deficient SNP emitted a joint result")
    for key,settings in (("hwe",{"hwe":.01}), ("population",{"population":"unsupported"}),
                         ("weights",{"weights":"unsupported"}), ("unconverged",{"maxit":0})):
        harness("one_hot",healthy,"guard_"+key,settings,expected=False); guard_results.append(key)
    for kind in ("phased","ploidy","multiallelic","bits","uncompressed"):
        bad=run/(kind+".bgen"); mutated(paths["high_confidence"],bad,kind)
        harness("one_hot",healthy,"guard_"+kind,expected=False,source=bad); guard_results.append(kind)
    # Ordinary CLI must still reject BGEN AD; no internal feature flag exposed.
    cli=execute([args.scanner,"--bgen",paths["one_hot"],"--gwas-model","additive_dominance","--out",run/"public_guard"],False)
    require(cli.returncode!=0 and "not supported for BGEN" in cli.stdout,"ordinary CLI guard removed")

    # Strict nulls used by the current contract. Global null uses global KA/KD;
    # matched LOCO null uses the outside-chromosome KA/KD for target chromosome 1.
    # The full-genome polygenic case below is retained only as a stress diagnostic.
    strict_reps = args.reps
    strict_qualities = list(paths) if not args.smoke else ["one_hot"]
    strict_rows=[]; strict_summary=[]
    strict_rng=np.random.default_rng(3301)
    ka_loco1,kd_loco1=[np.fromfile(parity_dirs[True]/("rep0.loco.1"+s),dtype="<f8").reshape(n,n)
                       for s in (".ka.bin",".kd.bin")]
    truth=.25,.20,.55
    def calibration_metrics(label, directory):
        for rep in range(strict_reps):
            table=associations(directory/f"rep{rep}.assoc.ad.txt")
            require(bool(table), f"strict {label} produced no associations")
            p=np.array([float(row["P_JOINT"]) for row in table.values()])
            strict_rows.append({"quality":label.split("|")[0],"case":label.split("|")[1],"replicate":rep,
                "type1":float(np.mean(p<.05)),
                "lambda_gc":float(np.median(stats.chi2.isf(np.clip(p,1e-300,1),2))/stats.chi2.ppf(.5,2))})
        rows=[r for r in strict_rows if r["quality"]==label.split("|")[0] and r["case"]==label.split("|")[1]]
        return {"quality":label.split("|")[0],"case":label.split("|")[1],"replicates":len(rows),
                "mean_pooled_type1":float(np.mean([r["type1"] for r in rows])),
                "median_replicate_lambda_gc":float(np.median([r["lambda_gc"] for r in rows]))}
    for quality in strict_qualities:
        global_y=np.linalg.cholesky(truth[0]*ka+truth[1]*kd+truth[2]*np.eye(n))@strict_rng.standard_normal((n,strict_reps))+.4*cov[:,1,None]
        global_dir=harness(quality,global_y,"strict_global_"+quality,{"loco":"false"})
        strict_summary.append(calibration_metrics(quality+"|strict_global",global_dir))
        matched_y=np.linalg.cholesky(truth[0]*ka_loco1+truth[1]*kd_loco1+truth[2]*np.eye(n))@strict_rng.standard_normal((n,strict_reps))+.4*cov[:,1,None]
        matched_dir=harness(quality,matched_y,"strict_matched_loco1_"+quality,{"loco":"true","chromosome":"1"})
        strict_summary.append(calibration_metrics(quality+"|strict_matched_loco1",matched_dir))
    if args.smoke:
        strict_ok=all(math.isfinite(r["mean_pooled_type1"]) and math.isfinite(r["median_replicate_lambda_gc"])
                      for r in strict_summary) and len(strict_summary)>=2
    else:
        strict_ok=all(.025<=r["mean_pooled_type1"]<=.075 and .85<=r["median_replicate_lambda_gc"]<=1.15
                      for r in strict_summary) and len(strict_summary)==2*len(paths)

    # End-to-end matrix: actual LOCO/writer runs, not engine-only re-scoring.
    xa,xd=columns(probabilities["one_hot"],missing)
    residual=rng.normal(size=(n,args.reps)); stress_y=np.linalg.cholesky(.25*ka+.20*kd+.55*np.eye(n))@residual+.4*cov[:,1,None]
    tasks=[]
    for quality in paths:
        for arm,beta in (("null",(0.,0.)),("additive",(.6,0.)),("dominance",(0.,.8)),("mixed",(.6,.8))):
            y=stress_y+beta[0]*xa[causal,:,None]+beta[1]*xd[causal,:,None]
            tasks.append((quality,arm,beta,y))
    result_rows=[]; uncertainty_rows=[]
    def run_task(t):
        quality,arm,beta,y=t
        return t,harness(quality,y,"matrix_"+quality+"_"+arm)
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures=[executor.submit(run_task,t) for t in tasks]
        for future in as_completed(futures):
            (quality,arm,beta,y),directory=future.result()
            for rep in range(args.reps):
                table=associations(directory/f"rep{rep}.assoc.ad.txt"); r=table.get("rs000016")
                require(r is not None,"causal target filtered unexpectedly")
                meta=read_meta(directory/f"rep{rep}.scan.meta.txt")
                result_rows.append({"quality":quality,"scenario":arm,"replicate":rep,
                    "beta_a":float(r["BETA_A"]),"beta_d":float(r["BETA_D"]),"se_a":float(r["SE_A"]),"se_d":float(r["SE_D"]),
                    "p_a":float(r["P_A"]),"p_d":float(r["P_D"]),"p_joint":float(r["P_JOINT"]),
                    "rank":1+sum(float(row["P_JOINT"])<float(r["P_JOINT"]) for row in table.values()),
                    "weak_nulls":int(meta["WEAK_INFORMATION_NULLS"]),"converged":True})
                if arm=="null":
                    p=np.array([float(row["P_JOINT"]) for row in table.values()])
                    uncertainty_rows.append({"quality":quality,"replicate":rep,"type1":float(np.mean(p<.05)),
                        "lambda_gc":float(np.median(stats.chi2.isf(np.clip(p,1e-300,1),2))/stats.chi2.ppf(.5,2))})
    # Four working-null variance-recovery/boundary scenarios use full pipeline
    # without LOCO, where the generating covariance coefficients are comparable.
    boundary=[]
    br=max(3,min(30,args.reps))
    for label,truth in (("a_zero",(0.,.5,.5)),("d_zero",(.5,0.,.5)),("ad_zero",(0.,0.,1.)),("mixed",(.3,.3,.4))):
        y=np.linalg.cholesky(truth[0]*ka+truth[1]*kd+truth[2]*np.eye(n))@rng.normal(size=(n,br))
        directory=harness("one_hot",y,"boundary_"+label,{"loco":"false"})
        fits=[]
        for rep in range(br):
            with (directory/f"rep{rep}.ad.null.tsv").open() as f: fits+=list(csv.DictReader(f,delimiter="\t"))
        estimates=np.array([[float(r[k]) for k in ("SIGMA_A2","SIGMA_D2","SIGMA_E2")] for r in fits])
        boundary.append({"scenario":label,"truth":truth,"replicates":br,"mean":estimates.mean(axis=0).tolist(),
                         "exact_a_zero":int(np.sum(estimates[:,0]==0)),"exact_d_zero":int(np.sum(estimates[:,1]==0)),
                         "converged":sum(r["CONVERGED"]=="true" for r in fits)})
    summaries=[]
    for quality in paths:
        for arm in ("null","additive","dominance","mixed"):
            selected=[r for r in result_rows if r["quality"]==quality and r["scenario"]==arm]
            summaries.append({"quality":quality,"scenario":arm,"replicates":len(selected),
                "mean_beta_a":float(np.mean([r["beta_a"] for r in selected])),"mean_beta_d":float(np.mean([r["beta_d"] for r in selected])),
                "mean_se_a":float(np.mean([r["se_a"] for r in selected])),"mean_se_d":float(np.mean([r["se_d"] for r in selected])),
                "add_rejection":float(np.mean([r["p_a"]<.05 for r in selected])),"dom_rejection":float(np.mean([r["p_d"]<.05 for r in selected])),
                "joint_rejection":float(np.mean([r["p_joint"]<.05 for r in selected])),"median_rank":float(np.median([r["rank"] for r in selected])),
                "weak_warning_fraction":float(np.mean([r["weak_nulls"]>0 for r in selected]))})
    partial_ok=all(r["dom_rejection"]<=.12 if r["scenario"]=="additive" else r["add_rejection"]<=.12
                   for r in summaries if r["scenario"] in ("additive","dominance"))
    stress_summary=[]
    for quality in paths:
        selected=[r for r in uncertainty_rows if r["quality"]==quality]
        require(len(selected)==args.reps, f"polygenic stress replicates missing for {quality}")
        stress_summary.append({"quality":quality,"mean_pooled_type1":float(np.mean([r["type1"] for r in selected])),
                              "median_replicate_lambda_gc":float(np.median([r["lambda_gc"] for r in selected]))})
    # The full-genome polygenic case is deliberately still measured and reported.
    # It is not a strict global-null FPR gate under dual-LOCO analysis.
    stress_ok=all(math.isfinite(r["mean_pooled_type1"]) and math.isfinite(r["median_replicate_lambda_gc"])
                  for r in stress_summary) and len(stress_summary)==len(paths)
    if args.smoke:
        partial_ok=True
    # Strict calibration is a CTest contract for the full integration run. The
    # smoke run verifies that the strict code path executes and emits finite metrics;
    # production-scale type-I calibration remains R31's frozen evidence.
    strict_gate=strict_ok
    if args.smoke: partial_ok=True
    write_tsv(run/"replicate_metrics.tsv",list(result_rows[0]),result_rows)
    write_tsv(run/"qtl_summary.tsv",list(summaries[0]),summaries)
    write_tsv(run/"null_calibration.tsv",list(uncertainty_rows[0]),uncertainty_rows)
    write_tsv(run/"strict_calibration.tsv",list(strict_summary[0]),strict_summary)
    summary={"verdict":("INTEGRATION_SMOKE_PASS" if args.smoke else
                        "READY FOR BGEN AD INTEGRATION") if partial_ok and stress_ok and strict_gate else "BLOCKED",
        "smoke":args.smoke,
        "n":n,"m":m,"replicates":args.reps,"source_workdir":str(run),"one_hot_parity":parity_results,
        "cli_guard":True,"guard_cases":guard_results,"subset_n":int(kept.sum()),"boundary":boundary,
        "qtl":summaries,"partial_null_no_clear_inflation":partial_ok,
        "strict_calibration":strict_summary,"strict_calibration_replicates":strict_reps,
        "strict_calibration_gate":strict_gate,
        "polygenic_stress":stress_summary,"polygenic_stress_not_strict_null":True,
        "null_calibration_generation_scope":"FULL_GENOME_POLYGENIC_COVARIANCE;ANALYSIS_DUAL_LOCO",
        "strict_null_generation_scope":{"strict_global":"GLOBAL_KA_KD","strict_matched_loco1":"LOCO1_KA_KD"}}
    (run/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    (work/"latest_run.txt").write_text(str(run)+"\n")
    print(json.dumps(summary,indent=2)); require(partial_ok and stress_ok and strict_gate,"BGEN AD integration validation blocked")


if __name__=="__main__": main()
