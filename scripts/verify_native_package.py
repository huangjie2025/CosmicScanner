#!/usr/bin/env python3
"""Run the same four packaged workflows on Windows, macOS, and Linux."""
from pathlib import Path
import csv,json,math,os,subprocess,tempfile
def main():
    pkg=Path(__file__).resolve().parents[1]
    exe=pkg/'bin'/('cosmicscanner.exe' if os.name=='nt' else 'cosmicscanner')
    data=pkg/'examples/dataset'
    out=Path(tempfile.mkdtemp(prefix='cosmicscanner-native-demo-',dir=os.environ.get('COSMICSCANNER_TEST_TMPDIR')))
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    common=['--bfile',str(data/'synth'),'--threads','1']
    commands={
      'lmm':['gwas','--model','lmm',*common,'--pheno',str(data/'quantitative.pheno'),'--mpheno','1'],
      'mv':['gwas','--model','multivariate',*common,'--pheno',str(data/'multivariate.pheno'),'--mpheno-list','1,2,3','--scan-tier','approx','--reml-maxit','1200','--reml-tol','1e-10','--no-loco'],
      'rep':['gwas','--model','repeatability',*common,'--pheno',str(data/'long.pheno'),'--tpheno','1','--mpheno','2','--repeatability-mode','exact','--scan-tier','approx','--no-loco','--reml-vg','.5','--reml-vpe','.2','--reml-ve','.3'],
      'grm':['grm',*common]}
    for name,flags in commands.items():
        with (out/(name+'.log')).open('w') as log:
            result=subprocess.run([str(exe),*flags,'--out',str(out/name)],env=env,stdout=log,stderr=subprocess.STDOUT,timeout=180)
        if result.returncode:raise RuntimeError((out/(name+'.log')).read_text(encoding="utf-8",errors="replace")[-6000:])
    counts={}
    for name,suffix in [('lmm','.assoc.txt'),('mv','.assoc.mv.txt'),('rep','.assoc.txt')]:
        with (out/(name+suffix)).open(encoding="utf-8") as f:rows=list(csv.DictReader(f,delimiter='\t'))
        assert len(rows)==180,(name,len(rows))
        pfields=[k for k in rows[0] if k in ('P','P_A','P_JOINT','P_MV','P_JOINT_GLS','P_JOINT_WALD') or k.lower() in ('p_wald','p_joint','pvalue')]
        assert pfields,(name,'No recognized P field',list(rows[0]))
        for row in rows:
            numeric=[float(row[k]) for k in pfields if row[k] not in ('NA','.','nan','')]
            assert numeric and all(math.isfinite(p) and 0<=p<=1 for p in numeric),(name,row)
        counts[name]=len(rows)
    assert len((out/'grm.grm.id').read_text().splitlines())==60
    lines=[s for s in (out/'lmm.cosmic.sumstats').read_text().splitlines() if s and not s.startswith('#')]
    assert len(lines)==181
    manifests=list(out.glob('*.run_manifest.json'))
    assert manifests,'Missing LMM run manifest'
    for manifest in manifests:json.loads(manifest.read_text(encoding='utf-8'))
    result={'status':'PASS','variant_counts':counts,'grm_samples':60,'results':str(out)}
    (out/'verification.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result),flush=True)
if __name__=='__main__':main()
