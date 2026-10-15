#!/usr/bin/env python3
"""R19-E: frozen binary, sequential real-data timings and output comparisons."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

import numpy as np


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def load(path):
    data = {}
    with open(path) as f:
        rows = csv.DictReader(f, delimiter='\t')
        suffix = '_A' if 'BETA_A' in rows.fieldnames else ''
        for row in rows:
            key = (row['CHR'], row['SNP'], row['BP'], row['A1'], row['A2'])
            if key in data:
                raise ValueError(f'duplicate association key: {key}')
            beta, se, p = (float(row[k+suffix]) for k in ('BETA','SE','P'))
            if not all(np.isfinite([beta,se,p])) or not 0 <= p <= 1 or se <= 0:
                raise ValueError(f'invalid association: {key}')
            # Use printed LOG10P for underflowed p, otherwise recompute -log10(P).
            logp = -np.log10(p) if p > 0 else float(row['LOG10P'+suffix])
            data[key] = (beta,se,logp)
    return data


def compare(reference, candidate):
    a, b = load(reference), load(candidate)
    keys = sorted(a.keys() & b.keys())
    if len(keys) < 2:
        raise ValueError('insufficient matched associations')
    av, bv = np.array([a[k] for k in keys]), np.array([b[k] for k in keys])
    result = {'reference_rows':len(a), 'candidate_rows':len(b), 'matched':len(keys)}
    for j, name in enumerate(('beta','se','logp')):
        result[name+'_corr'] = float(np.corrcoef(av[:,j],bv[:,j])[0,1])
        result[name+'_max_abs_difference'] = float(np.max(np.abs(av[:,j]-bv[:,j])))
    # Top sets come from the full outputs, not just the overlap.
    topa = set(sorted(a, key=lambda k:(-a[k][2],k))[:10])
    topb = set(sorted(b, key=lambda k:(-b[k][2],k))[:10])
    result['top10_overlap'] = len(topa & topb)
    result['accuracy_pass'] = (len(keys)==len(a)==len(b) and
                               result['beta_corr']>0.999 and result['logp_corr']>0.999)
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--binary', default='build/cosmicscanner')
    parser.add_argument('--out', required=True)
    parser.add_argument('--baseline', default='/tmp/r19d')
    parser.add_argument('--threads', type=int, default=32)
    args=parser.parse_args()
    out=Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    binary=out/'cosmicscanner.frozen'
    if binary.exists():
        raise RuntimeError('use a new evidence directory; frozen binary already exists')
    shutil.copy2(args.binary, binary)
    sources={str(p):sha256(p) for p in (Path('src/kinship_operator.cpp'),
             Path('src/kinship_operator.h'),Path('src/lmm_engine.cpp'),
             Path('src/gwas_scanner.cpp'),Path('src/matrix_free_vce.cpp'))}
    env=os.environ.copy()
    env.update({k:str(args.threads) for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS')})
    env.update({'MKL_DYNAMIC':'FALSE','OMP_DYNAMIC':'FALSE'})
    summary={'binary_sha256':sha256(binary),'source_sha256':sources,
             'threads':args.threads,'thread_env':{k:env[k] for k in env if k.startswith(('OMP_','MKL_','OPENBLAS_'))},
             'runs':{}}
    def save():
        (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    save()
    for name,prefix,pheno in (
        ('5k','tests/data/gwas_speed_n5000_m50000/speed','tests/data/gwas_speed_n5000_m50000/speed.pheno'),
        ('20k','benchmarks/r18_1_data/r18_1_20k_1m_bed','benchmarks/r18_1_data/r18_1_20k_1m.pheno')):
        stem=out/('mf_'+name)
        cmd=[str(binary),'--bfile',prefix,'--pheno',pheno,'--fastlmm','--fastlmm-matrix-free',
             '--fastlmm-max-rank','256','--fastlmm-no-fallback','--threads',str(args.threads),'--out',str(stem)]
        started=time.time()
        print(f'START {name}: '+ ' '.join(cmd),flush=True)
        with open(str(stem)+'.log','w') as log:
            rc=subprocess.run(['/usr/bin/time','-v','-o',str(stem)+'.time']+cmd,
                              stdout=log,stderr=subprocess.STDOUT,env=env).returncode
        timing=Path(str(stem)+'.time').read_text()
        log=Path(str(stem)+'.log').read_text()
        run={'command':cmd,'exit_code':rc,'driver_wall_seconds':time.time()-started,
             'time_v':timing,
             'forbidden_log_matches':re.findall(r'global GRM build|Eigen Decomposition of projected G|EigenVce|Calculating Global[^\n]*GRM',log,re.I)}
        summary['runs'][name]=run
        save()
        if rc: raise RuntimeError(f'{name} failed: {rc}')
        manifest=json.loads(Path(str(stem)+'.run_manifest.json').read_text())
        run['performance']=manifest['performance']
        with open(str(stem)+'.low_rank.meta.tsv') as f:
            run['low_rank']=list(csv.DictReader(f,delimiter='\t'))
        reference=Path(args.baseline)/('mf_'+name+'.assoc.txt')
        candidate=Path(str(stem)+'.assoc.txt')
        run['reference']=str(reference)
        run['reference_sha256']=sha256(reference)
        run['candidate_sha256']=sha256(candidate)
        run['accuracy_vs_old_matrix_free']=compare(reference,candidate)
        save()
        print(f'DONE {name}: '+json.dumps(run['performance'])+' '+json.dumps(run['accuracy_vs_old_matrix_free']),flush=True)
    summary['sources_unchanged_during_run']={p:sha256(p)==h for p,h in sources.items()}
    save()


if __name__=='__main__':
    main()
