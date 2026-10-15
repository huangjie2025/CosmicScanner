#!/usr/bin/env python3
"""Targeted regressions for the 2026-10-13 review; all outputs are isolated.

Uses the existing tiny R33 synthetic fixture. Does not rebuild or write under
the repository unless --out explicitly selects such a directory.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile


def main():
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bin', type=Path, default=root/'build/cosmicscanner')
    ap.add_argument('--data', type=Path, default=root/'examples/data/r33_synthetic/dataset')
    ap.add_argument('--out', type=Path)
    args = ap.parse_args()
    binary = args.bin.resolve()
    data = args.data.resolve()
    out = args.out.resolve() if args.out else Path(tempfile.mkdtemp(prefix='cosmicscanner-review-regression-'))
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    records = []

    def run(name, flags, success=True, message=None):
        argv = [str(binary), '--threads', '1', *map(str, flags), '--out', str(out/name)]
        p = subprocess.run(argv, cwd=out, env=env, text=True, capture_output=True, encoding="utf-8", errors="replace", timeout=90)
        log = p.stdout+'\n'+p.stderr
        (out/(name+'.log')).write_text(log, encoding="utf-8")
        records.append({'case': name, 'argv': argv, 'exit': p.returncode})
        assert (p.returncode == 0) == success, (name, p.returncode, log[-5000:])
        if message:
            assert message in log, (name, message, log[-5000:])
        print('PASS:', name, flush=True)
        return log

    def text(path):
        raw = path.read_bytes()
        return (gzip.decompress(raw) if raw.startswith(b'\x1f\x8b') else raw).decode()

    def rows(name, suffix):
        return list(csv.DictReader(io.StringIO(text(out/(name+suffix))), delimiter='\t'))

    def summary_check(name):
        assoc = rows(name, '.assoc.txt')
        content = text(out/(name+'.cosmic.sumstats'))
        lines = [s for s in content.splitlines() if s and not s.startswith('#')]
        assert sum(s.startswith('CHR\tSNP\tPOS\t') for s in lines) == 1
        summary = list(csv.DictReader(lines, delimiter='\t'))
        assert len(assoc) == len(summary) == 180
        assert {r['SNP'] for r in assoc} == {r['SNP'] for r in summary}
        assert {r['CHR'] for r in summary} == {'1', '2'}

    basic = ['--bfile', data/'synth', '--pheno', data/'quantitative.pheno', '--gwas']
    for name, extra in [('loco_plain', []), ('loco_gzip', ['--out-gz'])]:
        run(name, basic+extra)
        summary_check(name)
        # Reusing the output prefix starts a fresh run, not a duplicate append.
        run(name, basic+extra)
        summary_check(name)

    extract = out/'extract40.txt'
    extract.write_text('\n'.join(s.split()[1] for s in (data/'synth.bim').read_text().splitlines()[:40])+'\n')
    rep = ['--bfile', data/'synth', '--pheno', data/'long.pheno', '--tpheno', '1', '--mpheno', '2',
           '--repeatability-model', '--repeatability-mode', 'exact', '--gwas', '--no-loco',
           '--reml-vg', '.5', '--reml-vpe', '.2', '--reml-ve', '.3']
    shared = ['--repeatability-grm-cache', out/'shared']
    run('cache_seed', rep+shared)
    run('cache_changed', rep+shared+['--extract', extract], message='identity/shape mismatch')
    run('cache_fresh', rep+['--repeatability-grm-cache', out/'fresh', '--extract', extract])
    a, b = rows('cache_changed', '.assoc.txt'), rows('cache_fresh', '.assoc.txt')
    assert a == b and len(a) == 40
    run('cache_reuse', rep+shared+['--extract', extract], message='loaded GRM cache')
    c = rows('cache_reuse', '.assoc.txt')
    assert len(c) == len(a)
    for ra, rc in zip(a, c):
        assert ra['SNP'] == rc['SNP']
        for col in ['BETA', 'SE', 'P', 'LOG10P']:
            assert abs(float(ra[col])-float(rc[col])) < 2e-6, (col, ra, rc)
    # Old metadata without identity must be invalidated, even with a valid binary.
    (out/'shared.global.meta.txt').write_text('scale 100\nchroms 1 2\n')
    run('cache_legacy', rep+shared+['--extract', extract], message='identity/shape mismatch')
    assert rows('cache_legacy', '.assoc.txt') == b

    mv = ['--bfile', data/'synth', '--pheno', data/'multivariate.pheno', '--mpheno-list', '1,2,3',
          '--multivariate', '--gwas', '--no-loco', '--reml-maxit', '1200', '--reml-tol', '1e-10']
    run('mv_native_semi', mv+['--scan-tier', 'semi'], success=False, message='semi is not implemented')
    run('mv_approx', mv+['--scan-tier', 'approx'])
    run('mv_legacy', mv+['--mv-method', 'emmv'])
    assert rows('mv_approx', '.assoc.mv.txt') == rows('mv_legacy', '.assoc.mv.txt')
    meta = dict(csv.reader((out/'mv_legacy.scan.meta.txt').read_text().splitlines(), delimiter='\t'))
    assert meta['resolved_scan_tier'] == 'approx' and meta['snp_vc_strategy'] == 'fixed_null_covariance'
    run('mv_compatible_semi', mv+['--scan-tier', 'semi', '--mv-compatible', 'emmvgwas'])
    assert len(rows('mv_compatible_semi', '.assoc.mv.txt')) == 180

    low_budget = mv[:-4]+['--reml-maxit', '1', '--reml-tol', '1e-12', '--mv-method', 'wald']
    run('mv_unconverged', low_budget, success=False, message='null model did not converge')
    failed_table = out/'mv_unconverged.assoc.mv.txt'
    assert not failed_table.exists() or not rows('mv_unconverged', '.assoc.mv.txt')
    run('batch_fail', ['--bfile', data/'synth', '--pheno', out/'absent.pheno', '--batch-models', 'lm,lmm'],
        success=False, message='2/2 models failed')
    run('batch_unknown', basic+['--batch-models', 'lm,unknown'], success=False, message='1/2 models failed')
    run('batch_success', basic+['--batch-models', 'lm,lmm'], message='2 models succeeded')

    result = {'status': 'PASS', 'binary': str(binary), 'binary_sha256': hashlib.sha256(binary.read_bytes()).hexdigest(),
              'output': str(out), 'commands': records}
    (out/'review_regressions.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
