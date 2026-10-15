#!/usr/bin/env python3
"""CLI compatibility, explicit task/model validation, and documented examples."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin', type=Path, default=root/'build/cosmicscanner')
    parser.add_argument('--out', type=Path)
    parser.add_argument('--baseline', type=Path)
    args = parser.parse_args()
    binary = args.bin.resolve()
    out = args.out.resolve() if args.out else Path(tempfile.mkdtemp(prefix='cosmicscanner-cli-test-'))
    out.mkdir(parents=True, exist_ok=True)
    data = root/'examples/data/r33_synthetic/dataset'
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    results = []

    def run(name, flags, ok=True, expected=None, executable=None):
        command = [str(executable or binary), *map(str, flags)]
        p = subprocess.run(command, cwd=out, env=env, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=90)
        log = p.stdout+'\n'+p.stderr
        (out/(name+'.log')).write_text(log, encoding="utf-8")
        assert (p.returncode == 0) == ok, (name, p.returncode, log[-2500:])
        if expected:
            assert expected in log, (name, expected, log[-2500:])
        results.append({'name':name,'argv':command,'exit':p.returncode})
        return p.stdout

    brief = run('help', ['--help'])
    assert len(brief.splitlines()) <= 60 and '--model NAME' in brief
    full = run('help_all', ['--help-all'])
    assert len(full.splitlines()) > 300
    for topic in ['gwas','models','input','grm','repeatability','rrm','ssgwas','multivariate']:
        assert run('help_'+topic, ['help',topic]) == run('flag_help_'+topic, ['--help',topic])
    run('invalid_help', ['help','invalid'], False, 'Unknown help topic')
    run('version', ['version'], expected='CosmicScanner v')
    run('extra_version', ['version','ignored'], False)
    run('unknown_task', ['unknown'], False, 'Unknown command')
    run('missing_model', ['gwas','--model'], False, '--model requires a value')
    run('unknown_model', ['--model','not-a-model'], False, 'Unknown model')
    run('duplicate_model', ['--model','lm','--model','lmm'], False, 'only once')
    run('model_conflict', ['gwas','--model','lmm','--multivariate'], False, 'Conflicting model')
    run('task_conflict', ['grm','--model','lm'], False, 'another task')
    run('gwas_qc_conflict', ['gwas','--geno','.1'], False, 'standalone')
    run('missing_traits', ['gwas','--model','multivariate'], False, 'requires --mpheno-list')
    run('pedigree_context', ['--pedigree','dummy'], False, 'requires an animal')
    run('ped_alias_conflict', ['--pedigree','a','--ped-genotypes','b'], False, 'not both')
    run('ped_repeated', ['--pedigree','a','--ped','b'], False, 'only once')

    common = ['--bfile',data/'synth','--threads','1']
    single = common+['--pheno',data/'quantitative.pheno','--no-loco']
    cases = [
        ('lmm',single+['--gwas'],single+['--model','lmm'],'.assoc.txt'),
        ('lm',single+['--lm'],single+['--model','lm'],'.lm.assoc.txt'),
        ('mv',common+['--pheno',data/'multivariate.pheno','--mpheno-list','1,2,3','--no-loco',
                      '--reml-maxit','1200','--reml-tol','1e-10','--multivariate','--gwas'],
              common+['--pheno',data/'multivariate.pheno','--mpheno-list','1,2,3','--no-loco',
                      '--reml-maxit','1200','--reml-tol','1e-10','--model','multivariate'],'.assoc.mv.txt'),
        ('rep',common+['--pheno',data/'long.pheno','--tpheno','1','--mpheno','2','--no-loco',
                       '--reml-vg','.5','--reml-vpe','.2','--reml-ve','.3','--repeatability-mode','exact',
                       '--repeatability-model','--gwas'],
               common+['--pheno',data/'long.pheno','--tpheno','1','--mpheno','2','--no-loco',
                       '--reml-vg','.5','--reml-vpe','.2','--reml-ve','.3','--repeatability-mode','exact',
                       '--model','repeatability'],'.assoc.txt')]
    for name, old, new, suffix in cases:
        run(name+'_old',old+['--out',out/(name+'_old')])
        run(name+'_new',['gwas',*new,'--out',out/(name+'_new')])
        assert (out/(name+'_old'+suffix)).read_bytes() == (out/(name+'_new'+suffix)).read_bytes(), name
        if args.baseline:
            run(name+'_baseline',old+['--out',out/(name+'_baseline')],executable=args.baseline.resolve())
            assert (out/(name+'_old'+suffix)).read_bytes() == (out/(name+'_baseline'+suffix)).read_bytes()
    run('flat_model',single+['--gwas','--model','lm','--out',out/'flat_model'])
    assert (out/'flat_model.lm.assoc.txt').read_bytes() == (out/'lm_old.lm.assoc.txt').read_bytes()
    run('grm_old',common+['--make-grm','--out',out/'grm_old'])
    run('grm_new',['grm',*common,'--out',out/'grm_new'])
    for suffix in ['.grm.bin','.grm.id']:
        assert (out/('grm_old'+suffix)).read_bytes() == (out/('grm_new'+suffix)).read_bytes()
    # Historical --model lmm was a modifier, not a task/model-family reset.
    legacy_modifier_cases = [
        ('modifier_mv', cases[2][1]+['--model','lmm'], 'mv_old', '.assoc.mv.txt'),
        ('modifier_grm', common+['--make-grm','--model','lmm'], 'grm_old', '.grm.bin')]
    for name, flags, reference, suffix in legacy_modifier_cases:
        run(name,flags+['--out',out/name])
        assert (out/(name+suffix)).read_bytes() == (out/(reference+suffix)).read_bytes()
        if args.baseline:
            run(name+'_baseline',flags+['--out',out/(name+'_baseline')],executable=args.baseline.resolve())
            assert (out/(name+suffix)).read_bytes() == (out/(name+'_baseline'+suffix)).read_bytes()
    run('native_semi_rejected',['gwas',*common,'--pheno',data/'multivariate.pheno',
        '--mpheno-list','1,2,3','--model','multivariate','--scan-tier','semi'],False,'not implemented')
    pedigree = out/'founders.ped'
    pedigree.write_text(''.join(row.split()[1]+' 0 0\n' for row in (data/'synth.fam').read_text().splitlines()))
    run('animal_pedigree',['gwas',*single,'--model','animal','--pedigree',pedigree,
                          '--reml-vg','.5','--reml-ve','.5','--out',out/'animal'])
    # The typed command, rather than frontend-injected flags, remains provenance.
    manifest=(out/'lmm_new.run_manifest.json').read_text()
    assert '"gwas"' in manifest and '"--model"' in manifest
    result={'status':'PASS','commands':len(results),'output':str(out),'cases':results}
    (out/'cli_verification.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({'status':'PASS','commands':len(results),'output':str(out)}))


if __name__ == '__main__':
    main()
