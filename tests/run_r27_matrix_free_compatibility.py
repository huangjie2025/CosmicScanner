#!/usr/bin/env python3
"""Reuse the R19-E 5k command and comparator for before/after compatibility.

Runs sequentially; does not extrapolate or rerun the 20k/50k benchmarks.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import subprocess

from run_matrix_free_kernel_benchmark import compare, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', type=Path, required=True)
    parser.add_argument('--after', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    prefix = root / 'tests/data/gwas_speed_n5000_m50000/speed'
    env = dict(os.environ, OMP_NUM_THREADS='32', MKL_NUM_THREADS='32',
               OPENBLAS_NUM_THREADS='32', OMP_DYNAMIC='FALSE', MKL_DYNAMIC='FALSE',
               MKL_INTERFACE_LAYER='LP64,GNU')
    report = {'scope': 'R19-E 5k x 50k, same LOCO command, before/after', 'runs': {}}
    metadata = {}
    for label, binary in (('before', args.before.resolve()), ('after', args.after.resolve())):
        stem = out / label
        command = [str(binary), '--bfile', str(prefix), '--pheno', str(prefix) + '.pheno',
                   '--fastlmm', '--fastlmm-matrix-free', '--fastlmm-max-rank', '256',
                   '--fastlmm-no-fallback', '--threads', '32', '--out', str(stem)]
        with Path(str(stem) + '.log').open('w') as log:
            subprocess.run(['/usr/bin/time', '-v', '-o', str(stem) + '.time', *command],
                           cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        text = Path(str(stem) + '.log').read_text()
        assert 'dense GRM not built' in text
        assert not re.search(r'global GRM build|Eigen Decomposition of projected G|EigenVce|Calculating Global[^\n]*GRM', text, re.I)
        with Path(str(stem) + '.low_rank.meta.tsv').open() as stream:
            metadata[label] = list(csv.DictReader(stream, delimiter='\t'))
        assert metadata[label] and all(row['FALLBACK'] == 'none' for row in metadata[label])
        report['runs'][label] = {'argv': command, 'binary_sha256': sha256(binary),
                                  'association_sha256': sha256(str(stem) + '.assoc.txt'),
                                  'chromosome_nulls': len(metadata[label])}
    report['comparison'] = compare(str(out / 'before.assoc.txt'), str(out / 'after.assoc.txt'))
    report['bitwise_identical'] = (out / 'before.assoc.txt').read_bytes() == (out / 'after.assoc.txt').read_bytes()
    fields = ('SCOPE', 'ALGORITHM', 'EFFECTIVE_RANK', 'SKETCH_SIZE', 'TRACE_FRACTION',
              'RESIDUAL_EIGENVALUE', 'RANK_SELECTION', 'FALLBACK')
    report['rank_metadata_identical'] = ([tuple(row[field] for field in fields) for row in metadata['before']] ==
                                       [tuple(row[field] for field in fields) for row in metadata['after']])
    (out / 'compatibility_summary.json').write_text(json.dumps(report, indent=2) + '\n')
    assert report['bitwise_identical'] and report['rank_metadata_identical'] and report['comparison']['accuracy_pass']
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
