#!/usr/bin/env python3
"""Build, regress, package, extract and verify one native CI snapshot."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tarfile
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
COMPRESSION = [
    ('zlib', '1.3.2', 'https://github.com/madler/zlib/archive/refs/tags/v1.3.2.tar.gz',
     'b99a0b86c0ba9360ec7e78c4f1e43b1cbdf1e6936c8fa0f6835c0cd694a495a1'),
    ('zstd', '1.5.7', 'https://github.com/facebook/zstd/archive/refs/tags/v1.5.7.tar.gz',
     '37d7284556b20954e56e1ca85b80226768902e2edabd3b649e9e72c0c9012ee3'),
]


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def capture(command):
    if command[0] == 'git' and platform.system() == 'Windows':
        # Checkout uses Git for Windows; MSYS2 Git has different global EOL
        # and executable-bit defaults. Compare content with checkout semantics.
        command = ['git', '-c', 'core.autocrlf=true', '-c', 'core.filemode=false', *command[1:]]
    return subprocess.check_output(list(map(str, command)), cwd=ROOT, text=True,
                                   encoding='utf-8', errors='replace').strip()


def run(command, env, cwd=ROOT):
    print('+', ' '.join(map(str, command)), flush=True)
    subprocess.run(list(map(str, command)), cwd=cwd, env=env, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--eigen', type=Path, required=True)
    parser.add_argument('--libomp-prefix', type=Path)
    parser.add_argument('--runtime-dir', type=Path)
    args = parser.parse_args()
    if sys.version_info < (3, 12):
        raise SystemExit('Python 3.12+ required for safe archive extraction')
    system = platform.system()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    artifacts = work/'artifacts'
    artifacts.mkdir(exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    env.pop('CONDA_PREFIX', None)
    commit = capture(['git', 'rev-parse', 'HEAD'])
    if os.environ.get('GITHUB_SHA', commit) != commit:
        raise RuntimeError('Checkout does not match the triggering source commit')
    dirty = capture(['git', 'status', '--porcelain', '--untracked-files=no'])
    if dirty:
        raise RuntimeError('Tracked checkout is modified before building:\n'+dirty)
    inventory = {'Windows': ['pacman', '-Q'], 'Darwin': ['brew', 'list', '--versions'],
                 'Linux': ['dpkg-query', '-W']}[system]
    identity = {
        'schema': 1, 'source_commit': commit,
        'source_tree': capture(['git', 'show', '-s', '--format=%T', 'HEAD']),
        'source_repository': os.environ.get('GITHUB_REPOSITORY', 'local'),
        'source_ref': os.environ.get('GITHUB_REF', 'local'),
        'workflow': os.environ.get('GITHUB_WORKFLOW_REF', 'local'),
        'run_url': (os.environ.get('GITHUB_SERVER_URL', 'https://github.com') + '/' +
                    os.environ.get('GITHUB_REPOSITORY', '') + '/actions/runs/' +
                    os.environ.get('GITHUB_RUN_ID', 'local')),
        'run_attempt': os.environ.get('GITHUB_RUN_ATTEMPT', '1'),
        'built_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'platform': platform.platform(), 'architecture': platform.machine(),
        'runner_image': {k: os.environ.get(k) for k in ['ImageOS', 'ImageVersion', 'RUNNER_OS', 'RUNNER_ARCH']},
        'python': sys.version, 'cmake': capture(['cmake', '--version']),
        'compiler': capture(['c++', '--version']), 'installed_packages': capture(inventory),
        'compression_sources': [{'name': n, 'version': v, 'url': u, 'sha256': h}
                                for n, v, u, h in COMPRESSION],
        'baseline': {'Windows': 'Windows Server 2022 / UCRT64',
                     'Linux': 'Ubuntu 24.04 / glibc 2.39+', 'Darwin': 'macOS 15+'}[system],
    }
    prefix = work/'deps'
    legal = work/'licenses'
    legal.mkdir(exist_ok=True)
    for name, version, url, digest in COMPRESSION:
        archive = work/(name+'.tar.gz')
        if not archive.exists():
            urllib.request.urlretrieve(url, archive)
        if sha(archive) != digest:
            raise RuntimeError('Dependency checksum mismatch: '+name)
        with tarfile.open(archive) as stream:
            stream.extractall(work, filter='data')
        source = work/(name+'-'+version)
        shutil.copy2(source/'LICENSE', legal/(name+'.txt'))
        flags = (['-DZLIB_BUILD_SHARED=OFF', '-DZLIB_BUILD_STATIC=ON', '-DZLIB_BUILD_TESTING=OFF']
                 if name == 'zlib' else ['-DZSTD_BUILD_SHARED=OFF', '-DZSTD_BUILD_STATIC=ON',
                                        '-DZSTD_BUILD_TESTS=OFF', '-DZSTD_BUILD_PROGRAMS=OFF'])
        run(['cmake', '-S', source if name == 'zlib' else source/'build/cmake',
             '-B', work/('build-'+name), '-G', 'Ninja', '-DCMAKE_BUILD_TYPE=Release',
             '-DCMAKE_INSTALL_PREFIX='+str(prefix), *flags], env)
        run(['cmake', '--build', work/('build-'+name), '--parallel', '2'], env)
        run(['cmake', '--install', work/('build-'+name)], env)
    # Preserve the license texts for compiler runtimes, OpenMP and Eigen.
    if system == 'Windows':
        license_root = args.runtime_dir.resolve().parent/'share/licenses'
        if not license_root.is_dir():
            raise RuntimeError('Missing MSYS2 dependency license directory')
        shutil.copytree(license_root, legal/'msys2', dirs_exist_ok=True)
    elif system == 'Linux':
        for name in ['libeigen3-dev', 'libgomp1', 'libgcc-s1', 'libstdc++6', 'gcc-13-base']:
            license_root = Path('/usr/share/doc')/name
            if not (license_root/'copyright').is_file():
                raise RuntimeError('Missing dependency copyright: '+name)
            (legal/name).mkdir(exist_ok=True)
            for path in license_root.iterdir():
                if path.is_file() and (path.name.startswith('copyright') or path.name.startswith('NEWS')):
                    shutil.copy2(path, legal/name/path.name)
        for path in Path('/usr/share/common-licenses').iterdir():
            if path.is_file():
                shutil.copy2(path, legal/path.name)
    else:
        for name, location in [('eigen', args.eigen.resolve().parents[1]),
                               ('libomp', args.libomp_prefix.resolve())]:
            paths = [*location.glob('*LICENSE*'), *location.glob('*COPYING*'),
                     *location.glob('share/**/LICENSE*')]
            paths = [p for p in paths if p.is_file()]
            if not paths:
                raise RuntimeError('Missing dependency license: '+name)
            for index, path in enumerate(paths):
                shutil.copy2(path, legal/(name+'-'+str(index)+'-'+path.name))
    build = work/'build'
    flags = ['-DCMAKE_BUILD_TYPE=Release', '-DCMAKE_CXX_COMPILER_LAUNCHER=ccache', '-DCOSMICSCANNER_BUILD_TESTS=OFF',
             '-DCMAKE_DISABLE_FIND_PACKAGE_BLAS=TRUE', '-DCMAKE_DISABLE_FIND_PACKAGE_LAPACK=TRUE',
             '-DCMAKE_PREFIX_PATH='+str(prefix), '-DEIGEN3_INCLUDE_DIR='+str(args.eigen.resolve()),
             '-DZLIB_USE_STATIC_LIBS=ON', '-DZLIB_INCLUDE_DIR='+str(prefix/'include')]
    # Avoid platform-dependent find-library ordering (zlib names differ on MinGW).
    zlib = list((prefix/'lib').glob('*z*.a'))
    zlib = [p for p in zlib if 'zstd' not in p.name]
    if len(zlib) != 1:
        raise RuntimeError('Expected exactly one static zlib archive: '+str(zlib))
    flags.append('-DZLIB_LIBRARY='+str(zlib[0]))
    if system == 'Darwin':
        omp = args.libomp_prefix.resolve()
        flags += ['-DOpenMP_CXX_FLAGS=-Xpreprocessor -fopenmp', '-DOpenMP_CXX_LIB_NAMES=omp',
                  '-DOpenMP_omp_LIBRARY='+str(omp/'lib/libomp.dylib'),
                  '-DCMAKE_CXX_FLAGS=-I'+str(omp/'include'), '-DCMAKE_OSX_DEPLOYMENT_TARGET=15.0']
    elif system == 'Linux':
        flags.append('-DCMAKE_EXE_LINKER_FLAGS=-static-libgcc -static-libstdc++')
    identity['configure_flags'] = flags
    provenance = artifacts/'BUILD_ENVIRONMENT.json'
    provenance.write_text(json.dumps(identity, indent=2)+'\n', encoding='utf-8')
    run(['cmake', '-S', ROOT, '-B', build, '-G', 'Ninja', *flags], env)
    run(['cmake', '--build', build, '--target', 'cosmicscanner', '--parallel', '2'], env)
    # Guard the record-field initialization regression independently of CLI tests.
    record_test = work/('test_pheno_record_fields.exe' if system == 'Windows' else 'test_pheno_record_fields')
    run(['c++', '-std=c++17', '-I'+str(args.eigen.resolve()), '-I'+str(ROOT/'libCosmicCore/src'),
         ROOT/'tests/test_pheno_record_fields.cpp', ROOT/'libCosmicCore/src/pheno_reader.cpp',
         '-o', record_test], env)
    run([record_test, work/'pheno-records.txt'], env)
    binary = build/('cosmicscanner.exe' if system == 'Windows' else 'cosmicscanner')
    for suite in ['cli', 'review']:
        run([sys.executable, ROOT/('tests/run_'+suite+'_regressions.py'), '--bin', binary,
             '--out', work/(suite+'-checks')], env)
    version = '2.0.0-dev.'+commit[:12]+'.'+os.environ.get('GITHUB_RUN_NUMBER', 'local')
    command = [sys.executable, ROOT/'scripts/package_native.py', '--binary', binary,
               '--out', artifacts, '--licenses', legal, '--version', version,
               '--provenance', provenance]
    if args.runtime_dir:
        command += ['--runtime-dir', args.runtime_dir.resolve()]
    if system == 'Darwin':
        command += ['--runtime-dir', args.libomp_prefix.resolve()/'lib']
    run(command, env)
    archives = [*artifacts.glob('*.zip'), *artifacts.glob('*.tar.gz')]
    if len(archives) != 1:
        raise RuntimeError('Expected exactly one native archive')
    archive = archives[0]
    unpack = work/'extracted'
    unpack.mkdir()
    if archive.suffix == '.zip':
        with zipfile.ZipFile(archive) as stream:
            stream.extractall(unpack)
    else:
        with tarfile.open(archive) as stream:
            stream.extractall(unpack, filter='data')
    package = next(unpack.iterdir())
    for line in (package/'SHA256SUMS').read_text().splitlines():
        digest, relative = line.split('  ', 1)
        if sha(package/relative) != digest:
            raise RuntimeError('Extracted file checksum mismatch: '+relative)
    package_id = json.loads((package/'PACKAGE_ID.json').read_text())
    if package_id['build']['source_commit'] != commit:
        raise RuntimeError('Package provenance differs from source commit')
    clean = dict(env, PATH=(str(Path(os.environ['SystemRoot'])/'System32')
                           if system == 'Windows' else '/usr/bin:/bin'),
                 COSMICSCANNER_TEST_TMPDIR=str(work))
    for key in list(clean):
        if key.startswith(('LD_', 'DYLD_')):
            clean.pop(key)
    run([sys.executable, package/'test/verify_native_package.py'], clean, cwd=package)
    report = {'status': 'PASS', 'source_commit': commit, 'run_url': identity['run_url'],
              'archive': archive.name, 'archive_sha256': sha(archive),
              'packaged_binary_sha256': package_id['binary_sha256'],
              'checks': ['phenotype record fields', 'CLI regressions', 'review regressions',
                         'packaged LMM/multivariate/repeatability/GRM',
                         'archive roundtrip file hashes', 'extracted package clean-PATH tests']}
    (artifacts/'CI_RESULT.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    if os.environ.get('GITHUB_STEP_SUMMARY'):
        with open(os.environ['GITHUB_STEP_SUMMARY'], 'a', encoding='utf-8') as stream:
            stream.write('### '+system+' '+platform.machine()+' — PASS\n\n'
                         'Source: `'+commit+'`\n\nPackage: `'+archive.name+'`\n\n'
                         'SHA-256: `'+sha(archive)+'`\n\n'
                         'Development snapshot; model support remains scoped and RRM experimental.\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
