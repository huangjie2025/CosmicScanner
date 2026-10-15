#!/usr/bin/env python3
"""Build and verify a native macOS package; requires CMake, Eigen and libomp."""
from pathlib import Path
import argparse,hashlib,os,platform,subprocess,sys,tarfile,urllib.request,shutil
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--work',type=Path,required=True)
    ap.add_argument('--eigen',type=Path,required=True)
    ap.add_argument('--libomp-prefix',type=Path,required=True)
    args=ap.parse_args()
    if platform.system()!='Darwin':raise SystemExit('Requires a native macOS builder.')
    root=Path(__file__).resolve().parents[1]
    if sys.version_info < (3,12):raise SystemExit('Use Python 3.12 or newer for safe source extraction.')
    if platform.machine().lower() not in ('arm64','x86_64'):raise SystemExit('Unsupported macOS architecture')
    work=args.work.resolve();work.mkdir(parents=True,exist_ok=True)
    prefix=work/'deps'
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    env['MACOSX_DEPLOYMENT_TARGET']='15.0'
    specs=[
      ('zlib','1.3.2','https://github.com/madler/zlib/archive/refs/tags/v1.3.2.tar.gz','b99a0b86c0ba9360ec7e78c4f1e43b1cbdf1e6936c8fa0f6835c0cd694a495a1'),
      ('zstd','1.5.7','https://github.com/facebook/zstd/archive/refs/tags/v1.5.7.tar.gz','37d7284556b20954e56e1ca85b80226768902e2edabd3b649e9e72c0c9012ee3')]
    for name,version,url,digest in specs:
        archive=work/(name+'.tar.gz')
        if not archive.exists():urllib.request.urlretrieve(url,archive)
        if hashlib.sha256(archive.read_bytes()).hexdigest()!=digest:raise RuntimeError('Dependency checksum mismatch: '+name)
        source=work/(name+'-'+version)
        if not source.exists():
            with tarfile.open(archive) as t:t.extractall(work,filter='data')
        cmake_source=source if name=='zlib' else source/'build/cmake'
        flags=['-DZLIB_BUILD_SHARED=OFF','-DZLIB_BUILD_STATIC=ON','-DZLIB_BUILD_TESTING=OFF'] if name=='zlib' else ['-DZSTD_BUILD_SHARED=OFF','-DZSTD_BUILD_STATIC=ON','-DZSTD_BUILD_TESTS=OFF','-DZSTD_BUILD_PROGRAMS=OFF']
        subprocess.run(['cmake','-S',str(cmake_source),'-B',str(work/('build-'+name)),'-DCMAKE_BUILD_TYPE=Release','-DCMAKE_INSTALL_PREFIX='+str(prefix),*flags],env=env,check=True)
        subprocess.run(['cmake','--build',str(work/('build-'+name)),'-j2'],env=env,check=True)
        subprocess.run(['cmake','--install',str(work/('build-'+name))],env=env,check=True)
    legal=work/'third-party-licenses';legal.mkdir(exist_ok=True)
    for name,source in [('zlib',work/'zlib-1.3.2/LICENSE'),('zstd',work/'zstd-1.5.7/LICENSE')]:
        if not source.exists():raise RuntimeError('Missing '+name+' copyright file')
        shutil.copy2(source,legal/(name+'.txt'))
    for name,license_prefix in [('eigen',args.eigen.resolve().parents[1]),('libomp',args.libomp_prefix.resolve())]:
        candidates=list(license_prefix.glob('*COPYING*'))+list(license_prefix.glob('*LICENSE*'))+list(license_prefix.glob('share/**/LICENSE*'))
        if not candidates:raise RuntimeError('Missing '+name+' license; install full dependency with licenses')
        for i,source in enumerate(candidates):
            if source.is_file():shutil.copy2(source,legal/(name+'-'+str(i)+'-'+source.name))
    if not (prefix/'lib/libz.a').is_file():raise RuntimeError('Static zlib archive was not installed')
    flags=['-DCMAKE_BUILD_TYPE=Release','-DCOSMICSCANNER_BUILD_TESTS=OFF','-DCMAKE_DISABLE_FIND_PACKAGE_BLAS=TRUE','-DCMAKE_DISABLE_FIND_PACKAGE_LAPACK=TRUE','-DCMAKE_PREFIX_PATH='+str(prefix),'-DEIGEN3_INCLUDE_DIR='+str(args.eigen.resolve()),'-DZLIB_LIBRARY='+str(prefix/'lib/libz.a'),'-DZLIB_INCLUDE_DIR='+str(prefix/'include'),'-DOpenMP_CXX_FLAGS=-Xpreprocessor -fopenmp','-DOpenMP_CXX_LIB_NAMES=omp','-DOpenMP_omp_LIBRARY='+str(args.libomp_prefix/'lib/libomp.dylib'),'-DCMAKE_CXX_FLAGS=-I'+str(args.libomp_prefix/'include'),'-DCMAKE_OSX_DEPLOYMENT_TARGET=15.0']
    subprocess.run(['cmake','-S',str(root),'-B',str(work/'build'),*flags],env=env,check=True)
    subprocess.run(['cmake','--build',str(work/'build'),'--target','cosmicscanner','-j2'],env=env,check=True)
    subprocess.run([sys.executable,str(root/'tests/run_cli_regressions.py'),'--bin',str(work/'build/cosmicscanner'),'--out',str(work/'cli-checks')],env=env,check=True)
    subprocess.run([sys.executable,str(root/'tests/run_review_regressions.py'),'--bin',str(work/'build/cosmicscanner'),'--out',str(work/'review-checks')],env=env,check=True)
    subprocess.run([sys.executable,str(root/'scripts/package_native.py'),'--binary',str(work/'build/cosmicscanner'),'--out',str(work/'artifacts'),'--runtime-dir',str(args.libomp_prefix/'lib'),'--licenses',str(legal)],env=env,check=True)
    archives=list((work/'artifacts').glob('CosmicScanner-*.tar.gz'))
    if len(archives)!=1:raise RuntimeError('Expected exactly one macOS archive')
    unpack=work/'archive-roundtrip';unpack.mkdir()
    with tarfile.open(archives[0]) as t:t.extractall(unpack,filter='data')
    package=next(unpack.iterdir())
    # Verify file integrity after transport, then run without Homebrew in PATH.
    for line in (package/'SHA256SUMS').read_text().splitlines():
        digest,relative=line.split('  ',1)
        if hashlib.sha256((package/relative).read_bytes()).hexdigest()!=digest:
            raise RuntimeError('Archive integrity mismatch: '+relative)
    clean=dict(env,PATH='/usr/bin:/bin',COSMICSCANNER_TEST_TMPDIR=str(work))
    for key in list(clean):
        if key.startswith('DYLD_'):clean.pop(key)
    subprocess.run([sys.executable,str(package/'test/verify_native_package.py')],env=clean,cwd=package,check=True)
    import json
    report={'status':'PASS','arch':platform.machine(),'host':platform.platform(),
            'minimum_macos':'15.0','archive':archives[0].name,
            'archive_sha256':hashlib.sha256(archives[0].read_bytes()).hexdigest(),
            'checks':['cli regressions','review regressions','packaged examples','archive hashes','extracted examples'],
            'signing':'ad-hoc; no Developer ID notarization'}
    (work/'artifacts/MACOS_ACCEPTANCE.json').write_text(json.dumps(report,indent=2))
if __name__=='__main__':main()
