from pathlib import Path
import argparse,hashlib,json,os,platform,shutil,subprocess,sys,tarfile
ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--source-run',required=True);args=ap.parse_args()
assert platform.system()=='Darwin'
root=Path(__file__).resolve().parents[1];out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
source_report=json.loads((args.input/'artifacts/MACOS_ACCEPTANCE.json').read_text());assert source_report['status']=='PASS'
archives=list((args.input/'artifacts').glob('CosmicScanner-*.tar.gz'));assert len(archives)==1
unpack=out/'staging';unpack.mkdir()
with tarfile.open(archives[0]) as t:t.extractall(unpack,filter='data')
pkg=next(unpack.iterdir());identity=json.loads((pkg/'PACKAGE_ID.json').read_text());exe=pkg/'bin/cosmicscanner'
assert identity['binary_sha256']==hashlib.sha256(exe.read_bytes()).hexdigest()
for name in ['USER_MANUAL_CN.md','USER_MANUAL_EN.md','NATIVE_INSTALL_CN.md','CLI_GUIDE_CN.md','CLI_GUIDE_EN.md','CosmicScanner_v2.0.0_User_Manual_CN.docx']:
    shutil.copy2(root/'docs'/name,pkg/'docs'/name)
identity['documentation_commit']=os.environ.get('GITHUB_SHA');identity['source_run']=args.source_run
(pkg/'PACKAGE_ID.json').write_text(json.dumps(identity,indent=2))
(pkg/'README.txt').write_text('CosmicScanner v2.0.0 macOS acceptance package\nmacOS 15.0 or later. Choose arm64 for Apple Silicon, x86_64 for Intel.\nUnpack the entire folder; keep bin/ and lib/ together.\nRun from the package folder: ./bin/cosmicscanner --help\nSee docs/NATIVE_INSTALL_CN.md for a complete runnable example.\nNo Python, Homebrew or compiler is required to run the software.\nPython 3 is only required for the optional test/verify_native_package.py.\nThis package passed the recorded macOS CLI and example tests.\nThis does not certify every scientific model as stable.\nAd-hoc signed; not Apple Developer ID signed or notarized.\n')
# Verify loader paths and signatures before producing the final archive.
loader=[]
for binary in [exe,*sorted((pkg/'lib').glob('*.dylib'))]:
    text=subprocess.check_output(['otool','-L',str(binary)],text=True);loader.append(text)
    for line in text.splitlines()[1:]:
        dep=line.strip().split(' (',1)[0]
        assert dep.startswith(('/usr/lib/','/System/Library/','@loader_path/','@rpath/')),(binary,dep)
    subprocess.run(['codesign','--verify','--strict',str(binary)],check=True)
(out/'runtime-dependencies.txt').write_text('\n'.join(loader))
sums={str(p.relative_to(pkg)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(pkg.rglob('*')) if p.is_file() and p.name!='SHA256SUMS'}
(pkg/'SHA256SUMS').write_text(''.join(h+'  '+p+'\n' for p,h in sums.items()))
archive=out/(pkg.name+'.tar.gz')
with tarfile.open(archive,'w:gz') as t:t.add(pkg,arcname=pkg.name)
digest=hashlib.sha256(archive.read_bytes()).hexdigest();(out/(archive.name+'.sha256')).write_text(digest+'  '+archive.name+'\n')
roundtrip=out/'roundtrip';roundtrip.mkdir()
with tarfile.open(archive) as t:t.extractall(roundtrip,filter='data')
final=next(roundtrip.iterdir())
for line in (final/'SHA256SUMS').read_text().splitlines():
    expected,relative=line.split('  ',1);assert hashlib.sha256((final/relative).read_bytes()).hexdigest()==expected
clean=dict(os.environ,PATH='/usr/bin:/bin',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',COSMICSCANNER_TEST_TMPDIR=str(out))
for key in list(clean):
    if key.startswith('DYLD_'):clean.pop(key)
subprocess.run([sys.executable,str(final/'test/verify_native_package.py')],cwd=final,env=clean,check=True)
report=dict(source_report,source_run=args.source_run,documentation_commit=os.environ.get('GITHUB_SHA'),archive_sha256=digest,final_archive=archive.name,final_archive_test='PASS',runtime_dependency_check='PASS',signature_check='PASS')
(out/'MACOS_ACCEPTANCE.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))

