#!/usr/bin/env python3
"""Assemble and smoke-test native packages. Never publishes a release."""
import argparse,csv,hashlib,json,math,os,platform,re,shutil,subprocess,tarfile,tempfile
from pathlib import Path

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--binary',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--runtime-dir',type=Path)
    ap.add_argument('--licenses',type=Path,required=True)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument('--version',default='2.0.0')
    ap.add_argument('--provenance',type=Path,help='Build environment JSON for a CI snapshot')
    args=ap.parse_args()
    system=platform.system();arch=platform.machine().lower()
    if system not in ('Windows','Darwin','Linux'):raise SystemExit('Unsupported native platform.')
    arch={'amd64':'x86_64','aarch64':'arm64'}.get(arch,arch)
    name=f'CosmicScanner-v{args.version}-'+{'Windows':'windows','Darwin':'macos','Linux':'linux'}[system]+'-'+arch
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    pkg=out/name
    if pkg.exists():raise SystemExit('Refusing to overwrite existing package: '+str(pkg))
    for d in ['bin','lib','licenses','docs','examples/dataset','test']:(pkg/d).mkdir(parents=True,exist_ok=True)
    exe=pkg/'bin'/('cosmicscanner.exe' if system=='Windows' else 'cosmicscanner')
    shutil.copy2(args.binary,exe)
    runtime=[]
    if system=='Windows':
        if not args.runtime_dir:raise SystemExit('--runtime-dir must identify the compiler runtime directory.')
        objdump=shutil.which('objdump') or str(args.runtime_dir/'objdump.exe')
        pending=[exe];copied=set()
        # Resolve non-system DLLs transitively from the toolchain; system DLLs
        # must resolve from Windows, never from an unrelated PATH installation.
        while pending:
            p=pending.pop()
            text=subprocess.check_output([objdump,'-p',str(p)],text=True)
            for dll in re.findall(r'DLL Name:\s*(\S+)',text):
                if dll.lower() in copied:continue
                candidate=args.runtime_dir/dll
                if candidate.is_file():
                    dst=pkg/'bin'/dll;shutil.copy2(candidate,dst)
                    copied.add(dll.lower());pending.append(dst);runtime.append(dll)
                elif not (Path(os.environ['SystemRoot'])/'System32'/dll).exists() and not dll.lower().startswith(('api-ms-','ext-ms-')):
                    raise RuntimeError('Unresolved DLL: '+dll)
        for legal in ['COPYING','COPYING.RUNTIME']:
            for candidate in [args.runtime_dir.parent/legal,args.runtime_dir.parent/'share/licenses/gcc'/legal]:
                if candidate.is_file():shutil.copy2(candidate,pkg/'licenses'/('gcc-'+legal));break
    elif system=='Linux':
        # Keep glibc and its loader on the host. Bundle other shared runtimes
        # recursively and resolve them relative to the extracted package.
        pending=[exe];copied=set()
        system_lib=re.compile(r'^(ld-linux.*|lib(c|m|pthread|rt|dl|util|resolv)\.so\..*)$')
        while pending:
            p=pending.pop()
            dependencies=subprocess.check_output(['ldd',str(p)],text=True)
            if 'not found' in dependencies:raise RuntimeError('Unresolved Linux dependency: '+dependencies)
            for lib,location in re.findall(r'^\s*(\S+)\s+=>\s+(/\S+)',dependencies,re.M):
                if system_lib.match(lib) or lib in copied:continue
                dst=pkg/'lib'/lib;shutil.copy2(location,dst)
                copied.add(lib);pending.append(dst);runtime.append(lib)
                subprocess.run(['patchelf','--set-rpath','$ORIGIN',str(dst)],check=True)
        subprocess.run(['patchelf','--set-rpath','$ORIGIN/../lib',str(exe)],check=True)
    else:
        # Homebrew libomp is the only expected non-system shared dependency.
        # zlib/zstd must be linked statically by the build job.
        pending=[exe];copied=set()
        while pending:
            p=pending.pop()
            for line in subprocess.check_output(['otool','-L',str(p)],text=True).splitlines()[1:]:
                dep=line.strip().split(' (',1)[0]
                if dep.startswith(('/usr/lib/','/System/Library/','@loader_path/')):continue
                if dep==str(p):continue
                candidate=Path(dep)
                if dep.startswith('@rpath/') and args.runtime_dir:candidate=args.runtime_dir/Path(dep).name
                if not candidate.is_file():raise RuntimeError('Unresolved macOS dylib: '+dep)
                name=Path(dep).name
                if name not in copied:
                    dst=pkg/'lib'/name;shutil.copy2(candidate,dst);copied.add(name);pending.append(dst);runtime.append(name)
                    subprocess.run(['install_name_tool','-id','@rpath/'+name,str(dst)],check=True)
                replacement=('@loader_path/../lib/' if p==exe else '@loader_path/')+name
                subprocess.run(['install_name_tool','-change',dep,replacement,str(p)],check=True)
        for p in [*list((pkg/'lib').glob('*.dylib')),exe]:
            subprocess.run(['codesign','--force','--sign','-',str(p)],check=True)
    for n in ['synth.bed','synth.bim','synth.fam','quantitative.pheno','multivariate.pheno','long.pheno']:
        shutil.copy2(args.root/'examples/data/r33_synthetic/dataset'/n,pkg/'examples/dataset'/n)
    shutil.copy2(args.root/'LICENSE',pkg/'LICENSE')
    if not args.licenses.is_dir():raise RuntimeError('Missing third-party license directory')
    shutil.copytree(args.licenses,pkg/'licenses/third-party')
    for n in ['CLI_GUIDE_CN.md','CLI_GUIDE_EN.md','NATIVE_INSTALL_CN.md',
              'USER_MANUAL_CN.md','USER_MANUAL_EN.md','SUPPORTED_MODES.md',
              'CosmicScanner_v2.0.0_User_Manual_CN.docx','CosmicScanner_v2.0.0_User_Manual_EN.docx']:
        shutil.copy2(args.root/'docs'/n,pkg/'docs'/n)
    shutil.copy2(args.root/'scripts/verify_native_package.py',pkg/'test/verify_native_package.py')
    (pkg/'README.txt').write_text(
        f'CosmicScanner {args.version}: '+('CI SNAPSHOT' if args.provenance else 'ACCEPTANCE PACKAGE')+', not a stable release.\n'
        'Unpack the whole folder, keeping the executable and libraries together.\n'
        'Windows: .\\bin\\cosmicscanner.exe --help\n'
        'Linux/macOS: ./bin/cosmicscanner --help\n'
        'Self-test (Python 3, only required for verification): python test/verify_native_package.py\n'
        'Detailed source-build examples in docs use different binary/data paths.\n'
        'Do not treat experimental models in help as stable functionality.\n'
        'macOS ad-hoc signature is not Apple Developer ID notarization.\n')
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    if system=='Windows':env['PATH']=str(Path(os.environ['SystemRoot'])/'System32')
    else:
        env['PATH']='/usr/bin:/bin'
        for key in list(env):
            if key.startswith(('DYLD_','LD_')):env.pop(key)
    subprocess.run([str(exe),'--version'],env=env,cwd=pkg,check=True)
    subprocess.run([__import__('sys').executable,str(pkg/'test/verify_native_package.py')],env=env,cwd=pkg,check=True)
    identity={'version':args.version,'status':'CI_SNAPSHOT' if args.provenance else 'ACCEPTANCE_ONLY','platform':system,'arch':arch,
              'host':platform.platform(),'runtimes':runtime,'binary_sha256':hashlib.sha256(exe.read_bytes()).hexdigest()}
    if args.provenance:
        identity['build']=json.loads(args.provenance.read_text(encoding='utf-8'))
        shutil.copy2(args.provenance,pkg/'BUILD_ENVIRONMENT.json')
    (pkg/'PACKAGE_ID.json').write_text(json.dumps(identity,indent=2))
    sums={str(p.relative_to(pkg)).replace(os.sep,'/'):hashlib.sha256(p.read_bytes()).hexdigest()
          for p in sorted(pkg.rglob('*')) if p.is_file()}
    (pkg/'SHA256SUMS').write_text(''.join(h+'  '+p+'\n' for p,h in sums.items()))
    if system=='Windows':archive=Path(shutil.make_archive(str(out/pkg.name),'zip',out,pkg.name))
    else:
        archive=out/(pkg.name+'.tar.gz')
        with tarfile.open(archive,'w:gz') as t:t.add(pkg,arcname=pkg.name)
    (out/(archive.name+'.sha256')).write_text(hashlib.sha256(archive.read_bytes()).hexdigest()+'  '+archive.name+'\n')
    print(archive,flush=True)

if __name__=='__main__':main()
