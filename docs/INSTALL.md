# Install CosmicScanner v2.0.0-rc1

Verify the archive SHA-256, extract the complete package, and open a terminal in
the extracted directory. Keep bin/, lib/, and libexec/ (where supplied) together.
The distribution is a release candidate. The unchanged program reports v2.0.0;
retain PACKAGE_ID.json and the archive checksum with your analysis.

## Requirements

- Windows x86_64: Windows 10 1903 or later recommended; keep DLLs beside the EXE.
- Linux x86_64: glibc 2.34 or later and /lib64/ld-linux-x86-64.so.2. Use the
  bin/cosmicscanner launcher, which selects the bundled runtime libraries.
- macOS arm64 / x86_64: macOS 15 or later. Use the package matching your CPU.
  Packages retain their original ad-hoc signature; they are not Developer ID
  signed or notarized.

Do not replace system libraries to work around an unsupported operating system.
No Python, Conda, Homebrew, or compiler is required for a normal analysis.

## First run

Linux / macOS: ./bin/cosmicscanner --version

Windows PowerShell: & '.\bin\cosmicscanner.exe' --version

Read ../README.txt for the complete bundled LMM command and expected output.
USER_MANUAL_EN.md and USER_MANUAL_CN.md are the complete manuals;
NATIVE_INSTALL_CN.md contains detailed Chinese platform guidance.
Optional verification (Python 3): python test/verify_native_package.py from the
package root. Use python3 on systems where that is the Python 3 command.
