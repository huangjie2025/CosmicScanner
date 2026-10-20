# CosmicScanner v2.0.0 平台安装与入门

## 第2章 安装与运行

### 2.1 选择程序包与核验文件

按操作系统和处理器选择下表中的程序包，同时下载同名 .sha256 校验文件。四种包均包含程序、配套运行库、示例数据和使用说明。分析程序本身不需要 Python、Conda、Homebrew 或编译器；只有可选自检脚本需要 Python 3。

| 平台 | 程序包名称 |
|---|---|
| Linux x86-64 | CosmicScanner-v2.0.0-rc1-linux-x86_64.tar.gz |
| Windows x86-64 | CosmicScanner-v2.0.0-rc1-windows-x86_64.zip |
| macOS Apple 芯片 | CosmicScanner-v2.0.0-rc1-macos-arm64.tar.gz |
| macOS Intel 芯片 | CosmicScanner-v2.0.0-rc1-macos-x86_64.tar.gz |

Linux 要求 glibc 2.34 或更新；Windows 建议 Windows 10 1903 或更新的64位系统；macOS 要求15或更新。Windows ARM64、Linux ARM64 和32位系统没有对应原生程序包。Apple 芯片用户选择 arm64，不需要使用 Intel 版本。

Linux 在下载文件所在目录核验：

```bash
sha256sum -c CosmicScanner-v2.0.0-rc1-linux-x86_64.tar.gz.sha256
```

macOS 使用系统自带 shasum。以下为 Apple 芯片包；Intel 用户将文件名中的 arm64 替换为 x86_64：

```bash
shasum -a 256 -c CosmicScanner-v2.0.0-rc1-macos-arm64.tar.gz.sha256
```

显示 OK 表示下载内容与校验文件一致。Windows 在 PowerShell 中计算校验值，并与同名 .sha256 文件的第一列比较：

```powershell
Get-FileHash -Algorithm SHA256 -LiteralPath CosmicScanner-v2.0.0-rc1-windows-x86_64.zip
Get-Content -LiteralPath CosmicScanner-v2.0.0-rc1-windows-x86_64.zip.sha256
```

校验值不区分字母大小写。若不一致，重新获取程序包与配套校验文件。保存 PACKAGE_ID.json 和压缩包校验值，便于记录实际使用的构建。

### 2.2 Linux 安装与首次运行

在选定的软件目录解压，进入解压生成的程序包根目录：

```bash
tar -xzf CosmicScanner-v2.0.0-rc1-linux-x86_64.tar.gz
cd CosmicScanner-v2.0.0-rc1-linux-x86_64
./bin/cosmicscanner --version
./bin/cosmicscanner --help
```

保留 bin、lib、libexec 等包内目录及启动脚本，不要只移动主程序。启动器负责加载配套依赖。若出现 glibc 版本不足，换用符合要求的系统，或按该版本源码构建说明编译；不要替换系统核心库。

### 2.3 Windows 安装与首次运行

将 ZIP 保存并解压到自己选择的目录，例如 E 盘的软件目录。不要在压缩包预览界面直接运行。以下命令在 ZIP 所在目录的 PowerShell 中执行：

```powershell
Expand-Archive -LiteralPath CosmicScanner-v2.0.0-rc1-windows-x86_64.zip -DestinationPath .
Set-Location .\CosmicScanner-v2.0.0-rc1-windows-x86_64
& '.\bin\cosmicscanner.exe' --version
& '.\bin\cosmicscanner.exe' --help
```

bin 内的 DLL 必须与 exe 一起保留。PowerShell 调用路径变量或带引号的程序路径时使用 &；含空格的文件路径要加引号。软件和数据均可放在选定磁盘，不需要复制到 C 盘。若目标目录已存在，先保留旧版及其运行记录，再为新包选择独立目录。

### 2.4 macOS 安装与首次运行

在“关于本机”查看芯片类型：Apple M 系列使用 arm64，Intel 处理器使用 x86_64。以下为 Apple 芯片操作；Intel 用户将两处 arm64 替换为 x86_64。

```bash
tar -xzf CosmicScanner-v2.0.0-rc1-macos-arm64.tar.gz
cd CosmicScanner-v2.0.0-rc1-macos-arm64
./bin/cosmicscanner --version
./bin/cosmicscanner --help
```

完整保留 bin 和 lib，程序会加载包内的 libomp.dylib。出现“bad CPU type”时，检查处理器架构与下载包是否对应。出现“Permission denied”时，先确认已完整解压到可执行文件的本地目录，并检查 bin/cosmicscanner 是否保留执行权限。

当前 macOS 包采用临时签名，尚未做 Apple Developer ID 签名与公证。通过浏览器下载后，系统可能提示无法验证开发者。确认下载来源和 SHA-256 后，按系统设置中的“隐私与安全性”对该程序进行处理，操作说明见 Apple 官方支持：https://support.apple.com/en-ie/102445。云端运行通过不能代替这一下载后的系统检查。

macOS 可使用第4.2节的完整终端命令，参数与 Linux 相同。

### 2.5 运行随包自检

首次使用、升级版本或更换计算机后，先按第4章运行四条示例。若已安装 Python 3，也可在程序包根目录运行随包自检；该脚本只使用 Python 标准库。

Linux 和 macOS：

```bash
python3 test/verify_native_package.py
```

Windows：

```powershell
python test/verify_native_package.py
```

成功时最后输出包含 status 为 PASS、lmm、mv、rep 各180个位点，以及 grm_samples 为60的记录，并给出自检结果目录。输出目录位于系统临时目录，不会覆盖自己的研究结果。未安装 Python 时，直接使用第4章命令即可，无需为常规分析安装 Python。

若自检失败，先检查系统版本、处理器架构、文件完整性和程序包目录。保存完整命令与错误信息，并附上 PACKAGE_ID.json 中的版本信息联系开发者。自检验证的是随包合成示例，自己的数据仍需检查 ID 匹配、缺失、协变量、拟合及结果数量。

### 2.6 命令组织与帮助

推荐顺序为任务、模型、输入、模型专用参数、资源和输出。gwas 与 grm 用于明确任务，原有平铺命令仍可使用。

```text
cosmicscanner gwas --model lmm --bfile PREFIX --pheno FILE --out PREFIX
cosmicscanner grm --bfile PREFIX --out PREFIX
cosmicscanner help models
cosmicscanner help input
cosmicscanner --help-all
```

在实际终端中，将 cosmicscanner 替换为本平台的程序路径。--model 选择模型家族；--gwas-model 选择单性状内部效应结构；--scan-tier 选择逐位点方差处理策略。完整帮助见附录A，中文常用参数索引见附录B。

### 2.7 从源码构建

源码构建适合开发者和需要自定义构建的用户。获取与主项目匹配的 libCosmicCore、第三方依赖和完整构建说明。仅下载主目录或混用不同 Core 版本，可能导致编译失败或行为差异。

源码构建的常见程序位置为 ./build/cosmicscanner，源码示例在 examples/data/r33_synthetic/dataset。平台程序包使用 ./bin/cosmicscanner 和 examples/dataset。请按获取方式选用对应路径，构建后同样运行第4章示例。


完整用户手册见 USER_MANUAL_CN.md 或 CosmicScanner_v2.0.0_User_Manual_CN.docx。
