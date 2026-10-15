# CosmicScanner v2.0.0 平台安装与入门
本说明适用于 v2.0.0 平台程序包。各包的构建与验收信息见 PACKAGE_ID.json 及随附记录。

## Windows x86-64
将整个ZIP解压到自己选择的目录。不要只复制EXE：bin内的DLL必须保留。
在解压目录打开PowerShell：
~~~powershell
.\bin\cosmicscanner.exe --version
.\bin\cosmicscanner.exe --help
New-Item -ItemType Directory -Force demo_results | Out-Null
.\bin\cosmicscanner.exe gwas --model lmm --bfile examples/dataset/synth --pheno examples/dataset/quantitative.pheno --mpheno 1 --threads 1 --out demo_results/lmm
(Import-Csv -LiteralPath demo_results/lmm.assoc.txt -Delimiter ([char]9)).Count
~~~
最后应显示180，表示180个SNP结果。软件本身无需Python、Conda或编译器。
结果包括lmm.assoc.txt、lmm.cosmic.sumstats和配置/运行记录。
若找不到DLL，确认整个ZIP已解压且bin文件完整；不要从不明网站单独下载DLL。
当前包为64位Windows程序，不是Windows ARM64原生包。
中文路径的UTF-8支持要求Windows 10 1903或更新系统；实际验收在当前Windows主机完成。

## Linux x86-64
解压后从包目录运行：
~~~bash
./bin/cosmicscanner --version
mkdir -p demo_results
./bin/cosmicscanner gwas --model lmm --bfile examples/dataset/synth --pheno examples/dataset/quantitative.pheno --mpheno 1 --threads 1 --out demo_results/lmm
~~~
示例应生成180个SNP结果。当前可移植包要求glibc>=2.34，暂未承诺全部Linux发行版。

## macOS Apple 芯片 / Intel 芯片
系统要求：macOS 15 或更新。Apple 芯片选择 arm64 包，Intel 芯片选择 x86_64 包。
解压整个 tar.gz，保留 bin 和 lib。进入解压目录后运行：
~~~bash
./bin/cosmicscanner --version
./bin/cosmicscanner --help
mkdir -p demo_results
./bin/cosmicscanner gwas --model lmm --bfile examples/dataset/synth --pheno examples/dataset/quantitative.pheno --mpheno 1 --threads 1 --out demo_results/lmm
wc -l demo_results/lmm.assoc.txt
~~~
应得到 181 行（表头加 180 个 SNP）。运行软件无需 Homebrew、Conda、Python 或编译器。
可选自检需要 Python 3：python3 test/verify_native_package.py。
若提示 CPU 架构不匹配，重新选择对应芯片的软件包。
当前采用临时签名，尚未完成 Apple Developer ID 签名与公证；系统可能显示安全提示。
包内 PACKAGE_ID.json 和随附验收记录用于核对构建版本及测试范围。

## 输出、模型和正式范围
--mpheno从FID/IID之后的第一列表型开始计数。
使用具名输出列；v2默认单性状效应/P值使用BETA_A/SE_A/P_A。
多性状/重复力的输出列与单性状不同，请查CLI_GUIDE_CN.md。
帮助中可见的研究模型不等于正式支持承诺。
软件启动成功、示例运行通过和正式统计验收是不同门槛。
