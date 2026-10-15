# CosmicScanner v2.0.0 平台安装与入门
更新：2026-10-14。当前为验收材料，尚未签署稳定发布或公开下载。
Linux/Windows已有可运行包；macOS包仍需原生构建和验收。

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

## macOS：待完成验收
需分别选择Apple芯片arm64和Intel芯片x86_64包，不能互相当作原生包。
当前已准备原生构建/自检流水线，但还没有经过Mac实机验收的可下载包。
目标首轮验收环境为macOS 15；更早系统尚未验收。
未来成功验收的包使用./bin/cosmicscanner，命令参数与Linux相同。
临时签名不等于Apple Developer ID签名或公证，发布页必须如实标明。

## 输出、模型和正式范围
--mpheno从FID/IID之后的第一列表型开始计数。
使用具名输出列；v2默认单性状效应/P值使用BETA_A/SE_A/P_A。
多性状/重复力的输出列与单性状不同，请查CLI_GUIDE_CN.md。
帮助中可见的研究模型不等于正式支持承诺。
软件启动成功、示例运行通过和正式统计验收是不同门槛。
