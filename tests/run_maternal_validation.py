#!/usr/bin/env python3
"""
母体效应模型验证脚本
====================
验证目标：CosmicScanner --maternal-model 与 R sommer 包在方差分量与 SNP 效应上的一致性。

设计要点：
1. 模拟数据：n_animals 个幼畜 + n_dams 头母畜（母畜也是基因型个体）
   - 5 种方差分量配置，覆盖典型/高母体/低直接/负相关/无母体效应场景
   - 3 个因果 SNP：纯直接效应 / 纯母体效应 / 直接+母体
2. 输出 PLINK bed/bim/fam + 表型文件（FID IID DAM_ID Pheno）
3. 调用 CosmicScanner --maternal-model 跑 GWAS
4. 调用 R sommer 跑参考实现（mvmme(u = ~animal + dam, ...))
5. 对比 σ²_ad, σ²_am, σ²_e, β, SE, -log10P

依赖：
- Python 3.8+, numpy, pandas
- R + sommer 包（CRAN: install.packages("sommer"))
- CosmicScanner 二进制（在 PATH 中或用 --bin 指定）

用法：
  python3 run_maternal_validation.py --bin /path/to/cosmicscanner --out maternal_val
  python3 run_maternal_validation.py --bin cosmicscanner --out maternal_val --n-replicates 10

输出：
  maternal_val/
    sim_1/                # 每个模拟 replicate
      plink.bed/bim/fam
      pheno.txt           # FID IID DAM_ID Pheno
      cosmic.assoc.txt    # CosmicScanner GWAS 结果
      cosmic_null.txt     # 方差分量
      sommer_result.RData # R 参考结果
      sommer_assoc.txt    # sommer SNP 效应
    validation_report.md  # 汇总对比报告
    vc_comparison.csv     # 方差分量对比表
    beta_comparison.csv   # β 对比表
"""

import argparse
import os
import subprocess
import sys
import shutil
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd


# ============================================================
# 1. 模拟数据生成
# ============================================================

def simulate_genotypes(n_animals, n_snps, maf_range=(0.05, 0.5), seed=0):
    """
    模拟 PLINK 格式基因型。
    返回 genotype matrix (n_animals x n_snps), maf array, snp_ids, positions
    """
    rng = np.random.default_rng(seed)
    mafs = rng.uniform(maf_range[0], maf_range[1], n_snps)
    # 0/1/2 编码
    G = np.zeros((n_animals, n_snps), dtype=np.int8)
    for j in range(n_snps):
        p = mafs[j]
        # Hardy-Weinberg
        dosages = rng.choice([0, 1, 2], size=n_animals,
                             p=[(1-p)**2, 2*p*(1-p), p**2])
        G[:, j] = dosages
    snp_ids = [f"rs{j+1}" for j in range(n_snps)]
    positions = np.arange(1, n_snps + 1) * 1000
    return G, mafs, snp_ids, positions


def build_grm(genotypes, mafs):
    """
    VanRaden GRM: G = (M-P)(M-P)' / (2*sum(p(1-p)))
    """
    M = genotypes.astype(np.float64)
    P = 2.0 * mafs[np.newaxis, :]  # 1 x n_snps
    Z = M - P
    S = 2.0 * np.sum(mafs * (1 - mafs))
    G = Z @ Z.T / S
    return G


def simulate_maternal_pheno(genotypes, grm, record_to_animal, record_to_dam,
                            vc_config, causal_idx, causal_effects, seed=0):
    """
    模拟母体效应表型
    y = Xβ + Z_a*a + Z_m*m + e
    其中 a ~ N(0, σ²_ad * G), m ~ N(0, σ²_am * G_m)
    G_m = G (母畜也是基因型个体)

    vc_config: dict with keys sigma2_ad, sigma2_am, sigma_adm, sigma2_e
    causal_idx: list of SNP indices
    causal_effects: list of (beta_direct, beta_maternal) tuples

    返回: y (n_records,), a (n_subjects,), m (n_subjects,)
    """
    rng = np.random.default_rng(seed + 1000)
    n_records = len(record_to_animal)
    n_subjects = grm.shape[0]

    # 直接加性效应 a ~ N(0, σ²_ad * G)
    # 通过 Cholesky 分解生成相关随机效应
    try:
        L_g = np.linalg.cholesky(grm + 1e-6 * np.eye(n_subjects))
    except np.linalg.LinAlgError:
        # 加入 jitter
        L_g = np.linalg.cholesky(grm + 1e-4 * np.eye(n_subjects))

    a = L_g @ rng.standard_normal(n_subjects) * np.sqrt(vc_config["sigma2_ad"])

    # 母体加性效应 m ~ N(0, σ²_am * G_m)
    # 注：此处简化假设 G_m = G，与 CosmicScanner 当前实现一致
    if vc_config["sigma2_am"] > 0:
        m = L_g @ rng.standard_normal(n_subjects) * np.sqrt(vc_config["sigma2_am"])
        # 可选：direct-maternal 协方差（通过相关噪声近似）
        if abs(vc_config.get("sigma_adm", 0.0)) > 0:
            # 生成相关结构
            rho = vc_config["sigma_adm"] / math.sqrt(
                vc_config["sigma2_ad"] * vc_config["sigma2_am"] + 1e-10)
            rho = np.clip(rho, -0.99, 0.99)
            z_shared = rng.standard_normal(n_subjects)
            m = rho * (a / np.sqrt(vc_config["sigma2_ad"] + 1e-10)) * \
                np.sqrt(vc_config["sigma2_am"]) + \
                np.sqrt(1 - rho**2) * m
    else:
        m = np.zeros(n_subjects)

    # 残差
    e = rng.standard_normal(n_records) * np.sqrt(vc_config["sigma2_e"])

    # SNP 固定效应
    snp_eff = np.zeros(n_records)
    for idx, (bd, bm) in zip(causal_idx, causal_effects):
        # 直接效应：幼畜基因型
        g_animal = genotypes[record_to_animal, idx].astype(np.float64)
        # 母体效应：母畜基因型
        g_dam = genotypes[record_to_dam, idx].astype(np.float64)
        snp_eff += bd * g_animal + bm * g_dam

    # 组装表型（标准化到均值 0、方差 ~1 便于 VCE 收敛）
    y = snp_eff + a[record_to_animal] + m[record_to_dam] + e
    y = (y - y.mean())
    return y, a, m


def write_plink_files(outdir, genotypes, mafs, snp_ids, positions, n_animals):
    """写入 PLINK bed/bim/fam
    注：所有个体使用相同 FID（"0"），确保 maternal model 的 dam_key="FID:dam_id"
    能正确匹配到 iid_to_subject_idx_ 中的 "0:dam_iid"
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = outdir / "plink"

    # .fam — 统一 FID="0"
    with open(str(prefix) + ".fam", "w") as f:
        for i in range(n_animals):
            f.write(f"0\tIND{i+1}\t0\t0\t1\t-9\n")

    # .bim
    with open(str(prefix) + ".bim", "w") as f:
        for j, (sid, pos, maf) in enumerate(zip(snp_ids, positions, mafs)):
            a1 = "A" if maf < 0.5 else "T"  # minor allele
            a2 = "T" if a1 == "A" else "A"
            f.write(f"1\t{sid}\t0\t{pos}\t{a1}\t{a2}\n")

    # .bed (PLINK 1.9 binary format)
    n_snps = len(snp_ids)
    # 每个 SNP 占 ceil(n_animals/4) 字节
    n_bytes_per_snp = (n_animals + 3) // 4
    bed = np.zeros((n_snps, n_bytes_per_snp), dtype=np.uint8)

    for j in range(n_snps):
        for i in range(n_animals):
            dosage = genotypes[i, j]
            # PLINK .bed 编码（SNP-major, 2-bit per individual）：
            #   00 = homozygous A1 (dosage 0)
            #   01 = missing
            #   10 = heterozygous   (dosage 1)
            #   11 = homozygous A2 (dosage 2)
            code = {0: 0b00, 1: 0b10, 2: 0b11}.get(int(dosage), 0b01)
            byte_idx = i // 4
            bit_offset = (i % 4) * 2
            bed[j, byte_idx] |= (code << bit_offset)

    with open(str(prefix) + ".bed", "wb") as f:
        f.write(b'\x6c\x1b\x01')  # PLINK magic number
        f.write(bed.tobytes())

    return str(prefix)


def write_pheno_file(outdir, record_to_animal, record_to_dam, y, n_animals):
    """
    写入表型文件，格式：FID IID DAM_ID Pheno
    所有个体 FID 统一为 "0"，DAM_ID 为母畜 IID（与 .fam 一致）
    注：不带 header 行（maternal model 的 loadPhenotype 不跳过 header）
    """
    outdir = Path(outdir)
    pheno_path = outdir / "pheno.txt"
    with open(pheno_path, "w") as f:
        for r in range(len(record_to_animal)):
            ai = record_to_animal[r]
            di = record_to_dam[r]
            iid = f"IND{ai+1}"
            dam_id = f"IND{di+1}"
            f.write(f"0\t{iid}\t{dam_id}\t{y[r]:.6f}\n")
    return str(pheno_path)


def write_covar_file(outdir, n_animals):
    """写入定量协变量文件（Sex 0/1），FID 统一为 "0"，无 header """
    outdir = Path(outdir)
    rng = np.random.default_rng(42)
    sex = rng.integers(0, 2, n_animals)
    covar_path = outdir / "covar.txt"
    with open(covar_path, "w") as f:
        for i in range(n_animals):
            f.write(f"0\tIND{i+1}\t{sex[i]}\n")
    return str(covar_path)


# ============================================================
# 2. CosmicScanner 调用
# ============================================================

def run_cosmicscanner(binary, plink_prefix, pheno_file, covar_file,
                      out_prefix, thread_num=4):
    """运行 CosmicScanner 母体模型
    注：--maternal-model 只设置模型标志，需配合 --gwas 触发 GWAS 调度
    Sex 作为定量协变量传入（0/1 编码）
    """
    cmd = [
        binary,
        "--bfile", plink_prefix,
        "--pheno", pheno_file,
        "--mpheno", "1",
        "--gwas",              # 触发 GWAS 调度
        "--maternal-model",
        "--maternal-id", "3",
        "--qcovar", covar_file,  # Sex 0/1 作为定量协变量
        "--thread-num", str(thread_num),
        "--out", out_prefix,
    ]
    print("  [CosmicScanner] " + " ".join(cmd))
    # 清理 LD_LIBRARY_PATH，避免 Trae 沙箱的旧 libstdc++ 被抢先加载
    clean_env = dict(os.environ)
    traec_paths = [p for p in clean_env.get("LD_LIBRARY_PATH", "").split(":")
                   if "trae-cn-server" not in p]
    clean_env["LD_LIBRARY_PATH"] = ":".join(traec_paths)
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800,
                            env=clean_env)
    if result.returncode != 0:
        print("  STDERR:", result.stderr[-2000:])
        raise RuntimeError(f"CosmicScanner failed: {result.returncode}")
    return out_prefix + ".assoc.txt", out_prefix + "_maternal_null.txt"


def parse_cosmic_null(null_file):
    """解析 CosmicScanner 方差分量输出"""
    df = pd.read_csv(null_file, sep="\t")
    return {
        "sigma2_ad": float(df["sigma2_ad"].iloc[0]),
        "sigma2_am": float(df["sigma2_am"].iloc[0]),
        "sigma_adm": float(df["sigma_adm"].iloc[0]),
        "sigma2_c": float(df["sigma2_c"].iloc[0]),
        "sigma2_e": float(df["sigma2_e"].iloc[0]),
    }


def parse_cosmic_assoc(assoc_file):
    """解析 CosmicScanner 关联结果"""
    df = pd.read_csv(assoc_file, sep="\t")
    return df


# ============================================================
# 3. R sommer 参考实现
# ============================================================

SOMMER_SCRIPT = r'''
# R 参考实现：sommer 母体效应模型
library(sommer)

args <- commandArgs(trailingOnly = TRUE)
plink_prefix <- args[1]
pheno_file   <- args[2]
covar_file   <- args[3]
out_prefix   <- args[4]

# 读取表型
pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
colnames(pheno) <- c("FID", "IID", "DAM_ID", "Pheno")
covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
colnames(covar) <- c("FID", "IID", "Sex")
dat <- merge(pheno, covar, by = c("FID", "IID"))

# 读取基因型（简化：用 plink --recode A 导出的 0/1/2 矩阵）
# 这里假设外部已用 plink --recode A 生成 .raw 文件
raw_file <- paste0(plink_prefix, ".raw")
raw <- read.table(raw_file, header = TRUE)
# raw 第一列是 FID_IID，后面是 SNP 列 (格式 rs1_A, rs2_T, ...)
snp_cols <- grep("_[A-Z]$", colnames(raw), value = TRUE)
iid_raw <- sapply(strsplit(as.character(raw$FID_IID), "_"), `[`, 2)
raw$IID <- iid_raw
dat <- merge(dat, raw[, c("IID", snp_cols)], by = "IID")

# 构造 GRM（VanRaden）
G <- as.matrix(dat[, snp_cols])
n_snps <- ncol(G)
p <- colMeans(G) / 2
Z <- sweep(G, 2, 2 * p, "-")
S <- 2 * sum(p * (1 - p))
K <- tcrossprod(Z) / S

# 动物 ID 与母畜 ID 映射到 1..n
animals <- dat$IID
dams <- dat$DAM_ID
animal_idx <- match(animals, animals)
dam_idx <- match(dams, animals)
# 母畜若不在 animals 列表（理论上都在），用 NA 处理
dam_idx[is.na(dam_idx)] <- animal_idx[is.na(dam_idx)]

dat$animal <- animal_idx
dat$dam <- dam_idx

# sommer mvmme：direct + maternal 双随机效应
# 使用 GRM 作为方差结构
ans <- mmer(
  Pheno ~ Sex,
  random = ~ vm(animal, K) + vm(dam, K),
  rcov = ~ units,
  data = dat,
  verbose = FALSE
)

# 提取方差分量
vc_ad <- ans$sigma$`vm(animal, K)`[1,1]
vc_am <- ans$sigma$`vm(dam, K)`[1,1]
vc_e  <- ans$sigma$$units[1,1]

vc_df <- data.frame(
  sigma2_ad = vc_ad, sigma2_am = vc_am, sigma_adm = 0,
  sigma2_c = 0, sigma2_e = vc_e
)
write.table(vc_df, paste0(out_prefix, "_sommer_null.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# SNP 扫描：对每个 SNP 跑 mmer（小样本可行，大样本需优化）
results <- data.frame(
  SNP = character(), BETA = numeric(), SE = numeric(), P = numeric()
)
for (j in seq_along(snp_cols)) {
  snp_name <- sub("_[A-Z]$", "", snp_cols[j])
  dat$g <- dat[, snp_cols[j]]
  tryCatch({
    ans2 <- mmer(
      Pheno ~ Sex + g,
      random = ~ vm(animal, K) + vm(dam, K),
      rcov = ~ units,
      data = dat,
      verbose = FALSE
    )
    beta <- ans2$Beta[ans2$Beta$Effect == "g", "Estimate"]
    se <- ans2$Beta[ans2$Beta$Effect == "g", "Std.Error"]
    p_val <- 2 * pnorm(-abs(beta / se))
    results <- rbind(results, data.frame(
      SNP = snp_name, BETA = beta, SE = se, P = p_val
    ))
  }, error = function(e) {
    results <- rbind(results, data.frame(
      SNP = snp_name, BETA = NA, SE = NA, P = NA
    ))
  })
}

write.table(results, paste0(out_prefix, "_sommer_assoc.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("sommer reference done\n")
'''


def prepare_raw_for_r(plink_prefix):
    """用 plink --recode A 生成 .raw 文件供 R 读取"""
    raw_file = plink_prefix + ".raw"
    if not os.path.exists(raw_file):
        cmd = ["plink", "--bfile", plink_prefix, "--recode", "A",
               "--out", plink_prefix]
        subprocess.run(cmd, capture_output=True, check=False)
    return raw_file


def run_sommer_reference(plink_prefix, pheno_file, covar_file, out_prefix):
    """运行 R sommer 参考实现"""
    r_script_path = Path(out_prefix).parent / "sommer_reference.R"
    with open(r_script_path, "w") as f:
        f.write(SOMMER_SCRIPT)

    prepare_raw_for_r(plink_prefix)

    cmd = ["Rscript", str(r_script_path), plink_prefix, pheno_file,
           covar_file, out_prefix]
    print("  [sommer] " + " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    if result.returncode != 0:
        print("  R STDERR:", result.stderr[-3000:])
        # sommer 失败不致命，记录为 NA
        return None, None
    return (out_prefix + "_sommer_null.txt",
            out_prefix + "_sommer_assoc.txt")


def parse_sommer_null(null_file):
    if not os.path.exists(null_file):
        return None
    df = pd.read_csv(null_file, sep="\t")
    return {
        "sigma2_ad": float(df["sigma2_ad"].iloc[0]),
        "sigma2_am": float(df["sigma2_am"].iloc[0]),
        "sigma_adm": float(df["sigma_adm"].iloc[0]),
        "sigma2_c": float(df["sigma2_c"].iloc[0]),
        "sigma2_e": float(df["sigma2_e"].iloc[0]),
    }


def parse_sommer_assoc(assoc_file):
    if not os.path.exists(assoc_file):
        return None
    return pd.read_csv(assoc_file, sep="\t")


# ============================================================
# 4. 对比与报告
# ============================================================

def corr(x, y):
    """安全相关系数计算"""
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return float("nan")
    return float(np.corrcoef(x[mask], y[mask])[0, 1])


def generate_report(replicate_results, outdir):
    """生成对比报告"""
    outdir = Path(outdir)
    report_path = outdir / "validation_report.md"
    vc_csv = outdir / "vc_comparison.csv"
    beta_csv = outdir / "beta_comparison.csv"

    vc_rows = []
    beta_rows = []

    for rep in replicate_results:
        rep_id = rep["rep_id"]
        config_name = rep["config_name"]
        vc_truth = rep["vc_truth"]
        vc_cosmic = rep["vc_cosmic"]
        vc_sommer = rep.get("vc_sommer")

        for comp in ["sigma2_ad", "sigma2_am", "sigma2_e"]:
            row = {
                "rep": rep_id, "config": config_name,
                "component": comp,
                "truth": vc_truth[comp],
                "cosmic": vc_cosmic[comp] if vc_cosmic else None,
                "sommer": vc_sommer[comp] if vc_sommer else None,
                "cosmic_rel_err": (abs(vc_cosmic[comp] - vc_truth[comp])
                                   / vc_truth[comp]) if (vc_cosmic and vc_truth[comp] > 1e-10) else None,
                "sommer_rel_err": (abs(vc_sommer[comp] - vc_truth[comp])
                                  / vc_truth[comp]) if (vc_sommer and vc_truth[comp] > 1e-10) else None,
            }
            vc_rows.append(row)

        # β 对比
        cosmic_assoc = rep["cosmic_assoc"]
        sommer_assoc = rep.get("sommer_assoc")
        if sommer_assoc is not None:
            merged = cosmic_assoc.merge(
                sommer_assoc, left_on="SNP", right_on="SNP", suffixes=("_cosmic", "_sommer"))
            for _, r in merged.iterrows():
                beta_rows.append({
                    "rep": rep_id, "config": config_name,
                    "SNP": r["SNP"],
                    "beta_cosmic": r.get("BETA_cosmic"),
                    "beta_sommer": r.get("BETA_sommer"),
                    "se_cosmic": r.get("SE_cosmic"),
                    "se_sommer": r.get("SE_sommer"),
                })

    vc_df = pd.DataFrame(vc_rows)
    beta_df = pd.DataFrame(beta_rows)
    vc_df.to_csv(vc_csv, index=False)
    beta_df.to_csv(beta_csv, index=False)

    # 汇总统计
    vc_summary = vc_df.groupby("component").agg(
        cosmic_mean_rel_err=("cosmic_rel_err", "mean"),
        sommer_mean_rel_err=("sommer_rel_err", "mean"),
        n_reps=("cosmic_rel_err", "count"),
    ).reset_index()

    # β 相关
    beta_corr_beta = corr(beta_df["beta_cosmic"].values,
                          beta_df["beta_sommer"].values) if len(beta_df) > 0 else float("nan")
    beta_corr_se = corr(beta_df["se_cosmic"].values,
                        beta_df["se_sommer"].values) if len(beta_df) > 0 else float("nan")

    # 写报告
    with open(report_path, "w") as f:
        f.write("# 母体效应模型验证报告\n\n")
        f.write(f"**Replicates**: {len(replicate_results)}\n")
        f.write(f"**配置数**: {len(set(r['config_name'] for r in replicate_results))}\n\n")

        f.write("## 1. 方差分量对比\n\n")
        f.write("| Component | CosmicScanner 平均相对误差 | sommer 平均相对误差 | N |\n")
        f.write("|-----------|------------------------:|------------------:|---:|\n")
        for _, r in vc_summary.iterrows():
            ce = f"{r['cosmic_mean_rel_err']:.4f}" if pd.notna(r['cosmic_mean_rel_err']) else "NA"
            se = f"{r['sommer_mean_rel_err']:.4f}" if pd.notna(r['sommer_mean_rel_err']) else "NA"
            f.write(f"| {r['component']} | {ce} | {se} | {int(r['n_reps'])} |\n")

        f.write("\n## 2. SNP 效应对比\n\n")
        f.write(f"| Metric | Value |\n|--------|------:|\n")
        f.write(f"| β correlation (Cosmic vs sommer) | {beta_corr_beta:.6f} |\n")
        f.write(f"| SE correlation (Cosmic vs sommer) | {beta_corr_se:.6f} |\n")
        f.write(f"| Total SNPs compared | {len(beta_df)} |\n")

        f.write("\n## 3. 通过/失败判定\n\n")
        # 判定标准
        ad_err = vc_summary.loc[vc_summary.component == "sigma2_ad", "cosmic_mean_rel_err"]
        am_err = vc_summary.loc[vc_summary.component == "sigma2_am", "cosmic_mean_rel_err"]
        ad_pass = (ad_err.values[0] < 0.05) if len(ad_err) > 0 and pd.notna(ad_err.values[0]) else False
        am_pass = (am_err.values[0] < 0.10) if len(am_err) > 0 and pd.notna(am_err.values[0]) else False
        beta_pass = beta_corr_beta > 0.99 if not math.isnan(beta_corr_beta) else False

        f.write(f"- σ²_ad 相对误差 < 5%: {'**PASS**' if ad_pass else '**FAIL**'}\n")
        f.write(f"- σ²_am 相对误差 < 10%: {'**PASS**' if am_pass else '**FAIL**'}\n")
        f.write(f"- β(Cosmic vs sommer) > 0.99: {'**PASS**' if beta_pass else '**FAIL**'}\n")

        f.write("\n## 4. 逐 Replicate 详情\n\n")
        f.write("| Rep | Config | σ²_ad truth | σ²_ad cosmic | σ²_ad sommer | σ²_am cosmic | σ²_am sommer |\n")
        f.write("|----:|--------|------------:|-------------:|-------------:|-------------:|-------------:|\n")
        for rep in replicate_results:
            vt = rep["vc_truth"]
            vc = rep["vc_cosmic"] or {}
            vs = rep.get("vc_sommer") or {}
            ad_c = vc.get("sigma2_ad", float("nan"))
            am_c = vc.get("sigma2_am", float("nan"))
            ad_s = vs.get("sigma2_ad", float("nan"))
            am_s = vs.get("sigma2_am", float("nan"))
            f.write(f"| {rep['rep_id']} | {rep['config_name']} | "
                    f"{vt['sigma2_ad']:.3f} | {ad_c:.3f} | "
                    f"{ad_s:.3f} | "
                    f"{am_c:.3f} | "
                    f"{am_s:.3f} |\n")

    print(f"\n报告已生成: {report_path}")
    print(f"方差分量对比: {vc_csv}")
    print(f"β 对比: {beta_csv}")
    return report_path


# ============================================================
# 5. 主流程
# ============================================================

# 5 种方差分量配置
VC_CONFIGS = [
    {"name": "typical",       "sigma2_ad": 0.30, "sigma2_am": 0.10, "sigma_adm": -0.05, "sigma2_e": 0.55},
    {"name": "high_maternal", "sigma2_ad": 0.20, "sigma2_am": 0.25, "sigma_adm": -0.02, "sigma2_e": 0.50},
    {"name": "low_direct",    "sigma2_ad": 0.10, "sigma2_am": 0.10, "sigma_adm":  0.00, "sigma2_e": 0.70},
    {"name": "neg_corr",      "sigma2_ad": 0.30, "sigma2_am": 0.15, "sigma_adm": -0.10, "sigma2_e": 0.50},
    {"name": "no_maternal",   "sigma2_ad": 0.40, "sigma2_am": 0.00, "sigma_adm":  0.00, "sigma2_e": 0.60},
]


def main():
    parser = argparse.ArgumentParser(description="母体效应模型验证")
    parser.add_argument("--bin", required=True, help="cosmicscanner 可执行文件路径")
    parser.add_argument("--out", default="maternal_val", help="输出目录")
    parser.add_argument("--n-animals", type=int, default=1000, help="幼畜数")
    parser.add_argument("--n-dams", type=int, default=500, help="母畜数")
    parser.add_argument("--n-snps", type=int, default=1000, help="SNP 数")
    parser.add_argument("--n-replicates", type=int, default=2,
                        help="每个配置的重复次数（默认 2，正式跑设 10）")
    parser.add_argument("--configs", default="all",
                        help="逗号分隔的配置名，或 all")
    parser.add_argument("--skip-sommer", action="store_true",
                        help="跳过 sommer 参考实现（仅与真值对比）")
    parser.add_argument("--thread-num", type=int, default=4)
    args = parser.parse_args()

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    # 选择配置
    if args.configs == "all":
        configs = VC_CONFIGS
    else:
        names = set(args.configs.split(","))
        configs = [c for c in VC_CONFIGS if c["name"] in names]

    # 因果 SNP：3 个，分别模拟 纯直接 / 纯母体 / 直接+母体
    # 索引自适应：均匀分布在 n_snps 范围内，避免越界
    causal_idx = [max(1, args.n_snps // 10),
                  args.n_snps // 2,
                  args.n_snps - max(1, args.n_snps // 10)]
    causal_effects = [
        (0.5, 0.0),   # 纯直接
        (0.0, 0.5),   # 纯母体
        (0.3, 0.3),   # 直接+母体
    ]

    # 总样本数 = 幼畜 + 母畜（假设母畜也是基因型个体）
    # 简化设计：n_animals 个个体中前 n_dams 个为母畜
    n_total = max(args.n_animals, args.n_dams)
    n_animals = args.n_animals
    n_dams = min(args.n_dams, n_animals)  # 母畜从动物池中抽取

    all_results = []
    rep_counter = 0

    for ci, cfg in enumerate(configs):
        vc_config = {k: v for k, v in cfg.items() if k != "name"}
        print(f"\n=== Config {ci+1}/{len(configs)}: {cfg['name']} ===")
        print(f"  VCs: {vc_config}")

        for rep in range(args.n_replicates):
            rep_counter += 1
            print(f"\n--- Replicate {rep_counter} (config={cfg['name']}, rep={rep+1}) ---")
            rep_outdir = outdir / f"sim_{rep_counter}"
            rep_outdir.mkdir(parents=True, exist_ok=True)

            seed = rep_counter * 1000 + ci
            rng = np.random.default_rng(seed)

            # 1. 模拟基因型
            print("  [1/5] Simulating genotypes...")
            genotypes, mafs, snp_ids, positions = simulate_genotypes(
                n_animals, args.n_snps, seed=seed)

            # 2. 构造 GRM
            print("  [2/5] Building GRM...")
            grm = build_grm(genotypes, mafs)

            # 3. 构造 record → animal / dam 映射
            # 每个幼畜 1 条记录，母畜从基因型池中随机抽
            n_records = n_animals
            record_to_animal = np.arange(n_animals)
            # 母畜从前 n_dams 个个体中随机抽取（模拟母畜也是基因型个体）
            record_to_dam = rng.integers(0, n_dams, n_records)

            # 4. 模拟表型
            print("  [3/5] Simulating phenotypes...")
            y, a, m = simulate_maternal_pheno(
                genotypes, grm, record_to_animal, record_to_dam,
                vc_config, causal_idx, causal_effects, seed=seed)

            # 5. 写文件
            print("  [4/5] Writing input files...")
            plink_prefix = write_plink_files(
                rep_outdir, genotypes, mafs, snp_ids, positions, n_animals)
            pheno_file = write_pheno_file(
                rep_outdir, record_to_animal, record_to_dam, y, n_animals)
            covar_file = write_covar_file(rep_outdir, n_animals)

            # 保存真值
            truth = {
                "vc_truth": vc_config,
                "causal_idx": causal_idx,
                "causal_effects": causal_effects,
            }
            with open(rep_outdir / "truth.json", "w") as f:
                json.dump(truth, f, indent=2, default=str)

            # 6. 运行 CosmicScanner
            print("  [5/5a] Running CosmicScanner...")
            cosmic_out = str(rep_outdir / "cosmic")
            try:
                cosmic_assoc_file, cosmic_null_file = run_cosmicscanner(
                    args.bin, plink_prefix, pheno_file, covar_file,
                    cosmic_out, thread_num=args.thread_num)
                vc_cosmic = parse_cosmic_null(cosmic_null_file)
                cosmic_assoc = parse_cosmic_assoc(cosmic_assoc_file)
                print(f"    σ²_ad={vc_cosmic['sigma2_ad']:.4f} (truth={vc_config['sigma2_ad']:.4f})")
                print(f"    σ²_am={vc_cosmic['sigma2_am']:.4f} (truth={vc_config['sigma2_am']:.4f})")
                print(f"    σ²_e ={vc_cosmic['sigma2_e']:.4f} (truth={vc_config['sigma2_e']:.4f})")
            except Exception as e:
                print(f"    CosmicScanner FAILED: {e}")
                vc_cosmic = None
                cosmic_assoc = pd.DataFrame()

            # 7. 运行 sommer 参考（可选）
            vc_sommer = None
            sommer_assoc = None
            if not args.skip_sommer:
                print("  [5/5b] Running sommer reference...")
                try:
                    sommer_null, sommer_assoc_file = run_sommer_reference(
                        plink_prefix, pheno_file, covar_file, cosmic_out)
                    if sommer_null:
                        vc_sommer = parse_sommer_null(sommer_null)
                        print(f"    sommer σ²_ad={vc_sommer['sigma2_ad']:.4f}")
                        print(f"    sommer σ²_am={vc_sommer['sigma2_am']:.4f}")
                    if sommer_assoc_file:
                        sommer_assoc = parse_sommer_assoc(sommer_assoc_file)
                except Exception as e:
                    print(f"    sommer reference FAILED: {e}")

            all_results.append({
                "rep_id": rep_counter,
                "config_name": cfg["name"],
                "vc_truth": vc_config,
                "vc_cosmic": vc_cosmic,
                "vc_sommer": vc_sommer,
                "cosmic_assoc": cosmic_assoc,
                "sommer_assoc": sommer_assoc,
            })

    # 生成报告
    print("\n=== Generating report ===")
    generate_report(all_results, outdir)


if __name__ == "__main__":
    main()
