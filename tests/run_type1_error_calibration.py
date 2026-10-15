#!/usr/bin/env python3
"""
CosmicScanner Type-I Error Calibration Framework
=================================================

目的：在零假设（无因果 SNP）下校准各 GWAS 模型的 p-value 分布。
      理想情况下 p ~ U(0,1)，λ_GC ≈ 1.0，type-I error rate ≈ α。

覆盖模型：
  1. LMM          —— 单性状线性混合模型
  2. mvLMM        —— 多性状多元线性混合模型
  3. Repeatability —— 重复力模型（多个体记录/个体）
  4. RRM          —— 随机回归模型（纵向轨迹）

核心指标（每个模型 × 每个 replicate）：
  - λ_GC         基因组膨胀因子 = median(χ²_obs) / 0.456
  - Type-I error 在 α = 0.05 / 0.01 / 0.001 / 5e-8 下的经验错误率
  - KS p-value   p-value 与 U(0,1) 的 Kolmogorov-Smirnov 检验

判定标准（across replicates mean）：
  - PASS: λ_GC ∈ [0.95, 1.05]  且  type-I error ≤ 2 × α
  - WARN: λ_GC ∈ [0.90, 1.10]  或  type-I error ∈ (2α, 3α]
  - FAIL: 其他

用法：
  python3 run_type1_error_calibration.py \
      --bin build/cosmicscanner --out type1_calib \
      --n-animals 500 --n-snps 1000 --n-replicates 10 \
      --models lmm,mvlmm,repeatability,rrm --thread-num 4

快速冒烟测试：
  python3 run_type1_error_calibration.py \
      --bin build/cosmicscanner --out /tmp/type1_smoke \
      --n-animals 100 --n-snps 200 --n-replicates 2 \
      --models lmm --thread-num 2
"""

import argparse
import csv
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

# ============================================================
# 1. 基因型模拟
# ============================================================

def simulate_genotypes(n_animals, n_snps, seed=0):
    """模拟基因型矩阵 (n_animals × n_snps)，0/1/2 编码"""
    rng = np.random.default_rng(seed)
    mafs = rng.uniform(0.05, 0.5, n_snps)
    # Binomial(2, maf) 基因型
    genotypes = np.zeros((n_animals, n_snps), dtype=np.int8)
    for j in range(n_snps):
        genotypes[:, j] = rng.binomial(2, mafs[j], n_animals)
    return genotypes, mafs


def build_grm(genotypes):
    """VanRaden GRM (第 I 类)"""
    G = genotypes.astype(np.float64)
    p = G.mean(axis=0) / 2.0
    Z = G - 2.0 * p[np.newaxis, :]
    denom = 2.0 * (p * (1.0 - p)).sum()
    if denom < 1e-12:
        return np.eye(G.shape[0])
    K = (Z @ Z.T) / denom
    # 对称化 + 对角线微调保证正定
    K = (K + K.T) / 2.0
    return K


def write_plink_files(outdir, genotypes, mafs, n_animals):
    """写入 PLINK bed/bim/fam（统一 FID="0"）"""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = outdir / "plink"
    n_snps = genotypes.shape[1]

    with open(str(prefix) + ".fam", "w") as f:
        for i in range(n_animals):
            f.write(f"0\tIND{i+1}\t0\t0\t1\t-9\n")

    with open(str(prefix) + ".bim", "w") as f:
        for j in range(n_snps):
            f.write(f"1\trs{j+1}\t0\t{(j+1)*1000}\tA\tT\n")

    def code(g):
        return 0b00 if g == 2 else (0b10 if g == 1 else 0b11)

    n_bytes = (n_animals + 3) // 4
    bed = bytearray(b"\x6c\x1b\x01")
    for j in range(n_snps):
        block = bytearray(n_bytes)
        for i in range(n_animals):
            block[i // 4] |= code(int(genotypes[i, j])) << ((i % 4) * 2)
        bed.extend(block)
    with open(str(prefix) + ".bed", "wb") as f:
        f.write(bytes(bed))
    return str(prefix)


# ============================================================
# 2. 零假设表型模拟（无 SNP 效应）
# ============================================================

def simulate_null_pheno_lmm(grm, sigma2_g=0.3, sigma2_e=0.7, seed=0):
    """LMM 零假设：y = a + e, a ~ N(0, σ²_g G), e ~ N(0, σ²_e I)"""
    rng = np.random.default_rng(seed)
    n = grm.shape[0]
    L = np.linalg.cholesky(grm + 1e-6 * np.eye(n))
    a = L @ rng.standard_normal(n) * math.sqrt(sigma2_g)
    e = rng.standard_normal(n) * math.sqrt(sigma2_e)
    y = a + e
    y = (y - y.mean())
    return y


def simulate_null_pheno_mvlmm(grm, sigma2_g=0.3, sigma2_e=0.7,
                               rg=0.5, re=0.3, seed=0):
    """mvLMM 零假设：2 个性状，遗传相关 rg、残差相关 re"""
    rng = np.random.default_rng(seed)
    n = grm.shape[0]
    L = np.linalg.cholesky(grm + 1e-6 * np.eye(n))

    # 遗传效应：a1, a2 相关 rg
    z1 = rng.standard_normal(n)
    z2 = rg * z1 + math.sqrt(1 - rg**2) * rng.standard_normal(n)
    a1 = L @ z1 * math.sqrt(sigma2_g)
    a2 = L @ z2 * math.sqrt(sigma2_g)

    # 残差：e1, e2 相关 re
    w1 = rng.standard_normal(n)
    w2 = re * w1 + math.sqrt(1 - re**2) * rng.standard_normal(n)
    e1 = w1 * math.sqrt(sigma2_e)
    e2 = w2 * math.sqrt(sigma2_e)

    y1 = a1 + e1
    y2 = a2 + e2
    y1 = (y1 - y1.mean())
    y2 = (y2 - y2.mean())
    return y1, y2


def simulate_null_pheno_repeatability(grm, n_records_per_animal=3,
                                      sigma2_g=0.3, sigma2_pe=0.2,
                                      sigma2_e=0.5, seed=0):
    """重复力模型零假设：y = a + pe + e，每个个体多条记录"""
    rng = np.random.default_rng(seed)
    n = grm.shape[0]
    L = np.linalg.cholesky(grm + 1e-6 * np.eye(n))
    a = L @ rng.standard_normal(n) * math.sqrt(sigma2_g)
    pe = rng.standard_normal(n) * math.sqrt(sigma2_pe)

    records = []
    for i in range(n):
        for r in range(n_records_per_animal):
            e = rng.standard_normal() * math.sqrt(sigma2_e)
            y = a[i] + pe[i] + e
            records.append((i, y))
    return records


def simulate_null_pheno_rrm(grm, time_points=(0.0, 0.25, 0.5, 0.75, 1.0),
                            sigma2_a0=0.3, sigma2_a1=0.05,
                            cov_a01=0.0, sigma2_e=0.4, seed=0):
    """RRM 零假设：y(t) = a_0 + a_1*t + e(t)
    a_0, a_1 为随机回归系数（Legendre 0/1 阶）
    """
    rng = np.random.default_rng(seed)
    n = grm.shape[0]
    L = np.linalg.cholesky(grm + 1e-6 * np.eye(n))

    # a0, a1 有协方差
    z0 = rng.standard_normal(n)
    z1 = rng.standard_normal(n)
    if sigma2_a1 > 0 and abs(cov_a01) > 0:
        rho = cov_a01 / math.sqrt(sigma2_a0 * sigma2_a1)
        rho = np.clip(rho, -0.99, 0.99)
        z1 = rho * z0 + math.sqrt(1 - rho**2) * z1
    a0 = L @ z0 * math.sqrt(sigma2_a0)
    a1 = L @ z1 * math.sqrt(sigma2_a1) if sigma2_a1 > 0 else np.zeros(n)

    records = []
    for i in range(n):
        for t in time_points:
            e = rng.standard_normal() * math.sqrt(sigma2_e)
            y = a0[i] + a1[i] * t + e
            records.append((i, t, y))
    return records


def simulate_null_pheno_maternal(grm, n_dams_ratio=0.5,
                                 sigma2_ad=0.3, sigma2_am=0.1,
                                 sigma2_e=0.6, seed=0):
    """母体模型零假设：y = a_d + a_m(dam) + e（无 SNP 效应）
    每个个体有一个母畜 ID，母畜从 n_dams 中随机分配
    """
    rng = np.random.default_rng(seed)
    n = grm.shape[0]
    L = np.linalg.cholesky(grm + 1e-6 * np.eye(n))
    # 直接加性效应
    a_d = L @ rng.standard_normal(n) * math.sqrt(sigma2_ad)
    # 母体加性效应（共享同一 GRM）
    a_m = L @ rng.standard_normal(n) * math.sqrt(sigma2_am)

    n_dams = max(1, int(n * n_dams_ratio))
    # 每个个体的母畜从种群中随机分配
    dam_indices = rng.integers(0, n_dams, n)

    records = []
    for i in range(n):
        e = rng.standard_normal() * math.sqrt(sigma2_e)
        y = a_d[i] + a_m[dam_indices[i]] + e
        records.append((i, dam_indices[i], y))
    return records


def simulate_null_pheno_threshold(grm, prevalence=0.3,
                                  sigma2_g=0.5, seed=0):
    """二元性状零假设（probit 连接）：
    隐患势 liability = a + e，y = 1 if liability > threshold else 0
    threshold 由 prevalence 决定（使 P(y=1) = prevalence）
    """
    rng = np.random.default_rng(seed)
    n = grm.shape[0]
    L = np.linalg.cholesky(grm + 1e-6 * np.eye(n))
    a = L @ rng.standard_normal(n) * math.sqrt(sigma2_g)
    e = rng.standard_normal(n) * math.sqrt(1.0 - sigma2_g)  # 总方差=1
    liability = a + e
    # 阈值：使 P(y=1) ≈ prevalence
    threshold = float(np.quantile(liability, 1.0 - prevalence))
    y_binary = (liability > threshold).astype(int)
    return y_binary


def simulate_null_pheno_gxe(grm, env_mean=20.0, env_sd=5.0,
                            sigma2_g=0.3, sigma2_e=0.7, seed=0):
    """GxE 零假设：y = a + e（无 SNP、无交互），环境变量独立
    表型格式：FID IID ENV VALUE
    """
    rng = np.random.default_rng(seed)
    n = grm.shape[0]
    L = np.linalg.cholesky(grm + 1e-6 * np.eye(n))
    a = L @ rng.standard_normal(n) * math.sqrt(sigma2_g)
    e = rng.standard_normal(n) * math.sqrt(sigma2_e)
    env = rng.normal(env_mean, env_sd, n)
    y = a + e
    y = (y - y.mean())
    return env, y


def simulate_null_pheno_survival(grm, event_rate=0.6, seed=0):
    """Survival 零假设：生存时间与基因型无关（指数分布 + 随机删失）
    表型格式：FID IID TIME EVENT
    event_rate: 目标事件率（约 60% 个体发生事件）
    返回 (time, event)
    """
    rng = np.random.default_rng(seed)
    n = grm.shape[0]
    # 生存时间 ~ Exp(1)，与基因型无关
    true_time = rng.exponential(scale=1.0, size=n)
    # 删失时间 ~ Uniform，使事件率 ≈ event_rate
    # P(T < C) = event_rate → C 的分位数 = T 的 event_rate 分位数
    c_threshold = float(np.quantile(true_time, event_rate))
    cens = rng.uniform(c_threshold, c_threshold * 2.0, n)
    observed_time = np.minimum(true_time, cens)
    event = (true_time <= cens).astype(int)
    return observed_time, event


def simulate_null_pheno_count(grm, family="poisson", mean_count=5.0,
                              binomial_trials=10, binomial_rate=0.3, seed=0):
    """Count/binomial PSR 零假设：表型与基因型无关
    Poisson: y ~ Poisson(λ)，λ = mean_count
    Binomial: y/N ~ Binomial(N, p)/N，N = binomial_trials, p = binomial_rate
    返回 (y_counts, trials) trials 为 None 表示 Poisson
    """
    rng = np.random.default_rng(seed)
    n = grm.shape[0]
    if family == "poisson":
        y = rng.poisson(lam=mean_count, size=n).astype(float)
        return y, None
    elif family == "binomial":
        # y_(idx) 存事件率 (event/trial)，trials 存 N
        events = rng.binomial(n=binomial_trials, p=binomial_rate, size=n)
        y_rate = events.astype(float) / binomial_trials
        trials = np.full(n, float(binomial_trials))
        return y_rate, trials
    else:
        raise ValueError(f"Unsupported count family: {family}")


def write_pheno_count(outdir, y, trials, n_animals, family="poisson"):
    """Count 模型表型：
    Poisson: FID IID COUNT（mpheno=1）
    Binomial: FID IID EVENT TRIAL（用 --binomial-pheno 单独文件）
    """
    outdir = Path(outdir)
    if family == "poisson":
        path = outdir / "pheno_count.txt"
        with open(path, "w") as f:
            for i in range(n_animals):
                f.write(f"0\tIND{i+1}\t{int(y[i])}\n")
        return str(path), None
    else:
        # Binomial: --pheno 提供事件率（mpheno=1），--binomial-pheno 提供 FID IID EVENT TRIAL
        pheno_path = outdir / "pheno_binom.txt"
        binom_path = outdir / "binom_pheno.txt"
        with open(pheno_path, "w") as f:
            for i in range(n_animals):
                f.write(f"0\tIND{i+1}\t{y[i]:.8f}\n")
        with open(binom_path, "w") as f:
            for i in range(n_animals):
                f.write(f"0\tIND{i+1}\t{int(y[i]*trials[i])}\t{int(trials[i])}\n")
        return str(pheno_path), str(binom_path)


# ============================================================
# 3. 表型文件写入
# ============================================================

def write_pheno_lmm(outdir, y, n_animals):
    path = Path(outdir) / "pheno_lmm.txt"
    with open(path, "w") as f:
        for i in range(n_animals):
            f.write(f"0\tIND{i+1}\t{y[i]:.8f}\n")
    return str(path)


def write_pheno_mvlmm(outdir, y1, y2, n_animals):
    path = Path(outdir) / "pheno_mv.txt"
    with open(path, "w") as f:
        for i in range(n_animals):
            f.write(f"0\tIND{i+1}\t{y1[i]:.8f}\t{y2[i]:.8f}\n")
    return str(path)


def write_pheno_repeatability(outdir, records, n_animals):
    """重复力模型表型：FID IID Pheno（多条记录/个体）"""
    path = Path(outdir) / "pheno_rep.txt"
    with open(path, "w") as f:
        for animal_idx, y in records:
            f.write(f"0\tIND{animal_idx+1}\t{y:.8f}\n")
    return str(path)


def write_pheno_rrm(outdir, records, n_animals):
    """RRM 表型：FID IID TIME Pheno"""
    path = Path(outdir) / "pheno_rrm.txt"
    with open(path, "w") as f:
        for animal_idx, t, y in records:
            f.write(f"0\tIND{animal_idx+1}\t{t:.4f}\t{y:.8f}\n")
    return str(path)


def write_pheno_maternal(outdir, records, n_animals):
    """母体模型表型：FID IID DAM_ID Pheno（无 header）
    注：maternal model 的 loadPhenotype 不跳过 header（已修复但保守起见仍无 header）
    """
    path = Path(outdir) / "pheno_maternal.txt"
    with open(path, "w") as f:
        for animal_idx, dam_idx, y in records:
            f.write(f"0\tIND{animal_idx+1}\tIND{dam_idx+1}\t{y:.8f}\n")
    return str(path)


def write_pheno_threshold(outdir, y_binary, n_animals):
    """二元性状表型：FID IID Pheno（0/1）"""
    path = Path(outdir) / "pheno_threshold.txt"
    with open(path, "w") as f:
        for i in range(n_animals):
            f.write(f"0\tIND{i+1}\t{int(y_binary[i])}\n")
    return str(path)


def write_pheno_gxe(outdir, env, y, n_animals):
    """GxE 表型：FID IID ENV VALUE（env-covar 指向第 3 列）"""
    path = Path(outdir) / "pheno_gxe.txt"
    with open(path, "w") as f:
        for i in range(n_animals):
            f.write(f"0\tIND{i+1}\t{env[i]:.4f}\t{y[i]:.8f}\n")
    return str(path)


def write_pheno_survival(outdir, time, event, n_animals):
    """Survival 表型：FID IID TIME EVENT（time_col=3, event_col=4）"""
    path = Path(outdir) / "pheno_survival.txt"
    with open(path, "w") as f:
        for i in range(n_animals):
            f.write(f"0\tIND{i+1}\t{time[i]:.8f}\t{int(event[i])}\n")
    return str(path)


def write_covar_file(outdir, n_animals):
    """写入定量协变量文件（Sex 0/1），用于 maternal 模型"""
    outdir = Path(outdir)
    rng = np.random.default_rng(42)
    sex = rng.integers(0, 2, n_animals)
    covar_path = outdir / "covar.txt"
    with open(covar_path, "w") as f:
        for i in range(n_animals):
            f.write(f"0\tIND{i+1}\t{int(sex[i])}\n")
    return str(covar_path)


# ============================================================
# 4. CosmicScanner 调用
# ============================================================

def run_cosmicscanner(binary, plink_prefix, model, pheno_file,
                      out_prefix, thread_num=4, extra_args=None, binom_file=None):
    """运行指定模型的 GWAS"""
    cmd = [binary, "--bfile", plink_prefix, "--maf", "0.0",
           "--no-loco", "--thread-num", str(thread_num), "--out", out_prefix]

    if model == "lmm":
        cmd += ["--pheno", pheno_file, "--mpheno", "1", "--gwas"]
    elif model == "mvlmm":
        cmd += ["--pheno", pheno_file, "--mpheno-list", "1,2",
                "--multivariate"]
    elif model == "repeatability":
        cmd += ["--pheno", pheno_file, "--mpheno", "1",
                "--repeatability-model"]
    elif model == "rrm":
        cmd += ["--long-pheno", pheno_file, "--tpheno", "1", "--mpheno", "2",
                "--rrm-order", "2", "--gwas"]
    elif model == "maternal":
        # 母体模型：需要 --maternal-id 指定母畜列、--qcovar 提供协变量
        # maternal_id_col=3 指向第 3 列（DAM_ID）
        covar_file = str(Path(pheno_file).parent / "covar.txt")
        cmd += ["--pheno", pheno_file, "--mpheno", "1",
                "--maternal-model", "--maternal-id", "3",
                "--qcovar", covar_file]
    elif model == "threshold":
        # 二元性状：probit 连接，启用 SPA 校准
        cmd += ["--pheno", pheno_file, "--mpheno", "1",
                "--threshold-model", "--link", "probit", "--spa"]
    elif model == "gxe":
        # GxE：表型格式 FID IID ENV VALUE
        # --mpheno N 读取 FID IID 之后的第 N 列：mpheno=1→ENV, mpheno=2→VALUE
        # --env-covar 3 读取第 3 列（ENV）作为环境协变量
        cmd += ["--pheno", pheno_file, "--mpheno", "2",
                "--gxe-model", "--env-covar", "3"]
    elif model == "survival":
        # Survival：表型格式 FID IID TIME EVENT
        # time_col=3, event_col=4（默认值，无需额外指定）
        cmd += ["--pheno", pheno_file, "--survival-model"]
    elif model == "count_poisson":
        # Poisson PSR：FID IID COUNT
        cmd += ["--pheno", pheno_file, "--mpheno", "1",
                "--count-model", "--family", "poisson"]
    elif model == "count_binomial":
        # Binomial PSR：--pheno 事件率 + --binomial-pheno FID IID EVENT TRIAL
        cmd += ["--pheno", pheno_file, "--mpheno", "1",
                "--count-model", "--family", "binomial",
                "--binomial-pheno", binom_file]
    else:
        raise ValueError(f"Unknown model: {model}")

    if extra_args:
        cmd += extra_args

    # 清理 LD_LIBRARY_PATH，避免 Trae 沙箱的旧 libstdc++ 被抢先加载
    clean_env = dict(os.environ)
    traec_paths = [p for p in clean_env.get("LD_LIBRARY_PATH", "").split(":")
                   if "trae-cn-server" not in p]
    clean_env["LD_LIBRARY_PATH"] = ":".join(traec_paths)

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600,
                            env=clean_env)
    if result.returncode != 0:
        print(f"  STDERR: {result.stderr[-1500:]}")
        raise RuntimeError(f"CosmicScanner failed for {model}: {result.returncode}")

    # 确定输出文件
    if model == "mvlmm":
        return out_prefix + ".assoc.mv.txt"
    elif model == "rrm":
        return out_prefix + ".L0.assoc.txt"
    elif model == "gxe":
        return out_prefix + ".gxe.assoc.txt"
    return out_prefix + ".assoc.txt"


# ============================================================
# 5. 结果解析
# ============================================================

def parse_assoc_results(assoc_file, model):
    """解析关联结果，返回 p-values 数组"""
    pvals = []
    # 不同模型的 P 值列名
    p_col_map = {
        "lmm": "P_A",        # 标准加性 LMM
        "gxe": "JOINT_P",      # GxE 2-df 联合检验
        "rrm": "P",            # L0 单点检验（已选 .L0.assoc.txt）
    }
    p_col = p_col_map.get(model, "P")
    with open(assoc_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                p = float(row[p_col])
                if math.isfinite(p) and 0 < p <= 1:
                    pvals.append(p)
            except (KeyError, ValueError):
                continue
    return np.array(pvals)


def compute_metrics(pvals, df=1):
    """计算 λ_GC、type-I error、KS 检验

    λ_GC (基因组膨胀因子) 定义：
        λ_GC = median(χ²_obs) / median(χ²_df)
    其中 χ²_obs 是 GWAS 检验统计量，由 P 反推 χ²_obs = F^{-1}_{χ²(df)}(1 - P)。
    若 P~U(0,1)（零假设成立），χ²_obs 的中位数 = χ²(df) 的中位数，λ_GC = 1.0。

    df: 检验自由度（1=标准 Wald/t; 2=GxE 2-df 联合检验）
    """
    n = len(pvals)
    if n == 0:
        return {"lambda_gc": float("nan"), "n_snps": 0}

    # 各自由度下 χ² 的中位数
    chi2_median = {1: 0.4549364, 2: 1.3862944, 3: 2.3659715, 4: 3.3566940}

    # 由 P 反推 χ²(df) 检验统计量：χ² = F^{-1}(1 - P)
    try:
        from scipy.stats import chi2 as chi2_dist
        chi2_obs = chi2_dist.ppf(1.0 - pvals, df=df)
        med_ref = float(chi2_dist.ppf(0.5, df=df))
    except ImportError:
        if df == 1:
            # Acklam 逆正态 CDF → χ²(1) = z²
            def inv_norm_cdf(p):
                a = [-3.969683028665376e+01, 2.209460984245205e+02,
                     -2.759285104469687e+02, 1.383577518672690e+02,
                     -3.066479806614716e+01, 2.506628277459239e+00]
                b = [-5.447609879822406e+01, 1.615858368580409e+02,
                     -1.556989798598866e+02, 6.680131188771972e+01,
                     -1.328068155288572e+01]
                c = [-7.784894002430293e-03, -3.223964580411365e-01,
                     -2.400758277161838e+00, -2.549732539343734e+00,
                     4.374664141464968e+00, 2.938163982698783e+00]
                d = [7.784695709041462e-03, 3.224671290700398e-01,
                     2.445134137142996e+00, 3.754408661907416e+00]
                plow = 0.02425; phigh = 1 - plow
                if p < plow:
                    q = math.sqrt(-2 * math.log(p))
                    x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                        ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
                elif p <= phigh:
                    q = p - 0.5; r = q*q
                    x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
                        (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
                else:
                    q = math.sqrt(-2 * math.log(1-p))
                    x = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                         ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
                return x
            chi2_obs = np.array([inv_norm_cdf(1 - p/2)**2 for p in pvals])
        elif df == 2:
            # χ²(2) = -2 ln(P)，中位数 = 2 ln(2) = 1.3863
            chi2_obs = -2.0 * np.log(pvals)
        else:
            # 退化：用 χ²(2) 近似
            chi2_obs = -2.0 * np.log(pvals)
        med_ref = chi2_median.get(df, 1.3862944)

    chi2_obs = chi2_obs[np.isfinite(chi2_obs)]
    lambda_gc = float(np.median(chi2_obs)) / med_ref if len(chi2_obs) > 0 else float("nan")

    # Type-I error 在各 α 水平
    alphas = [0.05, 0.01, 0.001, 5e-8]
    type1 = {}
    for a in alphas:
        type1[f"type1_{a}"] = float(np.mean(pvals < a))

    # KS 检验 vs U(0,1)
    ks_stat, ks_p = ks_test_uniform(pvals)

    return {
        "lambda_gc": lambda_gc,
        "n_snps": n,
        "ks_stat": float(ks_stat),
        "ks_p": float(ks_p),
        **type1,
    }


def ks_test_uniform(pvals):
    """简化版 KS 检验 vs U(0,1)"""
    n = len(pvals)
    if n == 0:
        return float("nan"), float("nan")
    sorted_p = np.sort(pvals)
    # CDF 理论值 = 排序后的 p 值本身（U(0,1) 的 CDF = x）
    # 经验 CDF = i/n
    emp_upper = np.arange(1, n + 1) / n
    emp_lower = np.arange(0, n) / n
    d_plus = np.max(emp_upper - sorted_p)
    d_minus = np.max(sorted_p - emp_lower)
    d = max(d_plus, d_minus)
    # KS p-value 近似
    en = math.sqrt(n)
    p_val = 0.0
    # Kolmogorov 分布级数
    for j in range(1, 101):
        term = 2.0 * ((-1) ** (j - 1)) * math.exp(-2.0 * j * j * (d * en) ** 2)
        p_val += term
        if abs(term) < 1e-10:
            break
    p_val = max(0.0, min(1.0, p_val))
    return d, p_val


# ============================================================
# 6. 报告生成
# ============================================================

def generate_report(results, outdir, models, n_replicates):
    """生成 Markdown 校准报告 + CSV 汇总"""
    outdir = Path(outdir)

    # CSV 汇总
    csv_path = outdir / "type1_summary.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["model", "rep", "lambda_gc", "n_snps", "ks_stat", "ks_p",
                         "type1_0.05", "type1_0.01", "type1_0.001", "type1_5e-8"])
        for model in models:
            for rep in range(1, n_replicates + 1):
                key = (model, rep)
                if key not in results:
                    continue
                m = results[key]
                writer.writerow([model, rep,
                                 f"{m['lambda_gc']:.6f}", m["n_snps"],
                                 f"{m['ks_stat']:.6f}", f"{m['ks_p']:.6f}",
                                 f"{m.get('type1_0.05', float('nan')):.6f}",
                                 f"{m.get('type1_0.01', float('nan')):.6f}",
                                 f"{m.get('type1_0.001', float('nan')):.6f}",
                                 f"{m.get('type1_5e-8', float('nan')):.8f}"])

    # QQ 数据
    qq_path = outdir / "qq_data.csv"
    with open(qq_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["model", "rep", "rank", "n", "expected_log10p", "observed_log10p"])
        for model in models:
            for rep in range(1, n_replicates + 1):
                key = (model, rep)
                if key not in results or "pvals" not in results[key]:
                    continue
                pvals = np.sort(results[key]["pvals"])
                n = len(pvals)
                for i, p in enumerate(pvals):
                    exp_p = (i + 0.5) / n
                    f.write(f"{model},{rep},{i+1},{n},"
                            f"{-math.log10(exp_p):.6f},{-math.log10(max(p,1e-300)):.6f}\n")

    # Markdown 报告
    md_path = outdir / "type1_calibration_report.md"
    with open(md_path, "w") as f:
        f.write("# Type-I Error 校准报告\n\n")
        f.write(f"**模型**: {', '.join(models)}\n")
        f.write(f"**重复次数**: {n_replicates}\n\n")

        f.write("## 1. λ_GC 基因组膨胀因子\n\n")
        f.write("理想值 = 1.0。判定：PASS [0.95, 1.05] / WARN [0.90, 1.10] / FAIL 其他\n\n")
        f.write("| Model | Mean | SD | Min | Max | Status |\n")
        f.write("|-------|-----:|---:|----:|----:|:------:|\n")
        for model in models:
            lambdas = [results[(model, r)]["lambda_gc"]
                       for r in range(1, n_replicates + 1)
                       if (model, r) in results]
            if not lambdas:
                continue
            lambdas = np.array(lambdas)
            mean = float(np.nanmean(lambdas))
            sd = float(np.nanstd(lambdas))
            mn = float(np.nanmin(lambdas))
            mx = float(np.nanmax(lambdas))
            if 0.95 <= mean <= 1.05:
                status = "PASS"
            elif 0.90 <= mean <= 1.10:
                status = "WARN"
            else:
                status = "FAIL"
            f.write(f"| {model} | {mean:.4f} | {sd:.4f} | {mn:.4f} | {mx:.4f} | {status} |\n")

        f.write("\n## 2. Type-I Error Rate\n\n")
        f.write("理想值 = α。判定：PASS ≤ 2α / WARN (2α, 3α] / FAIL > 3α\n\n")
        for alpha in [0.05, 0.01, 0.001, 5e-8]:
            f.write(f"### α = {alpha}\n\n")
            f.write("| Model | Mean | SE | Expected | Ratio | Status |\n")
            f.write("|-------|-----:|---:|---------:|------:|:------:|\n")
            for model in models:
                rates = [results[(model, r)].get(f"type1_{alpha}", float("nan"))
                         for r in range(1, n_replicates + 1)
                         if (model, r) in results]
                if not rates:
                    continue
                rates = np.array(rates)
                mean = float(np.nanmean(rates))
                se = float(np.nanstd(rates) / math.sqrt(len(rates)))
                ratio = mean / alpha if alpha > 0 else float("nan")
                if mean <= 2 * alpha:
                    status = "PASS"
                elif mean <= 3 * alpha:
                    status = "WARN"
                else:
                    status = "FAIL"
                f.write(f"| {model} | {mean:.6f} | {se:.6f} | {alpha:.2e} | "
                        f"{ratio:.3f} | {status} |\n")
            f.write("\n")

        f.write("## 3. KS 检验（p-value 均匀性）\n\n")
        f.write("H0: p-values ~ U(0,1)。P > 0.05 表示不拒绝均匀性\n\n")
        f.write("| Model | Mean KS-P | Min KS-P | % P>0.05 |\n")
        f.write("|-------|----------:|---------:|---------:|\n")
        for model in models:
            ks_ps = [results[(model, r)].get("ks_p", float("nan"))
                     for r in range(1, n_replicates + 1)
                     if (model, r) in results]
            if not ks_ps:
                continue
            ks_ps = np.array(ks_ps)
            mean_p = float(np.nanmean(ks_ps))
            min_p = float(np.nanmin(ks_ps))
            pct = float(np.nanmean(ks_ps > 0.05))
            f.write(f"| {model} | {mean_p:.4f} | {min_p:.4f} | {pct:.2%} |\n")

        f.write("\n## 4. 总结\n\n")
        all_pass = True
        for model in models:
            lambdas = [results[(model, r)]["lambda_gc"]
                       for r in range(1, n_replicates + 1)
                       if (model, r) in results]
            if not lambdas:
                f.write(f"- **{model}**: 无数据\n")
                continue
            mean_lambda = float(np.nanmean(lambdas))
            lambda_ok = 0.95 <= mean_lambda <= 1.05

            type1_ok = True
            for alpha in [0.05, 0.01, 0.001]:
                rates = [results[(model, r)].get(f"type1_{alpha}", 0)
                         for r in range(1, n_replicates + 1)
                         if (model, r) in results]
                if rates and np.nanmean(rates) > 2 * alpha:
                    type1_ok = False
                    break

            status = "PASS" if (lambda_ok and type1_ok) else "FAIL"
            if status == "FAIL":
                all_pass = False
            f.write(f"- **{model}**: λ_GC={mean_lambda:.4f} {'✓' if lambda_ok else '✗'}  "
                    f"Type-I {'✓' if type1_ok else '✗'}  → **{status}**\n")

        f.write(f"\n**总体判定**: {'✓ ALL PASS' if all_pass else '✗ 需检查'}\n")
        f.write(f"\nQQ 数据见: {qq_path}\n")
        f.write(f"原始汇总见: {csv_path}\n")

    return str(md_path), str(csv_path), str(qq_path)


# ============================================================
# 7. 主流程
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="CosmicScanner Type-I Error Calibration")
    parser.add_argument("--bin", required=True, help="Path to cosmicscanner binary")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--n-animals", type=int, default=500)
    parser.add_argument("--n-snps", type=int, default=1000)
    parser.add_argument("--n-replicates", type=int, default=10)
    parser.add_argument("--models", default="lmm,mvlmm,repeatability,rrm",
                        help="Comma-separated: lmm,mvlmm,repeatability,rrm")
    parser.add_argument("--thread-num", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20260626)
    args = parser.parse_args()

    models = [m.strip() for m in args.models.split(",")]
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"=== Type-I Error 校准 ===")
    print(f"  Models: {models}")
    print(f"  N animals: {args.n_animals}, N SNPs: {args.n_snps}")
    print(f"  Replicates: {args.n_replicates}")
    print(f"  Output: {outdir}\n")

    results = {}

    for rep in range(1, args.n_replicates + 1):
        print(f"--- Replicate {rep}/{args.n_replicates} ---")
        rep_dir = outdir / f"sim_{rep}"
        rep_dir.mkdir(parents=True, exist_ok=True)

        # 1. 基因型（每个 replicate 重新生成）
        seed = args.seed + rep * 100
        print("  [1] 模拟基因型...")
        genotypes, mafs = simulate_genotypes(args.n_animals, args.n_snps, seed=seed)
        plink_prefix = write_plink_files(rep_dir, genotypes, mafs, args.n_animals)

        # 2. GRM
        print("  [2] 构建 GRM...")
        grm = build_grm(genotypes)

        for model in models:
            print(f"  [3] 运行 {model}...")
            binom_file = None  # 仅 count_binomial 模型会设置
            try:
                # 3. 零假设表型
                if model == "lmm":
                    y = simulate_null_pheno_lmm(grm, seed=seed + 1)
                    pheno_file = write_pheno_lmm(rep_dir, y, args.n_animals)
                elif model == "mvlmm":
                    y1, y2 = simulate_null_pheno_mvlmm(grm, seed=seed + 2)
                    pheno_file = write_pheno_mvlmm(rep_dir, y1, y2, args.n_animals)
                elif model == "repeatability":
                    records = simulate_null_pheno_repeatability(grm, seed=seed + 3)
                    pheno_file = write_pheno_repeatability(rep_dir, records, args.n_animals)
                elif model == "rrm":
                    records = simulate_null_pheno_rrm(grm, seed=seed + 4)
                    pheno_file = write_pheno_rrm(rep_dir, records, args.n_animals)
                elif model == "maternal":
                    records = simulate_null_pheno_maternal(grm, seed=seed + 5)
                    pheno_file = write_pheno_maternal(rep_dir, records, args.n_animals)
                    # 母体模型需要协变量文件
                    write_covar_file(rep_dir, args.n_animals)
                elif model == "threshold":
                    y_bin = simulate_null_pheno_threshold(grm, prevalence=0.3, seed=seed + 6)
                    pheno_file = write_pheno_threshold(rep_dir, y_bin, args.n_animals)
                elif model == "gxe":
                    env, y = simulate_null_pheno_gxe(grm, seed=seed + 7)
                    pheno_file = write_pheno_gxe(rep_dir, env, y, args.n_animals)
                elif model == "survival":
                    s_time, s_event = simulate_null_pheno_survival(grm, event_rate=0.6, seed=seed + 8)
                    pheno_file = write_pheno_survival(rep_dir, s_time, s_event, args.n_animals)
                elif model == "count_poisson":
                    y_cnt, _ = simulate_null_pheno_count(grm, family="poisson",
                                                         mean_count=5.0, seed=seed + 9)
                    pheno_file, _ = write_pheno_count(rep_dir, y_cnt, None,
                                                      args.n_animals, family="poisson")
                elif model == "count_binomial":
                    y_rate, trials = simulate_null_pheno_count(
                        grm, family="binomial", binomial_trials=10,
                        binomial_rate=0.3, seed=seed + 10)
                    pheno_file, binom_file = write_pheno_count(
                        rep_dir, y_rate, trials, args.n_animals, family="binomial")
                else:
                    print(f"    跳过未知模型: {model}")
                    continue

                # 4. 运行 CosmicScanner
                out_prefix = str(rep_dir / f"cosmic_{model}")
                assoc_file = run_cosmicscanner(
                    args.bin, plink_prefix, model, pheno_file,
                    out_prefix, thread_num=args.thread_num,
                    binom_file=binom_file if model == "count_binomial" else None)

                # 5. 解析结果（GxE 2-df 检验用 df=2 计算 λ_GC）
                pvals = parse_assoc_results(assoc_file, model)
                df_test = 2 if model == "gxe" else 1
                metrics = compute_metrics(pvals, df=df_test)
                metrics["pvals"] = pvals
                results[(model, rep)] = metrics

                print(f"    λ_GC={metrics['lambda_gc']:.4f}  "
                      f"KS-P={metrics.get('ks_p', float('nan')):.4f}  "
                      f"type1(0.05)={metrics.get('type1_0.05', float('nan')):.4f}  "
                      f"N={metrics['n_snps']}")

            except Exception as e:
                print(f"    ERROR: {e}")
                results[(model, rep)] = {
                    "lambda_gc": float("nan"), "n_snps": 0,
                    "ks_stat": float("nan"), "ks_p": float("nan"),
                    "type1_0.05": float("nan"), "type1_0.01": float("nan"),
                    "type1_0.001": float("nan"), "type1_5e-8": float("nan"),
                }

    print("\n=== 生成报告 ===")
    md_path, csv_path, qq_path = generate_report(
        results, outdir, models, args.n_replicates)
    print(f"报告: {md_path}")
    print(f"汇总: {csv_path}")
    print(f"QQ 数据: {qq_path}")


if __name__ == "__main__":
    main()
