#!/usr/bin/env python3
"""
CosmicScanner --maternal-snp-test 模拟验证套件
================================================
对指定场景模拟数据，调用 :
  --maternal-model --maternal-snp-test direct | maternal | joint
验证：
  A. Null（beta_d=0, beta_m=0）—— direct/maternal/joint 均无系统膨胀
  B. Direct-only causal           —— direct 强 / maternal 弱 / joint 强
  C. Maternal-only causal         —— 反过来
  D. Direct+maternal 同向
  E. Direct+maternal 反向（joint 发现 marginal 被削弱的关联）
  F. random direct-maternal covariance（sigma_dm != 0）不产生虚假 SNP 信号
  G. Missing dam genotype (10/30/50%) 安全处理与 usable N 报告
  H. near-identical offspring/dam SNP 向量 -> joint rank/conditioning 诊断

用法：
  python3 run_maternal_snp_effect_simulation.py --bin /path/cosmicscanner --out dir [--n 250] [--m 1200] [--seed 7]

输出：out/summary.tsv（每场景 PASS/FAIL + 指标）
退出码：全部 PASS = 0，否则 1
"""
import argparse
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np


# ============================================================
# 1. 数据生成（含孟德尔遗传：offspring 基因型由父/母单倍型抽样得到）
# ============================================================

DOSE_TO_CODE = {0: 0b11, 1: 0b10, 2: 0b00}  # 读入器: 00->2.0 (A1同型), 10->1.0, 11->0.0


def generate_subjects(n_founders, n_off, n_snps, seed):
    """生成全部基因型个体 (founder 亲本 + offspring 子代)，孟德尔遗传。

    返回 dict:
      geno: (n_total, n_snps) int8 dose（A1 alt 等位基因拷贝数 0/1/2）
      mafs: (n_snps,)
      ov:   每个 offspring 的 dam 在 subject 列表中的索引
    """
    rng = np.random.default_rng(seed)
    mafs = rng.uniform(0.05, 0.5, n_snps)

    # founder 每个位点的两条单倍型：等位基因 ~ Bernoulli(p_j)
    founder_hap = rng.binomial(1, mafs, size=(2, n_founders, n_snps))  # (2,nf,m)
    founder_geno = founder_hap[0] + founder_hap[1]  # (nf, m) dose

    # 每个 offspring: dam = 随机 founder, sire = 随机 founder
    dam_idx = rng.integers(0, n_founders, n_off)
    sire_idx = rng.integers(0, n_founders, n_off)

    off_geno = np.zeros((n_off, n_snps), dtype=np.int8)
    for i in range(n_off):
        d = dam_idx[i]
        s = sire_idx[i]
        g_d = rng.integers(0, 2, 2)  # 从 dam 两条单倍型中抽 1 条
        g_s = rng.integers(0, 2, 2)  # 从 sire 两条单倍型中抽 1 条
        off_hap0 = founder_hap[g_d[0], d]
        off_hap1 = founder_hap[g_s[0], s]
        off_geno[i] = off_hap0 + off_hap1

    geno = np.vstack([founder_geno, off_geno]).astype(np.int8)
    return {
        "geno": geno,
        "mafs": mafs,
        "n_founders": n_founders,
        "n_off": n_off,
        "dam_idx": dam_idx,
    }


def build_grm(genotypes, mafs):
    M = genotypes.astype(np.float64)
    Z = M - 2.0 * mafs[np.newaxis, :]
    S = 2.0 * np.sum(mafs * (1 - mafs))
    return Z @ Z.T / S


def simulate_pheno(data, grm, causal, effects, vc, seed):
    """y 记录级（每个 offspring 1 条记录）：
    y_r = beta_d * g_off[r,c] + beta_m * g_dam[r,c]
          + a[off_r] + m[dam_r] + Xb + e_r
    a ~ N(0, sigma_ad^2 G), m ~ N(0, sigma_am^2 G)，二者按 sigma_adm 相关。
    """
    rng = np.random.default_rng(seed + 2024)
    n_total = data["geno"].shape[0]
    n_off = data["n_off"]
    n_snps = data["geno"].shape[1]
    L = np.linalg.cholesky(grm + 1e-6 * np.eye(n_total))
    z = rng.standard_normal(n_total)
    a = (L @ z) * math.sqrt(vc["sigma2_ad"])

    if vc["sigma2_am"] > 0:
        if abs(vc.get("sigma_adm", 0.0)) > 0 and vc["sigma2_ad"] > 0 and vc["sigma2_am"] > 0:
            rho = np.clip(vc["sigma_adm"] / math.sqrt(vc["sigma2_ad"] * vc["sigma2_am"] + 1e-12),
                          -0.99, 0.99)
            z2 = rng.standard_normal(n_total)
            m = rho * (a / math.sqrt(vc["sigma2_ad"] + 1e-12)) * math.sqrt(vc["sigma2_am"]) \
                + math.sqrt(max(1 - rho ** 2, 0.0)) * z2 * math.sqrt(vc["sigma2_am"])
        else:
            z2 = rng.standard_normal(n_total)
            m = (L @ z2) * math.sqrt(vc["sigma2_am"])
    else:
        m = np.zeros(n_total)

    off_idx = np.arange(data["n_founders"], n_total)
    # dam 的 subject 索引就是 founder 内的索引（founder 排在最前）
    beta_d, beta_m = effects
    causal_off = data["geno"][off_idx, causal]
    causal_dam = data["geno"][data["dam_idx"], causal]
    y = beta_d * causal_off + beta_m * causal_dam + a[off_idx] + m[data["dam_idx"]]
    e = rng.standard_normal(n_off) * math.sqrt(vc["sigma2_e"])
    y = y + e
    y = y - y.mean()
    return y, off_idx, data["dam_idx"]


def write_plink(outdir, data, prefix="plink"):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = outdir / prefix
    n, m = data["geno"].shape
    with open(str(prefix) + ".fam", "w") as f:
        for i in range(n):
            f.write(f"0\tIND{i+1}\t0\t0\t1\t-9\n")
    with open(str(prefix) + ".bim", "w") as f:
        for j in range(m):
            maf = data["mafs"][j]
            a1 = "A" if maf <= 0.5 else "T"
            a2 = "T" if a1 == "A" else "A"
            f.write(f"1\trs{j+1}\t0\t{(j+1)*100}\t{a1}\t{a2}\n")
    nbytes = (n + 3) // 4
    bed = np.zeros((m, nbytes), dtype=np.uint8)
    for j in range(m):
        for i in range(n):
            d = int(data["geno"][i, j])
            code = DOSE_TO_CODE.get(d, 0b01)
            bi = i // 4
            off = (i % 4) * 2
            bed[j, bi] |= (code << off)
    with open(str(prefix) + ".bed", "wb") as f:
        f.write(b"\x6c\x1b\x01")
        f.write(bed.tobytes())
    return str(prefix)


def write_pheno(outdir, off_idx, dam_idx, y, missing_frac=0.0, self_dam=False, seed=0):
    """off_idx: subject 索引；dam_idx: founder subject 索引"""
    outdir = Path(outdir)
    rng = np.random.default_rng(seed + 99)
    pheno = outdir / "pheno.txt"
    n_off = len(off_idx)
    rows = []
    for r in range(n_off):
        iid = f"IND{off_idx[r]+1}"
        di = dam_idx[r]
        if self_dam:
            dam_id = iid
        elif missing_frac > 0 and rng.random() < missing_frac:
            dam_id = f"MISS{r}"
        else:
            dam_id = f"IND{di+1}"
        rows.append(f"0\t{iid}\t{dam_id}\t{y[r]:.6f}")
    with open(pheno, "w") as f:
        f.write("\n".join(rows) + "\n")
    return str(pheno)


def write_covar(outdir, n_total, seed=0):
    outdir = Path(outdir)
    rng = np.random.default_rng(seed + 7)
    cov = outdir / "covar.txt"
    with open(cov, "w") as f:
        for i in range(n_total):
            f.write(f"0\tIND{i+1}\t{rng.integers(0, 2)}\n")
    return str(cov)


# ============================================================
# 2. 运行 CosmicScanner
# ============================================================

def run_cosmic(binary, plink, pheno, covar, mode, out_prefix):
    cmd = [
        binary, "--bfile", plink, "--pheno", pheno, "--mpheno", "1", "--gwas",
        "--maternal-model", "--maternal-id", "3",
        "--maternal-snp-test", mode,
        "--qcovar", covar,
        "--no-loco",
        "--out", out_prefix,
    ]
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    r = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=3600)
    if r.returncode != 0:
        raise RuntimeError(f"cosmicscanner failed ({mode}):\n{r.stderr[-3000:]}")
    return out_prefix


def read_assoc(path):
    return np.genfromtxt(path, names=True, dtype=None, encoding="utf-8")


def read_null(path):
    with open(path) as f:
        header = f.readline().strip().split("\t")
        vals = f.readline().strip().split("\t")
    return dict(zip(header, [float(v) for v in vals]))


# ============================================================
# 3. 统计工具
# ============================================================

def lambda_gc(pvals):
    """Genomic-control lambda from 1-df p-values: median(qchisq(1-p, df=1)) / 0.4549."""
    from scipy.stats import chi2 as _chi2
    p = np.asarray(pvals, dtype=float)
    p = p[np.isfinite(p) & (p > 0)]
    p = np.clip(p, 1e-300, 1 - 1e-300)
    chi = _chi2.isf(p, 1)
    return float(np.median(chi) / 0.4549364) if len(p) else float("nan")


def qq_max_dev(pvals, alpha=0.05):
    p = np.asarray(pvals, dtype=float)
    p = np.sort(p[np.isfinite(p) & (p > 0)])
    n = len(p)
    if n == 0:
        return float("nan")
    expected = (np.arange(1, n + 1) - 0.5) / n
    return float(np.max(np.abs(p / alpha - expected / alpha)))

# ============================================================
# 4. 场景
# ============================================================

VC_DEFAULT = {"sigma2_ad": 0.25, "sigma2_am": 0.12, "sigma_adm": 0.0, "sigma2_e": 0.60}

SCENARIOS = ["A_null", "B_direct_only", "C_maternal_only", "D_same_dir",
             "E_opposite", "F_dm_cov", "G_missing10", "G_missing30", "G_missing50", "H_self_dam"]


def run_scenario(binary, outroot, name, n_off, n_snps, n_founders, seed):
    data = generate_subjects(n_founders, n_off, n_snps, seed)
    grm = build_grm(data["geno"], data["mafs"])
    vc = dict(VC_DEFAULT)
    causal = 17
    effects = (0.0, 0.0)
    missing = 0.0
    self_dam = False

    if name == "B_direct_only":
        effects = (0.8, 0.0)
    elif name == "C_maternal_only":
        effects = (0.0, 0.8)
    elif name == "D_same_dir":
        effects = (0.8, 0.5)
    elif name == "E_opposite":
        effects = (0.8, -1.2)
    elif name == "F_dm_cov":
        effects = (0.0, 0.0)
        vc["sigma_adm"] = -0.06
    elif name.startswith("G_missing"):
        missing = int(name.replace("G_missing", "")) / 100.0
        effects = (0.0, 0.8)
    elif name == "H_self_dam":
        effects = (0.0, 0.8)
        self_dam = True

    y, off_idx, dam_idx = simulate_pheno(data, grm, causal, effects, vc, seed)
    d = Path(outroot) / name
    d.mkdir(parents=True, exist_ok=True)
    plink = write_plink(d, data)
    pheno = write_pheno(d, off_idx, dam_idx, y, missing_frac=missing,
                        self_dam=self_dam, seed=seed)
    covar = write_covar(d, data["geno"].shape[0], seed=seed)

    # 运行三个模式
    out = {}
    out["direct"] = read_assoc(run_cosmic(binary, plink, pheno, covar, "direct", str(d / "direct")) + ".assoc.txt")
    out["maternal"] = read_assoc(run_cosmic(binary, plink, pheno, covar, "maternal", str(d / "maternal")) + ".maternal.assoc.txt")
    out["joint"] = read_assoc(run_cosmic(binary, plink, pheno, covar, "joint", str(d / "joint")) + ".joint.assoc.txt")

    snp_sid = np.char.strip(np.asarray(out["direct"]["SNP"], dtype=str))
    causal_sid = f"rs{causal+1}"

    def pcol(t, col):
        return np.asarray(t[col], dtype=float)

    res = {"scenario": name, "causal_snp": causal_sid}
    res["lambda_direct"] = lambda_gc(pcol(out["direct"], "P"))
    res["lambda_maternal"] = lambda_gc(pcol(out["maternal"], "P"))
    res["lambda_joint"] = lambda_gc(pcol(out["joint"], "P_JOINT"))

    # 因果位点信号
    def snp_row(t, sid, col):
        mask = snp_sid == sid
        if not mask.any():
            return None
        return float(np.asarray(t[mask][col], dtype=float)[0])

    direct_sig = snp_row(out["direct"], causal_sid, "P")
    maternal_sig = snp_row(out["maternal"], causal_sid, "P")
    joint_sig = snp_row(out["joint"], causal_sid, "P_JOINT")
    res["p_direct_causal"] = direct_sig
    res["p_maternal_causal"] = maternal_sig
    res["p_joint_causal"] = joint_sig

    # top-k recovery
    def topk_rank(t, sid, col, k=10):
        r = np.argsort(pcol(t, col))[:k]
        return int((snp_sid[r] == sid).sum())

    res["top10_direct"] = topk_rank(out["direct"], causal_sid, "P", 10)
    res["top10_maternal"] = topk_rank(out["maternal"], causal_sid, "P", 10)
    res["top10_joint"] = topk_rank(out["joint"], causal_sid, "P_JOINT", 10)

    # missing/self-dam 相关指标
    if name.startswith("G_missing"):
        meta = {}
        with open(str(d / "joint.maternal.snp.meta.txt")) as f:
            for line in f:
                if "\t" in line:
                    k, v = line.rstrip("\n").split("\t", 1)
                    meta[k.strip()] = v.strip()
        res["meta_n_total"] = int(meta["n_total"])
        res["meta_n_with_dam"] = int(meta["n_with_dam_genotype"])
        res["meta_n_dropped"] = int(meta["n_dropped_missing_dam"])
    elif name == "H_self_dam":
        statuses = np.asarray([str(s) for s in out["joint"]["STATUS"]])
        res["joint_status_counts"] = {}
        for s in sorted(set(statuses)):
            res["joint_status_counts"][s] = int((statuses == s).sum())
        ranks = pcol(out["joint"], "RANK")
        res["joint_rank1"] = int((ranks == 1).sum())
        res["joint_cond"] = float(np.nanmax(pcol(out["joint"], "COND")))

    return res


# ============================================================
# 5. 判定
# ============================================================

def evaluate(res):
    name = res["scenario"]
    checks = []

    if name == "A_null":
        checks.append(("null lambda_direct<1.35", res["lambda_direct"] < 1.35, res["lambda_direct"]))
        checks.append(("null lambda_maternal<1.35", res["lambda_maternal"] < 1.35, res["lambda_maternal"]))
        checks.append(("null lambda_joint<1.35", res["lambda_joint"] < 1.35, res["lambda_joint"]))
    elif name == "B_direct_only":
        l_d = -math.log10(res["p_direct_causal"] or 1e-300)
        l_m = -math.log10(res["p_maternal_causal"] or 1e-300)
        checks.append(("direct causal p<1e-3", (res["p_direct_causal"] or 1) < 1e-3, res["p_direct_causal"]))
        checks.append(("maternal weaker than direct (-log10)", l_m - l_d < -1.0, (l_m, l_d)))
        checks.append(("joint causal p<1e-3", (res["p_joint_causal"] or 1) < 1e-3, res["p_joint_causal"]))
        checks.append(("direct top10 hit", res["top10_direct"] >= 1, res["top10_direct"]))
    elif name == "C_maternal_only":
        l_m = -math.log10(res["p_maternal_causal"] or 1e-300)
        l_d = -math.log10(res["p_direct_causal"] or 1e-300)
        checks.append(("maternal causal p<1e-3", (res["p_maternal_causal"] or 1) < 1e-3, res["p_maternal_causal"]))
        checks.append(("direct weaker than maternal (-log10)", l_d - l_m < -1.0, (l_d, l_m)))
        checks.append(("joint causal p<1e-3", (res["p_joint_causal"] or 1) < 1e-3, res["p_joint_causal"]))
        checks.append(("maternal top10 hit", res["top10_maternal"] >= 1, res["top10_maternal"]))
    elif name == "D_same_dir":
        checks.append(("direct causal p<1e-3", (res["p_direct_causal"] or 1) < 1e-3, res["p_direct_causal"]))
        checks.append(("maternal causal p<1e-3", (res["p_maternal_causal"] or 1) < 1e-3, res["p_maternal_causal"]))
        checks.append(("joint causal p<1e-3", (res["p_joint_causal"] or 1) < 1e-3, res["p_joint_causal"]))
    elif name == "E_opposite":
        l_d = -math.log10(res["p_direct_causal"] or 1e-300)
        l_m = -math.log10(res["p_maternal_causal"] or 1e-300)
        l_j = -math.log10(res["p_joint_causal"] or 1e-300)
        checks.append(("joint causal p<1e-3", (res["p_joint_causal"] or 1) < 1e-3, res["p_joint_causal"]))
        checks.append(("joint strongest (-log10 joint > max marginals)",
                       l_j > max(l_d, l_m) + 0.5, (l_j, l_d, l_m)))
    elif name == "F_dm_cov":
        checks.append(("null(lambda, direct)<1.45 (cov not signal)", res["lambda_direct"] < 1.45, res["lambda_direct"]))
        checks.append(("null(lambda, maternal)<1.45", res["lambda_maternal"] < 1.45, res["lambda_maternal"]))
        checks.append(("null(lambda, joint)<1.45", res["lambda_joint"] < 1.45, res["lambda_joint"]))
    elif name.startswith("G_missing"):
        frac = int(name.replace("G_missing", "")) / 100.0
        checks.append(("meta reports dropped count", res["meta_n_dropped"] >= 1, res["meta_n_dropped"]))
        # 粗略: with_dam ≈ total*(1-frac)（允许采样波动）
        checks.append(("meta with-dam < total", res["meta_n_with_dam"] < res["meta_n_total"],
                       (res["meta_n_with_dam"], res["meta_n_total"])))
        checks.append(("maternal causal p<0.05", res["p_maternal_causal"] is None or res["p_maternal_causal"] < 0.05, res["p_maternal_causal"]))
    elif name == "H_self_dam":
        checks.append(("joint COLLINEAR flagged", res["joint_status_counts"].get("COLLINEAR", 0) > 0, res["joint_status_counts"]))
        checks.append(("joint (some) rank<2", res["joint_rank1"] > 0, res["joint_rank1"]))
    return checks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--out", default="maternal_snp_sim")
    ap.add_argument("--n", type=int, default=250, help="offspring 数")
    ap.add_argument("--m", type=int, default=1200, help="SNP 数")
    ap.add_argument("--nf", type=int, default=120, help="founder 数")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--scenarios", default="all")
    ap.add_argument("--thread-num", type=int, default=1)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    scen = SCENARIOS if args.scenarios == "all" else [s.strip() for s in args.scenarios.split(",")]

    results = []
    for name in scen:
        print(f"\n=== Scenario {name} ===", flush=True)
        res = run_scenario(args.bin, str(out), name, args.n, args.m, args.nf, args.seed)
        checks = evaluate(res)
        results.append((res, checks))
        for c, ok, val in checks:
            print("   [%s] %s (got %s)" % ("PASS" if ok else "FAIL", c, val))
        if not all(ok for _, ok, _ in checks):
            print("   -> scenario FAIL")

    # 汇总
    lines = ["scenario\tcausal_snp\tlambda_direct\tlambda_maternal\tlambda_joint\tp_direct\tp_maternal\tp_joint\ttop10_d\ttop10_m\ttop10_j\tverdict"]
    all_ok = True
    for res, checks in results:
        ok_all = all(ok for _, ok, _ in checks)
        all_ok &= ok_all
        lines.append("\t".join([
            res["scenario"], res["causal_snp"],
            f'{res["lambda_direct"]:.3f}', f'{res["lambda_maternal"]:.3f}', f'{res["lambda_joint"]:.3f}',
            f'{res["p_direct_causal"]}', f'{res["p_maternal_causal"]}', f'{res["p_joint_causal"]}',
            str(res["top10_direct"]), str(res["top10_maternal"]), str(res["top10_joint"]),
            "PASS" if ok_all else "FAIL",
        ]))
    (out / "summary.tsv").write_text("\n".join(lines) + "\n")
    print("\n" + "\n".join(lines))
    print("\nOVERALL:", "PASS" if all_ok else "FAIL")
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
