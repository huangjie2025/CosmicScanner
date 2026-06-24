import os
import subprocess
import math

def run_cmd(cmd, env_vars=None):
    print(f"Running: {cmd}")
    env = os.environ.copy()
    if env_vars:
        env.update(env_vars)
    res = subprocess.run(cmd, shell=True, env=env, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"Error running command: {cmd}")
        print(res.stderr)
        exit(1)
    return res

def pearsonr(x, y):
    n = len(x)
    if n == 0: return 0
    sum_x = sum(x)
    sum_y = sum(y)
    sum_x2 = sum(xi*xi for xi in x)
    sum_y2 = sum(yi*yi for yi in y)
    sum_xy = sum(xi*yi for xi, yi in zip(x, y))
    
    numerator = n * sum_xy - sum_x * sum_y
    denominator = math.sqrt((n * sum_x2 - sum_x**2) * (n * sum_y2 - sum_y**2))
    if denominator == 0: return 0
    return numerator / denominator

def main():
    data_dir = "/disk228/huangjie/PanOmics/CosmicScanner/examples"
    res_dir = "/disk228/huangjie/PanOmics/CosmicScanner/tests/results/raw_vs_gemma"
    grm_dir = "/disk228/huangjie/PanOmics/CosmicScanner/tests/results/grm"
    bfile = f"{data_dir}/mouse_hs1940/mouse_hs1940"
    
    conda_env_bin = "/disk228/huangjie/software/miniconda3/envs/benchmark/bin"
    gemma_bin = f"{conda_env_bin}/gemma"
    cosmic_bin = "/disk228/huangjie/PanOmics/Cosmic/build/CosmicScanner/cosmicscanner"

    print("=== Step 1: Running GEMMA ===")
    os.chdir(res_dir)
    run_cmd(f"{gemma_bin} -bfile {bfile} -gk 1 -o mouse")
    run_cmd(f"{gemma_bin} -bfile {bfile} -k output/mouse.cXX.txt -lmm 4 -o gemma_lmm")

    print("=== Step 2: Running CosmicScanner ===")
    run_cmd(f"{cosmic_bin} --bfile {bfile} --pheno {bfile}.pheno --gwas --no-loco --out {res_dir}/cosmic_lmm")

    print("=== Step 3: Comparing Results ===")
    gemma_out = f"{res_dir}/output/gemma_lmm.assoc.txt"
    cosmic_out = f"{res_dir}/cosmic_lmm.assoc.txt"

    if not os.path.exists(gemma_out) or not os.path.exists(cosmic_out):
        print("Output files not found!")
        exit(1)

    gemma_data = {}
    with open(gemma_out, 'r') as f:
        header = next(f).strip().split()
        idx_rs = header.index('rs')
        idx_beta = header.index('beta')
        idx_se = header.index('se')
        idx_p = header.index('p_wald')
        
        for line in f:
            parts = line.strip().split()
            if len(parts) > idx_p:
                rs = parts[idx_rs]
                if parts[idx_beta] != 'NA' and parts[idx_p] != 'NA':
                    gemma_data[rs] = {
                        'beta': float(parts[idx_beta]),
                        'se': float(parts[idx_se]),
                        'p': float(parts[idx_p])
                    }

    matched_beta, matched_se, matched_p = [], [], []
    g_beta, g_se, g_p = [], [], []

    with open(cosmic_out, 'r') as f:
        header = next(f).strip().split('\t')
        idx_snp = header.index('SNP')
        idx_beta = header.index('BETA')
        idx_se = header.index('SE')
        idx_logp = header.index('LOG10P') if 'LOG10P' in header else header.index('P')
        is_log = 'LOG10P' in header
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > idx_logp:
                snp = parts[idx_snp]
                if snp in gemma_data and parts[idx_beta] != 'nan':
                    c_beta = float(parts[idx_beta])
                    c_se = float(parts[idx_se])
                    if is_log:
                        c_logp = float(parts[idx_logp])
                    else:
                        pval = float(parts[idx_logp])
                        c_logp = -math.log10(pval + 1e-300)
                    
                    g = gemma_data[snp]
                    g_logp = -math.log10(g['p'] + 1e-300)
                    
                    matched_beta.append(c_beta)
                    matched_se.append(c_se)
                    matched_p.append(c_logp)
                    
                    g_beta.append(g['beta'])
                    g_se.append(g['se'])
                    g_p.append(g_logp)

    print(f"Matched {len(matched_beta)} SNPs.")

    if len(matched_beta) > 0:
        beta_r2 = pearsonr(matched_beta, g_beta) ** 2
        se_r2 = pearsonr(matched_se, g_se) ** 2
        p_r2 = pearsonr(matched_p, g_p) ** 2

        print(f"Beta R^2:      {beta_r2:.6f}")
        print(f"SE R^2:        {se_r2:.6f}")
        print(f"-log10(P) R^2: {p_r2:.6f}")
        
        if beta_r2 > 0.99 and p_r2 > 0.99:
            print("SUCCESS: CosmicScanner perfectly matches GEMMA!")
        else:
            print("WARNING: Discrepancy detected. Might be due to LOCO vs non-LOCO.")

    # Clean up intermediate files
    print("=== Step 4: Cleaning Up ===")
    for f in [f"{res_dir}/output/mouse.cXX.txt", f"{res_dir}/output/mouse.log.txt", f"{res_dir}/output/gemma_lmm.log.txt"]:
        if os.path.exists(f):
            os.remove(f)
    print(f"Cleanup complete. You can inspect the results at:\n - GEMMA: {gemma_out}\n - CosmicScanner: {cosmic_out}")

if __name__ == "__main__":
    main()

