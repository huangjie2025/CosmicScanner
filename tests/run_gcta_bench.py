import os, subprocess, math

def run_cmd(cmd):
    print(f"Running: {cmd}")
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"Error: {res.stderr}")
        exit(1)

def pearsonr(x, y):
    n = len(x)
    if n == 0: return 0
    sum_x = sum(x)
    sum_y = sum(y)
    sum_x2 = sum(xi*xi for xi in x)
    sum_y2 = sum(yi*yi for yi in y)
    sum_xy = sum(xi*yi for xi, yi in zip(x, y))
    numerator = n * sum_xy - sum_x * sum_y
    denominator = math.sqrt(max(0, (n * sum_x2 - sum_x**2) * (n * sum_y2 - sum_y**2)))
    if denominator == 0: return 0
    return numerator / denominator

root = os.environ.get("COSMICSCANNER_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
res_dir = os.environ.get("BENCHMARK_RESULT_DIR", os.path.join(root, "tests/results"))
data_dir = os.environ.get("BENCHMARK_DATA_DIR", os.path.join(root, "examples"))
gcta_bin = os.environ.get("GCTA_BIN", "gcta64")

os.chdir(res_dir)
# 1. Generate GCTA GRM
run_cmd(f"{gcta_bin} --bfile {data_dir}/mouse_hs1940/mouse_hs1940 --make-grm --out mouse_gcta")
# 2. Run GCTA MLMA (Mixed Linear Model Association)
run_cmd(f"{gcta_bin} --bfile {data_dir}/mouse_hs1940/mouse_hs1940 --grm mouse_gcta --pheno {data_dir}/mouse_hs1940/mouse_hs1940.pheno --mlma --out gcta_mlma")

# Read GCTA output
gcta_data = {}
with open("gcta_mlma.mlma") as f:
    header = next(f).strip().split()
    for line in f:
        p = line.strip().split()
        if len(p) > 8 and p[6] != 'NA' and p[6] != 'nan':
            # GCTA format: Chr SNP bp A1 A2 Freq b se p
            gcta_data[p[1]] = (float(p[6]), float(p[7]), float(p[8]))

# Read CosmicScanner output
cosmic_data = {}
with open("cosmic_lmm.assoc.txt") as f:
    header = next(f).strip().split('\t')
    idx_snp = header.index('SNP')
    idx_beta = header.index('BETA_A')
    idx_se = header.index('SE_A')
    idx_p = header.index('P_A')
    for line in f:
        p = line.strip().split('\t')
        if len(p) > idx_p and p[idx_beta] != 'nan':
            cosmic_data[p[idx_snp]] = (float(p[idx_beta]), float(p[idx_se]), float(p[idx_p]))

matched_b, matched_se, matched_p = [], [], []
g_b, g_se, g_p = [], [], []

for snp, (cb, cse, cp) in cosmic_data.items():
    if snp in gcta_data:
        gb, gse, gp = gcta_data[snp]
        matched_b.append(cb)
        matched_se.append(cse)
        matched_p.append(-math.log10(cp + 1e-300))
        g_b.append(gb)
        g_se.append(gse)
        g_p.append(-math.log10(gp + 1e-300))

print(f"Matched {len(matched_b)} SNPs with GCTA.")
print(f"Beta R^2:      {pearsonr(matched_b, g_b)**2:.6f}")
print(f"SE R^2:        {pearsonr(matched_se, g_se)**2:.6f}")
print(f"-log10(P) R^2: {pearsonr(matched_p, g_p)**2:.6f}")
