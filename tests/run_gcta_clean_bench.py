import os, subprocess, math

def run_cmd(cmd):
    print(f"Running: {cmd}")
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"Error: {res.stderr}")
        exit(1)

root = os.environ.get("COSMICSCANNER_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
data_dir = os.environ.get("BENCHMARK_DATA_DIR", os.path.join(root, "examples"))
res_dir = os.environ.get("BENCHMARK_RESULT_DIR", os.path.join(root, "tests/results"))
plink2 = os.environ.get("PLINK2_BIN", "plink2")
gcta64 = os.environ.get("GCTA_BIN", "gcta64")
gemma = os.environ.get("GEMMA_BIN", "gemma")
cosmic = os.environ.get("COSMICSCANNER_BIN", os.path.join(root, "build/cosmicscanner"))

# 1. Clean data for GCTA
os.chdir(f"{data_dir}/mouse_hs1940")
with open("mouse_hs1940.fam") as f:
    with open("remove.txt", "w") as fout:
        for line in f:
            p = line.strip().split()
            if "NA" in p[0:2]:
                fout.write(f"{p[0]}\t{p[1]}\n")

run_cmd(f"{plink2} --bfile mouse_hs1940 --remove remove.txt --make-bed --out {data_dir}/mouse_clean/mouse_clean")

os.chdir(f"{data_dir}/mouse_clean")
# Recreate pheno
pheno_dict = {}
with open(f"{data_dir}/mouse_hs1940/mouse_hs1940.pheno") as f:
    for line in f:
        p = line.strip().split()
        if len(p) >= 3:
            pheno_dict[f"{p[0]}_{p[1]}"] = p[2]
            
with open("mouse_clean.fam") as f, open("mouse_clean.pheno", "w") as fout:
    for line in f:
        p = line.strip().split()
        key = f"{p[0]}_{p[1]}"
        val = pheno_dict.get(key, "NA")
        fout.write(f"{p[0]}\t{p[1]}\t{val}\n")

# 2. Run GCTA
os.chdir(res_dir)
run_cmd(f"{gcta64} --bfile {data_dir}/mouse_clean --make-grm --out mouse_clean_gcta")
run_cmd(f"{gcta64} --bfile {data_dir}/mouse_clean --grm mouse_clean_gcta --pheno {data_dir}/mouse_clean.pheno --mlma --out gcta_clean_mlma")

# 3. Run CosmicScanner on cleaned data
run_cmd(f"{cosmic} --bfile {data_dir}/mouse_clean --pheno {data_dir}/mouse_clean.pheno --gwas --no-loco --out cosmic_clean_lmm")

# 4. Run GEMMA on cleaned data
run_cmd(f"{gemma} -bfile {data_dir}/mouse_clean -gk 1 -o mouse_clean")
run_cmd(f"{gemma} -bfile {data_dir}/mouse_clean -k output/mouse_clean.cXX.txt -lmm 4 -o gemma_clean_lmm")

# 5. Read & Compare
# Read GCTA
gcta_data = {}
with open("gcta_clean_mlma.mlma") as f:
    header = next(f).strip().split()
    for line in f:
        p = line.strip().split()
        if len(p) > 8 and p[6] != 'NA' and p[6] != 'nan':
            gcta_data[p[1]] = (float(p[6]), float(p[7]), float(p[8]))

# Read Cosmic
cosmic_data = {}
with open("cosmic_clean_lmm.assoc.txt") as f:
    header = next(f).strip().split('\t')
    idx_snp = header.index('SNP')
    idx_beta = header.index('BETA_A')
    idx_se = header.index('SE_A')
    idx_p = header.index('P_A')
    for line in f:
        p = line.strip().split('\t')
        if len(p) > idx_p and p[idx_beta] != 'nan':
            cosmic_data[p[idx_snp]] = (float(p[idx_beta]), float(p[idx_se]), float(p[idx_p]))

# Read GEMMA
gemma_data = {}
with open("output/gemma_clean_lmm.assoc.txt") as f:
    header = next(f).strip().split()
    idx_rs = header.index('rs')
    idx_beta = header.index('beta')
    idx_se = header.index('se')
    idx_p = header.index('p_wald')
    for line in f:
        p = line.strip().split()
        if len(p) > idx_p and p[idx_beta] != 'NA':
            gemma_data[p[idx_rs]] = (float(p[idx_beta]), float(p[idx_se]), float(p[idx_p]))

def calc_diff(d1, d2):
    max_db, max_dse, max_dp = 0, 0, 0
    for snp in d1:
        if snp in d2:
            b1, se1, p1 = d1[snp]
            b2, se2, p2 = d2[snp]
            max_db = max(max_db, abs(b1-b2))
            max_dse = max(max_dse, abs(se1-se2))
            max_dp = max(max_dp, abs(p1-p2))
    return max_db, max_dse, max_dp

print("=== CosmicScanner vs GCTA ===")
cb, cse, cp = calc_diff(cosmic_data, gcta_data)
print(f"Max Diff BETA: {cb:.2e}")
print(f"Max Diff SE:   {cse:.2e}")
print(f"Max Diff P:    {cp:.2e}")

print("=== GEMMA vs GCTA ===")
gb, gse, gp = calc_diff(gemma_data, gcta_data)
print(f"Max Diff BETA: {gb:.2e}")
print(f"Max Diff SE:   {gse:.2e}")
print(f"Max Diff P:    {gp:.2e}")

print("=== CosmicScanner vs GEMMA ===")
cgb, cgse, cgp = calc_diff(cosmic_data, gemma_data)
print(f"Max Diff BETA: {cgb:.2e}")
print(f"Max Diff SE:   {cgse:.2e}")
print(f"Max Diff P:    {cgp:.2e}")
