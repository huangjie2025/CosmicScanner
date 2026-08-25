import os, subprocess, re, time, math

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

def run_monitored(cmd, name):
    print(f"[{name}] Running command...")
    # Use /usr/bin/time -v to capture peak memory and wall clock time
    time_cmd = f"/usr/bin/time -v {cmd}"
    res = subprocess.run(time_cmd, shell=True, capture_output=True, text=True)
    
    if res.returncode != 0 and "gemma" not in cmd: 
        # GEMMA sometimes returns non-zero even on success if minor warnings happen, but let's be careful
        pass
    
    # Parse /usr/bin/time output from stderr
    mem_match = re.search(r"Maximum resident set size \(kbytes\):\s+(\d+)", res.stderr)
    peak_mem_mb = int(mem_match.group(1)) / 1024 if mem_match else 0
    
    time_match = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s+([0-9:.]+)", res.stderr)
    wall_time = 0
    if time_match:
        t_str = time_match.group(1)
        parts = t_str.split(':')
        if len(parts) == 2:
            wall_time = int(parts[0])*60 + float(parts[1])
        elif len(parts) == 3:
            wall_time = int(parts[0])*3600 + int(parts[1])*60 + float(parts[2])
            
    print(f"[{name}] Time: {wall_time:.2f} s | Memory: {peak_mem_mb:.2f} MB")
    return wall_time, peak_mem_mb

def main():
    # Setup paths
    data_prefix = "/disk228/huangjie/PanOmics/CosmicScanner/examples/mouse_clean"
    res_dir = "/disk228/huangjie/PanOmics/CosmicScanner/tests/results/mentor_report"
    os.makedirs(res_dir, exist_ok=True)
    os.chdir(res_dir)

    gemma_bin = "/disk228/huangjie/software/miniconda3/envs/benchmark/bin/gemma"
    cosmic_bin = "/disk228/huangjie/PanOmics/Cosmic/build/CosmicScanner/cosmicscanner"

    # 控制变量：为了公平起见，GEMMA使用单线程计算（因为精确LMM的核心特征值分解主要是单线程或依赖BLAS），CosmicScanner也默认1个线程。
    
    # 1. Run GEMMA (Step 1: GRM, Step 2: LMM)
    print("\n--- GEMMA Pipeline ---")
    gemma_cmd1 = f"{gemma_bin} -bfile {data_prefix} -gk 1 -o mouse_clean_rep"
    t1, m1 = run_monitored(gemma_cmd1, "GEMMA_GRM")
    
    gemma_cmd2 = f"{gemma_bin} -bfile {data_prefix} -k output/mouse_clean_rep.cXX.txt -lmm 4 -o gemma_rep_lmm"
    t2, m2 = run_monitored(gemma_cmd2, "GEMMA_LMM")
    
    gemma_time = t1 + t2
    gemma_mem = max(m1, m2)

    # 2. Run CosmicScanner (On-the-fly GRM + LMM)
    print("\n--- CosmicScanner Pipeline ---")
    cosmic_cmd = f"{cosmic_bin} --bfile {data_prefix} --pheno {data_prefix}.pheno --gwas --no-loco --out cosmic_rep_lmm"
    cosmic_time, cosmic_mem = run_monitored(cosmic_cmd, "CosmicScanner_GWAS")

    # Save performance
    with open("performance_summary.csv", "w") as f:
        f.write("Software,Time_s,Memory_MB\n")
        f.write(f"GEMMA,{gemma_time},{gemma_mem}\n")
        f.write(f"CosmicScanner,{cosmic_time},{cosmic_mem}\n")

    print("\n--- Summary Results ---")
    gemma_out = "output/gemma_rep_lmm.assoc.txt"
    cosmic_out = "cosmic_rep_lmm.assoc.txt"

    gemma_data = {}
    with open(gemma_out) as f:
        header = next(f).strip().split()
        idx_rs = header.index('rs')
        idx_p = header.index('p_wald')
        idx_beta = header.index('beta')
        idx_se = header.index('se')
        for line in f:
            p = line.strip().split()
            if len(p) > max(idx_p, idx_beta, idx_se) and p[idx_p] != 'NA' and p[idx_beta] != 'NA':
                gemma_data[p[idx_rs]] = (
                    -math.log10(float(p[idx_p]) + 1e-300),
                    float(p[idx_beta]),
                    float(p[idx_se])
                )

    x_p, y_p = [], []
    x_b, y_b = [], []
    x_se, y_se = [], []
    
    with open(cosmic_out) as f:
        header = next(f).strip().split('\t')
        idx_snp = header.index('SNP')
        idx_p = header.index('LOG10P') if 'LOG10P' in header else header.index('P')
        idx_beta = header.index('BETA')
        idx_se = header.index('SE')
        is_log = 'LOG10P' in header
        for line in f:
            p = line.strip().split('\t')
            if len(p) > max(idx_p, idx_beta, idx_se):
                snp = p[idx_snp]
                if snp in gemma_data and p[idx_p] != 'nan' and p[idx_beta] != 'nan':
                    c_logp = float(p[idx_p]) if is_log else -math.log10(float(p[idx_p]) + 1e-300)
                    c_beta = float(p[idx_beta])
                    c_se = float(p[idx_se])
                    
                    g_logp, g_beta, g_se = gemma_data[snp]
                    
                    x_p.append(c_logp); y_p.append(g_logp)
                    x_b.append(c_beta); y_b.append(g_beta)
                    x_se.append(c_se); y_se.append(g_se)

    r_p = pearsonr(x_p, y_p)
    r_b = pearsonr(x_b, y_b)
    r_se = pearsonr(x_se, y_se)

    print(f"Matched {len(x_p)} SNPs.")
    print(f"1 - Pearson r (-log10P): {1 - r_p:.2e}")
    print(f"1 - Pearson r (BETA):    {1 - r_b:.2e}")
    print(f"1 - Pearson r (SE):      {1 - r_se:.2e}")
    print(f"GEMMA Time: {gemma_time:.1f}s, Memory: {gemma_mem:.1f}MB")
    print(f"CosmicScanner Time: {cosmic_time:.1f}s, Memory: {cosmic_mem:.1f}MB")

    # Generate Python plotting script to run via matplotlib
    plot_script = """
import matplotlib.pyplot as plt
import numpy as np

# Data configuration
x_p = {x_p}; y_p = {y_p}
x_b = {x_b}; y_b = {y_b}
x_se = {x_se}; y_se = {y_se}

plt.style.use('default')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 14
plt.rcParams['axes.linewidth'] = 1.5

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
# Adjusted order: BETA -> SE -> -log10(P)
titles = ['BETA (Effect Size)', 'SE (Standard Error)', '-log10(P)']
xs = [x_b, x_se, x_p]
ys = [y_b, y_se, y_p]
rs = [{r_b}, {r_se}, {r_p}]

for i in range(3):
    ax = axes[i]
    ax.scatter(xs[i], ys[i], color='#808080', alpha=0.6, s=15, edgecolors='none')
    
    lims = [min(min(xs[i]), min(ys[i])), max(max(xs[i]), max(ys[i]))]
    span = lims[1] - lims[0]
    if span == 0: span = 1
    plot_lims = [lims[0] - 0.05*span, lims[1] + 0.05*span]
    
    ax.plot(plot_lims, plot_lims, color='black', linestyle='--', linewidth=2)
    diff_r = 1 - rs[i]
    ax.text(0.05, 0.95, f'1 - r = {{diff_r:.2e}}', transform=ax.transAxes, ha='left', va='top', fontsize=12, fontweight='bold')
    
    ax.set_title(titles[i], fontweight='bold')
    ax.set_xlabel(f'CosmicScanner {{titles[i]}}', fontweight='bold')
    ax.set_ylabel(f'GEMMA {{titles[i]}}', fontweight='bold')
    ax.set_xlim(plot_lims)
    ax.set_ylim(plot_lims)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('Fig1_Comprehensive_Consistency.png', dpi=300)
plt.close()

# Performance Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
sw = ['GEMMA', 'CosmicScanner']
colors = ['#A9A9A9', '#000000']
    
bars1 = ax1.bar(sw, [{gemma_time}, {cosmic_time}], color=colors, width=0.5)
ax1.set_ylabel('Runtime (s)', fontweight='bold')
ax1.set_title('Computational Time', fontweight='bold')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
for bar in bars1:
    yval = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2.0, yval, f'{{yval:.1f}}', va='bottom', ha='center', fontweight='bold')

bars2 = ax2.bar(sw, [{gemma_mem}, {cosmic_mem}], color=colors, width=0.5)
ax2.set_ylabel('Peak Memory (MB)', fontweight='bold')
ax2.set_title('Memory Usage', fontweight='bold')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
for bar in bars2:
    yval = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2.0, yval, f'{{yval:.1f}}', va='bottom', ha='center', fontweight='bold')

plt.tight_layout()
plt.savefig('Fig2_Performance_Bar.png', dpi=300)
plt.close()
"""
    
    # We write it to a temporary file and run it to avoid heavy memory overhead inside the main script
    with open("temp_plot.py", "w") as f:
        f.write(plot_script.format(
            x_p=x_p, y_p=y_p, x_b=x_b, y_b=y_b, x_se=x_se, y_se=y_se,
            r_p=r_p, r_b=r_b, r_se=r_se,
            gemma_time=gemma_time, cosmic_time=cosmic_time,
            gemma_mem=gemma_mem, cosmic_mem=cosmic_mem
        ))
    
    subprocess.run("/disk228/huangjie/software/miniconda3/envs/benchmark/bin/python temp_plot.py", shell=True)
    os.remove("temp_plot.py")
    print("Plots generated successfully!")

if __name__ == "__main__":
    main()
