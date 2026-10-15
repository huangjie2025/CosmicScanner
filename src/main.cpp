#ifdef _WIN32
#include <windows.h>
#endif
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <regex>
#include <cmath>
#include <stdexcept>
#include <vector>
#include "config.h"
#include "gwas_app.h"
#include "solver_planner.h"
#include "cli_frontend.h"

using namespace std;
using namespace cosmic;

static string json_escape(const string& s) {
    string out;
    for (char c : s) {
        if (c == '"' || c == '\\') { out.push_back('\\'); out.push_back(c); }
        else if (c == '\n') out += "\\n";
        else out.push_back(c);
    }
    return out;
}

static string read_text_file(const string& path) {
    ifstream in(path);
    if (!in) throw runtime_error("INPUT_CONFIG_NOT_READABLE: Cannot open input config: " + path);
    stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

static string config_string(const string& text, const string& key, const string& current = "") {
    regex re("\\\"" + key + "\\\"\\s*:\\s*\\\"([^\\\"]*)\\\"");
    smatch m;
    return regex_search(text, m, re) ? m[1].str() : current;
}

static bool config_bool(const string& text, const string& key, bool current = false) {
    regex re("\\\"" + key + "\\\"\\s*:\\s*(true|false)");
    smatch m;
    return regex_search(text, m, re) ? (m[1].str() == "true") : current;
}

static string config_array_as_csv(const string& text, const string& key, const string& current = "") {
    regex re("\\\"" + key + "\\\"\\s*:\\s*\\[([^\\]]*)\\]");
    smatch m;
    if (!regex_search(text, m, re)) return current;
    string body = m[1].str();
    regex item_re("\\\"([^\\\"]*)\\\"");
    string out;
    for (sregex_iterator it(body.begin(), body.end(), item_re), end; it != end; ++it) {
        if (!out.empty()) out += ",";
        out += (*it)[1].str();
    }
    return out;
}

static void apply_input_config(GwasConfig& cfg, const string& path) {
    string text = read_text_file(path);
    cfg.input_config = path;
    cfg.input_file = config_string(text, "input_file", cfg.input_file);
    cfg.input_format = config_string(text, "input_format", cfg.input_format);
    cfg.input_model = config_string(text, "input_model", cfg.input_model);
    cfg.prepared_input_dir = config_string(text, "prepared_input_dir", cfg.prepared_input_dir);
    cfg.id_col = config_string(text, "id_col", cfg.id_col);
    cfg.iid_col = config_string(text, "iid_col", cfg.iid_col);
    cfg.time_col = config_string(text, "time_col", cfg.time_col);
    cfg.input_trait = config_string(text, "input_trait", cfg.input_trait);
    cfg.trait_cols = config_array_as_csv(text, "trait_cols", cfg.trait_cols);
    cfg.record_qcovars = config_array_as_csv(text, "record_qcovars", cfg.record_qcovars);
    cfg.record_ccovars = config_array_as_csv(text, "record_ccovars", cfg.record_ccovars);
    cfg.subject_qcovars = config_array_as_csv(text, "subject_qcovars", cfg.subject_qcovars);
    cfg.subject_ccovars = config_array_as_csv(text, "subject_ccovars", cfg.subject_ccovars);
    cfg.time_parser = config_string(text, "time_parser", cfg.time_parser);
    cfg.time_unique_policy = config_string(text, "time_unique_policy", cfg.time_unique_policy);
    cfg.time_index_multiplier = config_string(text, "time_index_multiplier", cfg.time_index_multiplier);
    cfg.preview_schema = config_bool(text, "preview_schema", cfg.preview_schema);
    cfg.infer_columns = config_bool(text, "infer_columns", cfg.infer_columns);
    cfg.wide_to_long = config_bool(text, "wide_to_long", cfg.wide_to_long);
}

static vector<string> csv_to_vec(const string& s) {
    vector<string> out;
    stringstream ss(s);
    string tok;
    while (getline(ss, tok, ',')) if (!tok.empty()) out.push_back(tok);
    return out;
}

static void write_json_array(ofstream& out, const string& key, const vector<string>& vals, bool comma = true) {
    out << "  \"" << key << "\": [";
    for (size_t i = 0; i < vals.size(); ++i) out << (i ? ", " : "") << "\"" << json_escape(vals[i]) << "\"";
    out << "]" << (comma ? "," : "") << "\n";
}

static void write_template_config(const GwasConfig& cfg, const string& path) {
    ofstream out(path);
    if (!out) throw runtime_error("TEMPLATE_CONFIG_NOT_WRITABLE: Cannot write template config: " + path);
    out << "{\n";
    out << "  \"input_file\": \"" << json_escape(cfg.input_file.empty() ? "data/semen.csv" : cfg.input_file) << "\",\n";
    out << "  \"input_format\": \"" << json_escape(cfg.input_format) << "\",\n";
    out << "  \"input_model\": \"" << json_escape(cfg.input_model) << "\",\n";
    out << "  \"prepared_input_dir\": \"" << json_escape(cfg.prepared_input_dir.empty() ? "prep/semen" : cfg.prepared_input_dir) << "\",\n";
    out << "  \"id_col\": \"" << json_escape(cfg.id_col.empty() ? "pig_id" : cfg.id_col) << "\",\n";
    out << "  \"iid_col\": \"" << json_escape(cfg.iid_col) << "\",\n";
    out << "  \"time_col\": \"" << json_escape(cfg.time_col.empty() ? "collect_date" : cfg.time_col) << "\",\n";
    out << "  \"input_trait\": \"" << json_escape(cfg.input_trait) << "\",\n";
    write_json_array(out, "trait_cols", csv_to_vec(cfg.trait_cols.empty() ? "volume,concentration,motility" : cfg.trait_cols));
    write_json_array(out, "record_qcovars", csv_to_vec(cfg.record_qcovars.empty() ? "age,temp" : cfg.record_qcovars));
    write_json_array(out, "record_ccovars", csv_to_vec(cfg.record_ccovars));
    write_json_array(out, "subject_qcovars", csv_to_vec(cfg.subject_qcovars));
    write_json_array(out, "subject_ccovars", csv_to_vec(cfg.subject_ccovars.empty() ? "farm,line" : cfg.subject_ccovars));
    out << "  \"time_parser\": \"" << json_escape(cfg.time_parser) << "\",\n";
    out << "  \"time_unique_policy\": \"" << json_escape(cfg.time_unique_policy) << "\",\n";
    out << "  \"time_index_multiplier\": \"" << json_escape(cfg.time_index_multiplier) << "\",\n";
    out << "  \"wide_to_long\": " << (cfg.wide_to_long ? "true" : "false") << ",\n";
    out << "  \"infer_columns\": " << (cfg.infer_columns ? "true" : "false") << "\n";
    out << "}\n";
}

void print_usage(const char* prog_name) {
    cout << "CosmicScanner v" << COSMICSCANNER_VERSION << "\n";
    cout << "  --version             Print software version and exit\n";
    cout << "Public scope: A=scoped official contract, B=experimental, C=internal/legacy.\n";
    cout << "See COSMICSCANNER_V2_SCOPE_FREEZE_PROPOSAL.md; final acceptance is pending.\n";
    cout << "Formats, engines and output schemas are model-specific; no universal equivalence.\n";
    cout << "Usage: " << prog_name << " --bfile <prefix> [options]\n";
    cout << "\n=== INPUT/OUTPUT ===\n";
    cout << "  --pca --pca-k <n>     Experimental BED GRM PCA; default k=10\n";
    cout << "  --geno <f> --mind <f> Maximum variant/sample missing fraction [0,1]; BED QC\n";
    cout << "  --keep <f> --remove <f> BED sample lists: FID IID; remove takes precedence\n";
    cout << "  --ld-prune            Experimental BED window LD pruning\n";
    cout << "  --ld-clump <f>        Experimental clumping; summary header SNP P CHR BP\n";
    cout << "  --ld-window <bp> --ld-step <n> --ld-r2 <f> Defaults: 250000, 1, 0.2\n";
    cout << "  --clump-p <f>         Lead P threshold (default 5e-8)\n";
    cout << "  --skat-o              Fixed-grid dependence-adjusted SKAT-O\n";
    cout << "  --acat-v              Experimental set-wise Cauchy combination (--set-file)\n";
    cout << "  --rare-max-maf <f>    Maximum set-test MAF in (0,0.5] (default 0.5)\n";
    cout << "  --rare-min-mac <n>    Minimum set-test minor allele count (default 1)\n";
    cout << "  --snp-model <m>       additive (default), dominant, recessive, heterozygote\n";
    cout << "                        Fixed-effect marker coding for standard LM/LMM;\n";
    cout << "                        non-additive --snp-model requires hard-call genotypes.\n";
    cout << "                        --gwas-model dominance/additive_dominance requires\n";
    cout << "                        additive --snp-model.\n";
    cout << "  --gwas-model <m>      additive (default) | additive_dominance | dominance\n";
    cout << "                        additive:          a + A\n";
    cout << "                        additive_dominance: a + d + A + D LMM scan (PLINK BED or PGEN hardcall)\n";
    cout << "                        dominance:         d + D LMM scan (BED/PGEN hardcall or BGEN probability)\n";
    cout << "  --bgen-dom-max-missing <f>  BGEN dominance variant missingness filter (default 0.05)\n";
    cout << "  --bgen-dom-min-info <f>     Minimum dominance INFO_D (default 0.01)\n";
    cout << "  --bgen-dom-entropy-warning <f>  Normalized entropy warning threshold (default 0.50)\n";
    cout << "  --farmcpu             Experimental native FarmCPU multi-locus GWAS\n";
    cout << "  --farmcpu-max-iterations <n>  Iteration cap (default 5)\n";
    cout << "  --farmcpu-pseudo-qtn <n>      Pseudo-QTN cap (default 10)\n";
    cout << "  --blink               Experimental native BLINK multi-locus GWAS\n";
    cout << "  --condition-snp <f>   Experimental conditional GWAS; one lead SNP ID per line\n";
    cout << "  --conditional-joint  Experimental LD/BIC stepwise joint selection for --condition-snp\n";
    cout << "  --firth-logistic     Experimental Firth penalized binary GWAS\n";
    cout << "  --logistic <test>    Ordinary logistic GWAS: score or wald\n";
    cout << "  --logistic-glmm      B: PQL mixed null with full projected working-score variance\n";
    cout << "                        Not an exact mixed-likelihood fit.\n";
    cout << "  --saige-style        Experimental binary mixed score/SPA scan; requires --biobank-mode\n";
    cout << "  --cojo               COJO-style conditional analysis from summary stats + LD reference\n";
    cout << "  --cojo-slct          COJO stepwise conditional analysis (forward selection)\n";
    cout << "  --cojo-joint         COJO joint model fit of selected SNPs\n";
    cout << "  --cojo-summary <f>   GWAS summary statistics file\n";
    cout << "  --cojo-ld-ref <f>    Dense signed-r LD: SNP CHR BP A1 A2 <ordered SNP IDs>\n";
    cout << "  --cojo-p <f>         P-value threshold for COJO selection (default: 5e-8)\n";
    cout << "  --cojo-r2 <f>        LD r^2 threshold for collinearity (default: 0.9)\n";
    cout << "  --cojo-wind <bp>     Window size in bp for LD (default: 10000000)\n";
    cout << "  --cojo-cond <f>      Condition on these SNP IDs (one per line)\n";
    cout << "  --fine-map           Native summary-statistic fine-mapping with real genotype LD\n";
    cout << "  --finemap-summary <f> Canonical .cosmic.sumstats input\n";
    cout << "  --finemap-region <r> Locus in chr:start-end form\n";
    cout << "  --credible-level <f> Credible-set posterior level (default: 0.95)\n";
    cout << "  --finemap-prior-var <f> Prior effect variance W (default: 0.04)\n";
    cout << "  --finemap-max-signals <n> Maximum conditional signals (default: 5)\n";
    cout << "  --finemap-signal-p <f> Conditional signal threshold (default: 5e-8)\n";
    cout << "  --meta-analysis      Summary-statistics meta-analysis (fixed/random effects)\n";
    cout << "  --meta-studies <f>   File listing study summary stat files\n";
    cout << "  --meta-method <m>    Meta-analysis method: fixed, random-dl, random-reml (default: fixed)\n";
    cout << "  --meta-het-p <f>     Heterogeneity threshold p-value (default: 0.05)\n";
    cout << "  --meta-min-studies <n> Min studies required (default: 2)\n";
    cout << "  --meta-gc            Apply genomic control correction\n";
    cout << "  --meta-gc-lambda <f> Genomic control lambda (default: 1.0)\n";
    cout << "  --meta-out-format <f> Output format: metal, simple (default: metal)\n";
    cout << "  --bfile <prefix>       PLINK 1.x binary prefix (.bed/.bim/.fam)\n";
    cout << "  --bgen <file>          BGEN file (Layout 2, zlib compressed)\n";
    cout << "  --sample <file>        Sample file for BGEN (Oxford format)\n";
    cout << "  --pfile <prefix>       PLINK 2.0 prefix (.pgen/.pvar/.psam)\n";
    cout << "  --pgen <file>          PLINK 2.0 genotype file\n";
    cout << "  --pvar <file>          PLINK 2.0 variant file\n";
    cout << "  --psam <file>          PLINK 2.0 sample file\n";
    cout << "  --vcf <file>           VCF file (bgzip compressed, tabix indexed)\n";
    cout << "  --vcf-tbi <file>       Tabix index file for VCF (optional if .tbi alongside)\n";
    cout << "  --vcf-sample <file>    Sample mapping file for VCF (optional)\n";
    cout << "  --hmp <file>           Standard HapMap file (11 metadata columns + samples)\n";
    cout << "  --geno-matrix <file>   Text genotype matrix (default: variant rows, sample columns)\n";
    cout << "  --geno-var <file>      Text variant metadata: ID CHROM POS REF ALT\n";
    cout << "  --geno-sample <file>   Optional text sample metadata: IID or FID IID\n";
    cout << "  --geno-orientation <o> variants-by-samples (default) or samples-by-variants\n";
    cout << "  --npy <file>           NPY genotype matrix (float32, little-endian)\n";
    cout << "  --npy-var <file>       NPY variant metadata: ID CHROM POS REF ALT\n";
    cout << "  --npy-sample <file>    NPY sample metadata: FID IID\n";
    cout << "  --npy-orientation <o>  variants-by-samples (default) or samples-by-variants\n";
    cout << "  --ped <file>           PLINK PED file\n";
    cout << "  --map <file>           PLINK MAP file (must be used with --ped)\n";
    cout << "  --convert <format>     Convert genotype format (input: --bfile/--pfile/--bgen/--vcf/--hmp/--geno-matrix/--npy/--ped+--map)\n";
    cout << "  --to <format>          Target format for conversion: text, npy, ped, vcf (requires --convert)\n";
    cout << "  --out <prefix>         Output prefix for conversion\n";
    cout << "  --batch-models <list>  C: developer batch runner: lm,lmm,fastlmm,saige,farmcpu,blink\n";
    cout << "  --batch-out <prefix>   Output prefix for batch run (each model writes its own output)\n";
    cout << "  --out/-o <prefix>      Output prefix (default: cosmic_gwas)\n";
    cout << "  --out-ssf              EBI GWAS-SSF output on supporting paths only\n";
    cout << "  --out-gz               Gzip on supporting paths; AD/dominance guards reject it\n";
    cout << "                          No all-model output schema; sumstats INFO=1 is not measured INFO.\n";
    cout << "  --extract <file>       File listing SNP IDs to include in analysis\n";
    cout << "  --exclude <file>       File listing SNP IDs to exclude from analysis\n";
    cout << "\n=== PHENOTYPE / COVARIATES ===\n";
    cout << "  --pheno/-p <file>      Phenotype file: FID IID [TIME] VALUE\n";
    cout << "  --mpheno/-n <n>        Phenotype column number, 1-indexed (default: 1)\n";
    cout << "  --tpheno <n>           Time column number for repeatability/RRM\n";
    cout << "  --qcovar <file>        Quantitative covariates: FID IID VAL1 VAL2 ...\n";
    cout << "  --covar <file>         Categorical covariates: FID IID CAT1 CAT2 ...\n";
    cout << "\n=== INPUT ADAPTATION MVP ===\n";
    cout << "  --input-file <file>    Raw csv/tsv/txt table to prepare before modeling\n";
    cout << "  --input-config <json>  Load input adaptation settings from JSON config\n";
    cout << "  --write-template-config <json> Write an editable input config and exit\n";
    cout << "  --input-format <fmt>   auto|csv|tsv|txt (default: auto; xlsx not in MVP)\n";
    cout << "  --input-model <model>  repeatability|rrm|lmm|multivariate (default: repeatability)\n";
    cout << "  --prepared-input-dir <dir> Directory for prepared manifests and model files\n";
    cout << "  --prepare-input-only   Write prepared inputs and exit before GWAS\n";
    cout << "  --id-col <name>        Raw ID column; matched to genotype IID or FID:IID\n";
    cout << "  --iid-col <name>       Optional IID column when --id-col is FID\n";
    cout << "  --time-col <name>      Record time/date/batch column\n";
    cout << "  --trait-cols <list>    Comma-separated trait columns\n";
    cout << "  --input-trait <name>   Trait to run after preparation (default: first trait)\n";
    cout << "  --record-qcovars <l>   Record-level quantitative covariates\n";
    cout << "  --record-ccovars <l>   Record-level categorical covariates\n";
    cout << "  --subject-qcovars <l>  Subject-level quantitative covariates\n";
    cout << "  --subject-ccovars <l>  Subject-level categorical covariates\n";
    cout << "  --time-parser <mode>   auto|date|integer|batch (default: auto)\n";
    cout << "  --time-unique-policy <p> strict|within-time-index (default)\n";
    cout << "  --time-index-multiplier <n|auto> Multiplier for duplicate-time keys\n";
    cout << "  --preview-schema       Write schema_preview.tsv and exit before preparation\n";
    cout << "  --infer-columns        Suggest ID/time/trait/covariate columns in preview\n";
    cout << "  --wide-to-long         Treat --trait-cols as wide trait-time columns\n";
    cout << "  --maf <float>          Minor allele frequency filter (default: 0.01)\n";
    cout << "  --hwe <float>          HWE exact test p-value filter (default: 0.0, disabled)\n";
    cout << "  --gwas                 Run standard GWAS; some model flags enable scans implicitly\n";
    cout << "  --model <name>         Model family: lmm, lm, multivariate, repeatability, rrm,\n";
    cout << "                          animal, maternal, ssgwas, common-env, threshold, count, ordinal, gxe\n";
    cout << "                          Legacy --model lmm-eb/eb-lmm remains an EB modifier.\n";
    cout << "                          Documented analysis commands should include it explicitly.\n";
    cout << "  --no-loco              Disable leave-one-chromosome-out (run global)\n";
    cout << "  --scan-tier <tier>     Per-SNP VC handling: approx, semi, exact\n";
    cout << "                        Multivariate semi requires --mv-compatible emmvgwas; native uses approx\n";
    cout << "                          approx=fixed null; semi=fixed ratio/profile scale;\n";
    cout << "                          exact=full per-SNP variance-component refit\n";
    cout << "                          Not finite-sample exact P. mv/repeatability exact tier needs no-LOCO;\n";
    cout << "                          repeatability exact tier also requires BED and no record weights.\n";
    cout << "  --fast-lmm             Legacy alias for fixed-null single-trait scan\n";
    cout << "  --fastlmm              B: BED Gaussian trace-ranked randomized low-rank LMM\n";
    cout << "  --fastlmm-trace-threshold <f> Minimum retained GRM trace fraction (default: 0.95)\n";
    cout << "  --fastlmm-min-rank <n> Minimum automatic rank (default: 1)\n";
    cout << "  --fastlmm-max-rank <n> Maximum automatic rank (default: 256)\n";
    cout << "  --fastlmm-matrix-free  B: BED matrix-free low-rank path; rank adequacy is data-dependent\n";
    cout << "  --fastlmm-no-fallback  Don't fall back to full rank if trace threshold not met\n";
    cout << "  --lmm-low-rank <k>     Opt-in randomized rank-k LMM (0 keeps exact Eigen default)\n";
    cout << "  --lmm-low-rank-oversample <n> Randomized eigensolver oversampling (default: 10)\n";
    cout << "  --lmm-low-rank-seed <n> Deterministic randomized eigensolver seed (default: 42)\n";
    cout << "  --lmm-low-rank-error-tol <f> Relative Frobenius error limit; exact fallback above it\n";
    cout << "  --lmm-gemma-mode/-lmm <n> GEMMA LMM mode: 1 Wald, 2 LRT, 3 Score, 4 all\n";
    cout << "  --lmm-gemma-refit      Per-SNP alternative REML refit for GEMMA LMM Wald/LRT\n";
    cout << "  --lm                   Ordinary linear model flag; use --lm-mode <n>\n";
    cout << "  -lm <n>                GEMMA-style short LM alias\n";
    cout << "  --lm-mode <n>          GEMMA LM mode: 1 Wald, 2 LRT, 3 Score, 4 all\n";
    cout << "  --model lmm-eb         Empirical Bayes random-SNP LMM GWAS (experimental; report standard LMM too)\n";
    cout << "  --gwas-method eb-lmm   Alias for --model lmm-eb\n";
    cout << "  --test eb-wald         Alias enabling EB-LMM Wald test\n";
    cout << "  --adjust eb-mbonf      Experimental threshold in .eb.summary.txt; do not report alone\n";
    cout << "  --spa                  Use Saddle Point Approximation for binary traits\n";
    cout << "  --threads <n>          Number of OpenMP threads (default: max)\n";
    cout << "  --memory-gib <value>   SolverPlanner diagnostic memory budget in GiB\n";
    cout << "  --chr <chr>            Restrict scan to one chromosome\n";
    cout << "  --from-bp <pos>        Restrict scan start position\n";
    cout << "  --to-bp <pos>          Restrict scan end position\n";
    cout << "\n=== VCE OPTIONS (model-dependent) ===\n";
    cout << "  --lmm-null-solver <s>  Null model solver: brent, idul, idul_plus, aidul\n";
    cout << "                          (default: brent)\n";
    cout << "  --reml-maxit <n>       Max VCE iterations (default: 50)\n";
    cout << "  --reml-tol <f>         VCE convergence tolerance (default: 1e-6)\n";
    cout << "  --reml-vg <f>          Fix genetic variance (skip VCE)\n";
    cout << "  --reml-vpe <f>         Fix permanent environment variance (skip VCE)\n";
    cout << "  --reml-ve <f>          Fix residual variance (skip VCE)\n";
    cout << "\n=== MODEL SELECTION (choose exactly one) ===\n";
    cout << "\n--- Standard LMM (default, no flag needed) ---\n";
    cout << "  Exact single-trait mixed model GWAS with Eigen decomposition.\n";
    cout << "  Standard for quantitative traits with one record per individual.\n";
    cout << "\n--- Repeatability Model ---\n";
    cout << "  --repeatability-model  Multi-record per individual (long format).\n";
    cout << "  --pheno format: FID IID TIME VALUE  (--tpheno <col>)\n";
    cout << "  VCE: additive genetic + permanent environment + residual.\n";
    cout << "  --repeatability-mode <m> Strategy: auto, exact, fast, stcg\n";
    cout << "                          A: BED record-level approx/semi; B: fast/STCG/full refit.\n";
    cout << "                          Auto may switch strategy; exact mode differs from exact tier.\n";
    cout << "  --repeatability-grm-cache <p> Save/load subject GRM\n";
    cout << "  --repeatability-weight <file> Record weights: FID IID WEIGHT or FID IID TIME WEIGHT\n";
    cout << "  --vce-mode <m>         Repeatability VCE refinement: auto, fast,\n";
    cout << "                          subject-level-exact, full-record-exact\n";
    cout << "                          auto uses scalable subject-level Woodbury VCE;\n";
    cout << "                          full-record-exact is validation-only.\n";
    cout << "  --repeatability-refit-validation Per-SNP VCE refit (small data)\n";
    cout << "  --repeatability-refit-limit <n> Refit limit (0 = all)\n";
    cout << "\n--- Common-Environment Model (Experimental) ---\n";
    cout << "  --common-env-model     y = Xb + x_j*beta_j + Z_a*a + Z_c*c + e\n";
    cout << "                         a ~ N(0, K*sigma2_a) (genomic GRM), c ~ N(0, I*sigma2_c)\n";
    cout << "  --pheno format: FID IID [TIME] VALUE  (long format, --tpheno <col>)\n";
    cout << "  --common-env-file <f>  FID IID GROUP per subject, or per record with\n";
    cout << "  --common-env-record-level\n";
    cout << "  --reml-vg/-vc/-ve <v>  Fix variance components (skip VCE)\n";
    cout << "\n--- Random Regression Model (RRM, B experimental) ---\n";
    cout << "  --rrm-order <n>        Basis order/df-1 (Legendre default, 1-4)\n";
    cout << "  --rrm-basis <m>        Basis: legendre (default) or bspline\n";
    cout << "  --rrm-curve-grid <n>   Time grid points for dynamic curves (default: 21)\n";
    cout << "  --rrm-top <n>          Number of top SNPs for RRM curve outputs (default: 100)\n";
    cout << "  --rrm-ai-iter <n>      RRM null AI refinement iterations (-1=auto, 0=fast)\n";
    cout << "  --rrm-sig-threshold <p> P threshold for dynamic time windows (default: 0.05)\n";
    cout << "  --rrm-unbalance         Unbalanced longitudinal GWAS (long-format --pheno)\n";
    cout << "  --rrm-unbalance        Unbalanced longitudinal RRM scan (default: coefficient method)\n";
    cout << "  --rrm-gwas-method <m>  Unbalanced RRM-GWAS: record (record-level x*Phi Method A) | coefficient (two-stage Method B, default)\n";
    cout << "  --rrm-loco-global-null  Fast approximate RRM LOCO using global null model\n";
    cout << "  --tpoint <time_col>     Column name for time points (unbalanced mode)\n";
    cout << "  --tclass <cols>         Time-dependent class covariates (comma-sep column names)\n";
    cout << "  --tcovar <cols>         Time-dependent q-covariates (comma-sep column names)\n";
    cout << "  --rrm-compatible <m>   Compatibility mode: gma-fixreg (balanced PLINK) | gmat-unbalance\n";
    cout << "  --rrm-compatible-varcom <file> GMA .var file for gma-fixreg mode\n";
    cout << "  --long-pheno <file>    Long-format pheno: FID IID TIME VALUE\n";
    cout << "  Fits subject-level Legendre coefficients with GRM-based\n";
    cout << "  random regression. Outputs joint multi-df chi-square test.\n";
    cout << "\n--- Multivariate Model ---\n";
    cout << "  --multivariate         Requires --mpheno-list <c1,c2,...>\n";
    cout << "                          A: scoped BED reference; PGEN/BGEN extensions B; VCF unsupported.\n";
    cout << "                          Hybrid requires BED; exact scan tier requires no-LOCO.\n";
    cout << "  --mpheno-list <list>   Comma-separated phenotype column numbers\n";
    cout << "  --mv-method <m>        emmv, wald, eb, gemma, gemma-full, gemma-exact, gemma-gxe, gemma-gxe-full, hybrid\n";
    cout << "  --mv-test <list>       GEMMA-full tests: wald,lrt,score\n";
    cout << "  --mv-gemma-mode <n>    GEMMA mvLMM mode: 1 Wald, 2 LRT, 3 Score, 4 all\n";
    cout << "  --mv-pnr <p>           GEMMA mvLMM refinement threshold (default 0.001)\n";
    cout << "  --mv-emi <n>           GEMMA exact EM iterations (default 10000)\n";
    cout << "  --mv-nri <n>           GEMMA exact refinement iterations (default 100)\n";
    cout << "  --mv-emp <f>           GEMMA exact EM precision (default 0.0001)\n";
    cout << "  --mv-nrp <f>           GEMMA exact refinement precision (default 0.0001)\n";
    cout << "  --mv-crt              GEMMA-style Edgeworth/PCRT p-value correction after refinement\n";
    cout << "  --mv-refine-discrepancy-top <n> Hybrid: add top Score/Wald discrepancy SNPs\n";
    cout << "  --mv-refine-discrepancy <x>     Hybrid: min |log10P score-wald diff| (1.0)\n";
    cout << "  --mv-refine-score-prepass <n>   Hybrid: Score-only audit top n Stage1 SNPs\n";
    cout << "  --mv-vce-solver <s>    em (default) or aireml\n";
    cout << "  --mv-em-min-iter <n>   Min EM iterations before AI-REML (10)\n";
    cout << "  --mv-ai-max-iter <n>   Max AI-REML iterations (50)\n";
    cout << "  --mv-output-vbeta      Add upper-triangle beta covariance columns\n";
    cout << "  --mv-emmvgwas-output   Also write EMmvGWAS-style G/R/beta/association CSV files\n";
    cout << "  --mv-emmvgwas-G <csv>  EMmvGWAS semi_gwas mode: external G.csv from vc_estimation\n";
    cout << "  --mv-emmvgwas-R <csv>  EMmvGWAS semi_gwas mode: external R.csv from vc_estimation\n";
    cout << "  --emmvgwas             Alias: --multivariate --mv-method emmv --mv-compatible emmvgwas --no-loco\n";
    cout << "  --emmvgwas-vc-only     EMmvGWAS vc_estimation route: write G/R/beta and stop before scan\n";
    cout << "  --emmvgwas-output      Alias for --mv-emmvgwas-output\n";
    cout << "  --emmvgwas-G/R <csv>   Aliases for EMmvGWAS semi_gwas external G/R inputs\n";
    cout << "  --env-covar <col>      Environment column for --mv-method gemma-gxe or --gxe-model\n";
    cout << "\n--- Maternal Effects Model (B experimental) ---\n";
    cout << "  --maternal-model       Direct + maternal + PE multi-component model\n";
    cout << "  --pheno format: FID IID DAM_ID PHENO  (4+ columns)\n";
    cout << "  --maternal-id <col>    Maternal ID column (default: 3, 1-indexed)\n";
    cout << "  VCE: direct additive (ad), maternal additive (am),\n";
    cout << "       direct-maternal covariance (adm), PE (c), residual (e)\n";
    cout << "  Scan uses effective VCs: Vg=ad+am+(c) with LOCO GRM.\n";
    cout << "\n--- Single-Step GWAS (ssGWAS) Model ---\n";
    cout << "  --single-step          B Method 1 by default: genotyped H-subset fixed-SNP LMM\n";
    cout << "                          Not full-pedigree MME association; no Method-1 EBV writer.\n";
    cout << "  --make-agh             Build A-inverse, G, and H-inverse matrices and exit\n";
    cout << "  --write-agh            Write A-inverse, G, and H-inverse during ssGWAS\n";
    cout << "  --ped <file>           Pedigree: ID SIRE DAM (no header)\n";
    cout << "  --genotyped-list <f>   File listing genotyped individual IDs\n";
    cout << "  --write-ebv            C: reserved Method-1 output switch, writer not implemented\n";
    cout << "  --h2-tau <f>          Genomic inverse scaling (default: 1.0)\n";
    cout << "  --h2-omega <f>        A22 inverse correction scaling (default: 1.0)\n";
    cout << "  --blend-g <f>         G blend weight with I (default: 0.05)\n";
    cout << "  Method 1 estimates VCs on the genotyped H subset, then scans fixed SNPs.\n";
    cout << "\n--- Method 2: full-pedigree ssGBLUP back-solving ---\n";
    cout << "  --single-step-method lmm|backsolve\n";
    cout << "                        lmm       = Method 1 fixed-SNP association (default, backward compatible)\n";
    cout << "                        backsolve = Method 2: ssGBLUP genetic evaluation -> GEBV -> SNP\n";
    cout << "                                    back-solving -> SNP variance sums/window summaries\n";
    cout << "                        A reference: BED, fixed positive Va/Ve, VanRaden, LDLT; P default off.\n";
    cout << "                        PGEN/BGEN backsolve rejected; VC/PCG/operator extensions are B.\n";
    cout << "  --reml-vg <f> --reml-ve <f>\n";
    cout << "                        Fixed variance components for Method 2 (canonical/deterministic path)\n";
    cout << "  --ssgwas-vce <mode>   VCE estimation mode for Method 2 (ai|emai|em|hi|exact|...).\n";
    cout << "                        If omitted, --reml-vg/--reml-ve are required (fixed-VC path)\n";
    cout << "  --ssgwas-window-snps <n>\n";
    cout << "                        SNP-count windows (default: 10); denominator total SNP variance,\n";
    cout << "                        not model Va or LD-adjusted window variance.\n";
    cout << "  --ssgwas-tune-g       Align G to A22 (G = b*G + a, preGSf90-style) before inversion\n";
    cout << "  --ssgwas-blend-a22    Blend G* with A22 instead of I (HIBLUP/preGSf90 reference; backsolve only)\n";
    cout << "  --ssgwas-solver ldlt|pcg\n";
    cout << "                        MME solver (default: ldlt reference; pcg = matrix-free PCG)\n";
    cout << "  --ssgwas-pcg-tol <x>  Matrix-free MME PCG tolerance (default: 1e-8)\n";
    cout << "  --ssgwas-pcg-maxit <n> Matrix-free MME PCG maximum iterations (default: 2000)\n";
    cout << "  --ssgwas-pvalue-method off|exact|apy-exact|apy-approx\n";
    cout << "                        Frequentist P-values for back-solved SNP effects (Phases C1/C2/D).\n";
    cout << "                        Default off. exact = C1 reference (selected MME inverse via LDLT,\n";
    cout << "                        O(n_genotyped^2) memory). apy-exact = APY route (Leite et al.\n";
    cout << "                        2024 GSE 56:59): APY G^-1 in the MME + exact core-slice selected\n";
    cout << "                        inverse (requires --ssgwas-solver ldlt; O(n_core^2) uncertainty\n";
    cout << "                        block; choose an explicit estimator below).\n";
    cout << "  --ssgwas-apy-exact-estimator <core-v2|legacy-hybrid-v1>\n";
    cout << "                        core-v2: matched projected/core SNP effects and exact selected-core\n";
    cout << "                        covariance (Bermann 2022 Eq. 21-22), kappa=1/S. Initial scope:\n";
    cout << "                        fixed VC, tau=omega=1, no tuning, identity blend 0.05, dense Gcc\n";
    cout << "                        solver and fixed core selection. B: limited validation, no universal calibration.\n";
    cout << "                        Default legacy-hybrid-v1 preserves full-point/core-variance outputs;\n";
    cout << "                        C compatibility: unmatched pairing, not supported calibrated inference.\n";
    cout << "                        apy-approx = Phase D Approx_GinvAPY (Leite Eq. 8/9):\n";
    cout << "                        APY G^-1 + approximate core PEV C_core = (D_c + lambda*Gcc^-1 +\n";
    cout << "                        Q*S*Q')^-1 (Schur/Woodbury form of Eq. 8; D = Eq. 9 diagonal\n";
    cout << "                        effective records); NO selected MME inverse and NO dense n_g^2\n";
    cout << "                        uncertainty block; works with --ssgwas-solver pcg; status\n";
    cout << "                        EXPERIMENTAL_APPROXIMATE: raw P-values, no genomic control,\n";
    cout << "                        strict type-I control NOT guaranteed. Use for research/exploration,\n";
    cout << "                        not as a calibrated inferential P-value path.\n";
    cout << "  --ssgwas-apy-core-size <n>\n";
    cout << "                        Core animals for the APY G^-1 (apy-* modes; >= 1, <= n_genotyped;\n";
    cout << "                        n = n_genotyped reproduces the C1 exact path)\n";
    cout << "                        'auto' = C developer-only sizing; excluded from v2 supported inference\n";
    cout << "                        due to the trace defect; no C98/C99 coverage guarantee.\n";
    cout << "  --ssgwas-apy-core-file <file>\n";
    cout << "                        Core-animal ID file (one ID per line, pedigree IDs); overrides\n";
    cout << "                        the seeded random core selection\n";
    cout << "  --ssgwas-apy-seed <s>  Seed of the deterministic random core selection (default 12345)\n";
    cout << "  --ssgwas-apy-core-target-ratio <x>\n";
    cout << "                        Harness-assigned label core_size / K98 (recorded in meta only; 0 = unset)\n";
    cout << "  --ssgwas-apy-core-variance <x>\n";
    cout << "                        Target proportion of genomic variance for automatic core sizing\n";
    cout << "                        (C98=0.98, C99=0.99; used when --ssgwas-apy-core-size auto)\n";
    cout << "  --ssgwas-apy-core-eigensolver <method>\n";
    cout << "                        Eigensolver for automatic core sizing: 'auto' (dense for small,\n";
    cout << "                        Lanczos for large), 'exact' (dense always), 'lanczos' (matrix-free)\n";
    cout << "  --ssgwas-grm-mode <mode>\n";
    cout << "                        VanRaden G construction: 'auto' (dense for small n_g, block-streamed\n";
    cout << "                        for large), 'dense' (full n_g^2 G), 'streamed' (block-streamed,\n";
    cout << "                        never materializes the full G; APY methods only)\n";
    cout << "  --ssgwas-apy-inverse-mode <dense|operator>\n";
    cout << "                        APY effect inverse application backend (default dense); operator\n";
    cout << "                        streams Gcn blocks and does not materialize Q in the operator\n";
    cout << "  --ssgwas-mme-apy-mode <dense|operator>\n";
    cout << "                        APY precision backend used by the MME (default dense); operator\n";
    cout << "                        requires --ssgwas-solver pcg and fixed variance components\n";
    cout << "  --ssgwas-pev-apy-mode <dense|operator>\n";
    cout << "                        APY Eq. 8 core PEV backend (default dense); operator streams\n";
    cout << "                        Gcn and avoids Q, explicit Gcc^-1, Kc^-1, and core Vhat\n";
    cout << "  --ssgwas-hinv-mode <dense|operator>\n";
    cout << "                        Pedigree A22^-1 backend (default dense); operator uses the\n";
    cout << "                        exact sparse Schur action and requires PCG + fixed VC\n";
    cout << "  --ssgwas-pev-scan-mode <scalar|batched>\n";
    cout << "                        SNP uncertainty scan backend (default batched; scalar is reference)\n";
    cout << "                        Exact/core-v2 use batched only when explicitly requested; omission\n";
    cout << "                        preserves historical scalar execution. Legacy hybrid stays scalar.\n";
    cout << "  --ssgwas-pev-scan-block-size <n>\n";
    cout << "                        Marker block width for batched uncertainty scan (default 256)\n";
    cout << "  --allow-large-exact-ssgwas-pvalue\n";
    cout << "                        Bypass the 8 GiB selected-inverse memory guard of\n";
    cout << "                        --ssgwas-pvalue-method exact (and apy-exact core slice)\n";
    cout << "  Method 2 outputs: <out>.ss_backsolve.ebv.txt / .snp.txt / .window.txt / .meta.txt\n";
    cout << "  With --ssgwas-pvalue-method exact additionally:\n";
    cout << "    <out>.ss_backsolve.pvalue.txt (EFFECT_ALLELE / SE_ALLELE / Z / P / LOG10P)\n";
    cout << "\n--- Pedigree Animal Model (B experimental) ---\n";
    cout << "  --animal-model        Pedigree-based Simple Animal Model GWAS: y = Xb + x*beta + Za + e,\n";
    cout << "                        a ~ N(0, A*sigma2_a) with full-pedigree A22 for the analysis sample set.\n";
    cout << "                        Requires --ped. A is built from the pedigree only (no genomic GRM);\n";
    cout << "                        LOCO does not apply (pedigree kernel is the same for all chromosomes).\n";
    cout << "                        Outputs *.assoc.txt (standard format) and *.animal_model.null.txt.\n";
    cout << "  --write-ped-kernel    Write the analysis A22 kernel (*.animal_model.k.txt + .k.id) for\n";
    cout << "                        external benchmarking (e.g. GEMMA -k).\n";
    cout << "\n--- Threshold/Probit Model (B experimental) ---\n";
    cout << "  --threshold-model      Binary trait GWAS (PQL + probit/logit)\n";
    cout << "  --threshold-dichotomize  Allow --threshold-model to dichotomise a non-binary\n";
    cout << "                         response at 0.5 (off by default; binary input is the contract)\n";
    cout << "  --link probit|logit    Link function (default: probit)\n";
    cout << "  PQL with fixed Ve=1 (probit) or Ve=pi^2/3 (logit).\n";
    cout << "  SPA is not a universal calibration guarantee; null-summary convergence flag is unreliable.\n";
    cout << "\n--- Count/Binomial PSR-GWAS Model (B experimental) ---\n";
    cout << "  --count-model          Count/binomial trait GWAS (PSR/PQL + weighted LMM)\n";
    cout << "  --family poisson|negbin|binomial Distribution (default: poisson)\n";
    cout << "  --offset <file>        One numeric log-offset/exposure value per sample for count traits\n";
    cout << "  --binomial-pheno <file> Binomial FID IID EVENT TRIAL file for event/trial traits\n";
    cout << "\n--- Ordinal PSR-GWAS Model (B experimental) ---\n";
    cout << "  --ordinal-model        Ordinal trait GWAS (cumulative-binary or stacked PSR)\n";
    cout << "  --ordinal-pheno <file> Ordinal FID IID ORDINAL phenotype file\n";
    cout << "  --ordinal-mode <name>  experimental-cumulative-binary (default) or stacked-psr\n";
    cout << "  --ordinal-stacked-engine <name> fixed-gls (default) or dense-zkz\n";
    cout << "  --ordinal-vc-method <name> fixed (default) or grid-reml\n";
    cout << "  --ordinal-test <name>  wald (default), score, robust-score, or all\n";
    cout << "  --ordinal-allow-fallback Fall back from stacked-psr to cumulative-binary\n";
    cout << "  --psr-null-only        Output null model (Vg, Ve, fixed effects, BLUP, quasi-likelihood)\n";
    cout << "                        then exit without SNP scan (available for threshold/count/ordinal)\n";
    cout << "  --psr-maxit <n>       Max PSR/PQL iterations (default: 50)\n";
    cout << "  --psr-tol <f>         PSR/PQL convergence tolerance (default: 1e-6)\n";
    cout << "  --psr-mu-eps <f>      Mu clipping epsilon for PSR (default: 1e-6)\n";
    cout << "  --psr-min-weight <f>  Minimum PSR working weight (default: 1e-8)\n";
    cout << "  --psr-max-weight <f>  Maximum PSR working weight (default: 1e8)\n";
    cout << "\n--- Survival prototype (C developer-only, not supported inference) ---\n";
    cout << "  --survival-model       Incomplete residual-like LMM prototype; no fitted Cox null\n";
    cout << "  --pheno format: FID IID TIME STATUS  (2 columns after FID/IID)\n";
    cout << "  runCoxNull only prints a message; do not use for screening or inference.\n";
    cout << "\n--- GxE Reaction Norm Model (B experimental) ---\n";
    cout << "  --gxe-model            SNP-by-environment interaction GWAS\n";
    cout << "  --env-covar <col>      Environment variable column in --pheno\n";
    cout << "  VCE: Vg0 (intercept) + Vg1 (slope on env) + residual.\n";
    cout << "  Scan: 2-df Wald test for beta_g + beta_gxe jointly.\n";
    cout << "  Output: BETA_G, SE_G, BETA_GxE, SE_GxE, JOINT_CHISQ, JOINT_P\n";
    cout << "\n=== INFRASTRUCTURE ===\n";
    cout << "  --make-grm             Calculate GRM and exit\n";
    cout << "  -gk [n]                GEMMA-compatible GRM alias; optional algorithm\n";
    cout << "  --make-grm-sparse      PLINK-only sparse VanRaden GRM and exit\n";
    cout << "  --sparse-cutoff <f>    Sparse GRM absolute-value cutoff (default: 0.05)\n";
    cout << "  --gs/-gs <n>           GEMMA-style S matrix for summary-stat VC (mode 1)\n";
    cout << "  --gq/-gq <n>           GEMMA-style q vector from --beta summary stats (mode 1)\n";
    cout << "  --gemma-eigen/-eigen   Calculate GRM eigen decomposition and write eigenD/eigenU\n";
    cout << "  --calccor/-calccor/-gw GEMMA-style SNP correlation windows; writes <out>.cor.txt\n";
    cout << "  --windowbp <n>         BP window for --calccor (default: 1000000; 0=no BP limit)\n";
    cout << "  --windowcm <f>         cM window for --calccor (default: 0; PLINK BIM only)\n";
    cout << "  --windowns <n>         SNP-count window for --calccor (default: 0=no SNP-count limit)\n";
    cout << "  --vc/-vc <n>           GEMMA-style VC estimate: 1 HE, 2 REML, 3 AI/IDUL REML\n";
    cout << "  --reml, --vce          Standalone production REML/VCE estimation\n";
    cout << "  --reml-method <m>      ai, em, brent, idul, idul_plus, or aidul\n";
    cout << "  --ci/-ci <n>           GEMMA-style summary-stat PVE confidence interval\n";
    cout << "  --beta/-beta <file>    Summary statistics for --gq or direct --vc with --cor\n";
    cout << "  --cor/-cor <file>      GEMMA .cor.txt for direct summary-stat --vc with --beta\n";
    cout << "  --cat/-cat <file>      SNP category file for summary-stat routes\n";
    cout << "  --mcat/-mcat <file>    List of SNP category files for summary-stat routes\n";
    cout << "  --wcat/-wcat <file>    SNP-by-category weights for summary-stat routes\n";
    cout << "  --wsnp/-wsnp <file>    SNP weights for summary-stat routes\n";
    cout << "  --study/-study <p>     Study prefix for summary-stat --vc (<p>.Vq/q/size.txt)\n";
    cout << "  --ref/-ref <p>         Reference prefix for summary-stat --vc (<p>.S/size.txt)\n";
    cout << "  --mstudy/-mstudy <f>   File listing study prefixes for summary-stat --vc/--ci\n";
    cout << "  --mref/-mref <f>       File listing reference prefixes for summary-stat --vc/--ci\n";
    cout << "  --make-grm-alg <n>     Method: 1=VanRaden, 2=Zeng, 3=Yang, 4=Vitezica,\n";
    cout << "                          5=TrueVar, 6=EMmvGWAS-compatible (default: 1)\n";
    cout << "                          Note: alg 6 is auto-set by --mv-compatible emmvgwas\n";
    cout << "  --mass-regress <file>  Multi-phenotype mass regressor (testing)\n";
    cout << "  --biobank-mode         B: native REGENIE-style two-step research path, not reference equivalence\n";
    cout << "  --regenie-model-snps <f> Model SNP list for Step 1\n";
    cout << "  --write-pred <file>    Write step 1 LOCO predictions\n";
    cout << "  --use-pred <file>      Use existing LOCO predictions\n";
    cout << "  --regenie-step <n>     0=both, 1=predictions only, 2=association only\n";
    cout << "  --predict-effect/-epm <f> GEMMA-style prediction from SNP effect table\n";
    cout << "  --predict-score <f>    Output prediction file (default: <out>.predict.txt)\n";
    cout << "  --predict-snp-col <c>  SNP ID column in effect table (default: SNP/rs/id)\n";
    cout << "  --predict-effect-col <c> Effect column (default: BETA/beta/alpha/effect)\n";
    cout << "  --bslmm/-bslmm <n>     B: native/compatibility modes: sparse ridge, ridge, probit, DAP-style\n";
    cout << "                          Native sampling is separate; ridge parity does not validate posterior parity.\n";
    cout << "  --bslmm-top <n>        Sparse top-SNP count for --bslmm 1/3 approximation (default: 100)\n";
    cout << "  --bslmm-h2 <f>         Ridge/GBLUP shrinkage h2 in (0,1), default 0.5\n";
    cout << "  --bslmm-mcmc           Native posterior-sampling BSLMM\n";
    cout << "  --bvsr                 Native posterior-sampling BVSR\n";
    cout << "  --mcmc-burnin <n>      MCMC burn-in iterations\n";
    cout << "  --mcmc-samples <n>     MCMC post-burn-in iterations\n";
    cout << "  --mcmc-thin <n>        MCMC thinning interval\n";
    cout << "  --seed <n>             Deterministic BSLMM/BVSR seed\n";
    cout << "  --bf/-bf <file>        DAP BF file for --bslmm 5\n";
    cout << "  --hyp/-hyp <file>      DAP hyperparameter file for --bslmm 5\n";
    cout << "  --bslmm-hgrid <n>      DAP h grid size for --bslmm 4 (default: 5)\n";
    cout << "  --bslmm-rho-grid <n>   DAP rho grid size for --bslmm 4 (default: 5)\n";
    cout << "  -hmin/-hmax <f>        GEMMA BSLMM/DAP h prior range\n";
    cout << "  -rmin/-rmax <f>        GEMMA BSLMM/DAP rho prior range\n";
    cout << "  -pmin/-pmax <f>        GEMMA BSLMM/DAP log10(pi) prior range\n";
    cout << "  -gmean <f>             GEMMA BSLMM gamma-count prior mean hint\n";
    cout << "  -hscale/-rscale/-pscale <f> GEMMA BSLMM proposal scale hints\n";
    cout << "  -mh <n>                GEMMA BSLMM MH step-count hint\n";
    cout << "  -smin/-smax <n>        GEMMA sparse gamma size bounds (-smax aliases --bslmm-top)\n";
    cout << "  --block-size <n>       Ridge block size (default: 1000)\n";
    cout << "  --skat                 Scoped set test with tail approximation, not exact-tail equivalence\n";
    cout << "  --burden               Run Burden rare-variant test (--set-file)\n";
    cout << "  --set-file <file>      Variant set definition for SKAT/Burden\n";
    cout << "\nModel notes:\n";
    cout << "  Models are mutually exclusive. Use exactly one --*-model flag.\n";
    cout << "  Default (no flag): standard univariate LMM (single-trait).\n";
    cout << "  Report bugs at: https://github.com/anomalyco/PanOmics/issues\n";
}

int main(int argc, char* argv[]) {
#ifdef _WIN32
    // Source messages are UTF-8; Windows consoles otherwise often use a legacy
    // code page. Redirected result files retain their existing byte format.
    SetConsoleOutputCP(CP_UTF8);
#endif

    try {
    const auto invocation = prepareCliInvocation(argc, argv);
    if (invocation.help) {
        if (invocation.help_topic == "all" || invocation.help_topic == "advanced") print_usage(argv[0]);
        else printCliHelp(invocation.help_topic, argv[0]);
        return argc < 2 ? 1 : 0;
    }
    auto normalized_args = invocation.args;
    std::vector<char*> normalized_argv;
    for (auto& value : normalized_args) normalized_argv.push_back(value.data());
    argc = static_cast<int>(normalized_argv.size());
    argv = normalized_argv.data();

    GwasConfig cfg;
    std::string selected_cli_model;
    bool bgen_dom_qc_requested = false;
    cfg.command_argv = invocation.original;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--input-config" && i + 1 < argc) {
            apply_input_config(cfg, argv[++i]);
        }
    }

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--pca") { cfg.pca = true; }
        else if (arg == "--pca-k" && i + 1 < argc) { cfg.pca_k = stoi(argv[++i]); }
        else if (arg == "--geno" && i + 1 < argc) { cfg.geno = stod(argv[++i]); cfg.qc_requested = true; }
        else if (arg == "--mind" && i + 1 < argc) { cfg.mind = stod(argv[++i]); cfg.qc_requested = true; }
        else if (arg == "--keep" && i + 1 < argc) { cfg.keep_file = argv[++i]; cfg.qc_requested = true; }
        else if (arg == "--remove" && i + 1 < argc) { cfg.remove_file = argv[++i]; cfg.qc_requested = true; }
        else if (arg == "--ld-prune") { cfg.ld_prune = true; }
        else if (arg == "--ld-clump" && i + 1 < argc) { cfg.ld_clump_file = argv[++i]; }
        else if (arg == "--ld-window" && i + 1 < argc) { cfg.ld_window = stoi(argv[++i]); }
        else if (arg == "--ld-step" && i + 1 < argc) { cfg.ld_step = stoi(argv[++i]); }
        else if (arg == "--ld-r2" && i + 1 < argc) { cfg.ld_r2 = stod(argv[++i]); }
        else if (arg == "--clump-p" && i + 1 < argc) { cfg.clump_p = stod(argv[++i]); }
        else if (arg == "--skat-o") { cfg.skat_o = true; cfg.run_gwas = true; }
        else if (arg == "--acat-v") { cfg.acat_v = true; cfg.run_gwas = true; }
        else if (arg == "--rare-max-maf" && i + 1 < argc) { cfg.rare_max_maf = stod(argv[++i]); }
        else if (arg == "--rare-min-mac" && i + 1 < argc) { cfg.rare_min_mac = stoi(argv[++i]); }
        else if (arg == "--snp-model" && i + 1 < argc) { cfg.snp_model = argv[++i]; }
        else if (arg == "--gwas-model" && i + 1 < argc) { cfg.gwas_model = argv[++i]; }
        else if (arg == "--farmcpu") { cfg.farmcpu = true; cfg.run_gwas = true; }
        else if (arg == "--farmcpu-max-iterations" && i + 1 < argc) { cfg.farmcpu_max_iterations = stoi(argv[++i]); }
        else if (arg == "--farmcpu-pseudo-qtn" && i + 1 < argc) { cfg.farmcpu_pseudo_qtn = stoi(argv[++i]); }
        else if (arg == "--blink") { cfg.blink = true; cfg.run_gwas = true; }
        else if (arg == "--condition-snp" && i + 1 < argc) { cfg.condition_snp_file = argv[++i]; cfg.run_gwas = true; }
        else if (arg == "--conditional-joint") { cfg.conditional_joint = true; cfg.run_gwas = true; }
        else if (arg == "--firth-logistic") { cfg.firth_logistic = true; cfg.run_gwas = true; }
        else if (arg == "--logistic" && i + 1 < argc) { cfg.logistic_test = argv[++i]; cfg.run_gwas = true; }
        else if (arg == "--logistic-glmm") { cfg.logistic_glmm = true; cfg.run_gwas = true; }
        else if (arg == "--saige-style") { cfg.saige_style = true; cfg.run_gwas = true; }
        else if (arg == "--cojo") { cfg.cojo = true; cfg.run_gwas = true; }
        else if (arg == "--cojo-slct") { cfg.cojo_slct = true; cfg.run_gwas = true; }
        else if (arg == "--cojo-joint") { cfg.cojo_joint = true; cfg.run_gwas = true; }
        else if (arg == "--cojo-summary" && i + 1 < argc) { cfg.cojo_summary = argv[++i]; }
        else if (arg == "--cojo-ld-ref" && i + 1 < argc) { cfg.cojo_ld_ref = argv[++i]; }
        else if (arg == "--cojo-p" && i + 1 < argc) { cfg.cojo_p = stod(argv[++i]); }
        else if (arg == "--cojo-r2" && i + 1 < argc) { cfg.cojo_r2 = stod(argv[++i]); }
        else if (arg == "--cojo-wind" && i + 1 < argc) { cfg.cojo_window = stoi(argv[++i]); }
        else if (arg == "--cojo-cond" && i + 1 < argc) { cfg.cojo_cond_file = argv[++i]; }
        else if (arg == "--fine-map") { cfg.fine_map = true; cfg.run_gwas = true; }
        else if (arg == "--finemap-summary" && i + 1 < argc) { cfg.finemap_summary = argv[++i]; }
        else if (arg == "--finemap-region" && i + 1 < argc) { cfg.finemap_region = argv[++i]; }
        else if (arg == "--credible-level" && i + 1 < argc) { cfg.credible_level = stod(argv[++i]); }
        else if (arg == "--finemap-prior-var" && i + 1 < argc) { cfg.finemap_prior_var = stod(argv[++i]); }
        else if (arg == "--finemap-max-signals" && i + 1 < argc) { cfg.finemap_max_signals = stoi(argv[++i]); }
        else if (arg == "--finemap-signal-p" && i + 1 < argc) { cfg.finemap_signal_p = stod(argv[++i]); }
        else if (arg == "--meta-analysis") { cfg.meta_analysis = true; cfg.run_gwas = true; }
        else if (arg == "--meta-studies" && i + 1 < argc) { cfg.meta_studies_file = argv[++i]; }
        else if (arg == "--meta-method" && i + 1 < argc) { cfg.meta_method = argv[++i]; }
        else if (arg == "--meta-het-p" && i + 1 < argc) { cfg.meta_het_p = stod(argv[++i]); }
        else if (arg == "--meta-min-studies" && i + 1 < argc) { cfg.meta_min_studies = stoi(argv[++i]); }
        else if (arg == "--meta-gc") { cfg.meta_gc = true; }
        else if (arg == "--meta-gc-lambda" && i + 1 < argc) { cfg.meta_gc_lambda = stod(argv[++i]); }
        else if (arg == "--meta-out-format" && i + 1 < argc) { cfg.meta_out_format = argv[++i]; }
        else if (arg == "--bfile" && i + 1 < argc) { cfg.bfile_prefix = argv[++i]; }
        else if (arg == "--bgen" && i + 1 < argc) { cfg.bgen_file = argv[++i]; }
        else if (arg == "--sample" && i + 1 < argc) { cfg.sample_file = argv[++i]; }
        else if (arg == "--pfile" && i + 1 < argc) { cfg.pfile_prefix = argv[++i]; }
        else if (arg == "--pgen" && i + 1 < argc) { cfg.pgen_file = argv[++i]; }
        else if (arg == "--pvar" && i + 1 < argc) { cfg.pvar_file = argv[++i]; }
        else if (arg == "--psam" && i + 1 < argc) { cfg.psam_file = argv[++i]; }
        else if (arg == "--vcf" && i + 1 < argc) { cfg.vcf_file = argv[++i]; }
        else if (arg == "--vcf-tbi" && i + 1 < argc) { cfg.vcf_tbi_file = argv[++i]; }
        else if (arg == "--vcf-sample" && i + 1 < argc) { cfg.vcf_sample_file = argv[++i]; }
        else if (arg == "--hmp" && i + 1 < argc) { cfg.hmp_file = argv[++i]; }
        else if (arg == "--geno-matrix" && i + 1 < argc) { cfg.geno_matrix_file = argv[++i]; }
        else if (arg == "--geno-var" && i + 1 < argc) { cfg.geno_var_file = argv[++i]; }
        else if (arg == "--geno-sample" && i + 1 < argc) { cfg.geno_sample_file = argv[++i]; }
        else if (arg == "--geno-orientation" && i + 1 < argc) { cfg.geno_orientation = argv[++i]; }
        else if (arg == "--npy" && i + 1 < argc) { cfg.npy_file = argv[++i]; }
        else if (arg == "--npy-var" && i + 1 < argc) { cfg.npy_var_file = argv[++i]; }
        else if (arg == "--npy-sample" && i + 1 < argc) { cfg.npy_sample_file = argv[++i]; }
        else if (arg == "--npy-orientation" && i + 1 < argc) { cfg.npy_orientation = argv[++i]; }
        else if (arg == "--ped" && i + 1 < argc) { cfg.ped_file = argv[++i]; }
        else if (arg == "--map" && i + 1 < argc) { cfg.map_file = argv[++i]; }
        else if (arg == "--convert" && i + 1 < argc) { cfg.convert_from = argv[++i]; }
        else if (arg == "--to" && i + 1 < argc) { cfg.convert_to = argv[++i]; }
        else if (arg == "--batch-models" && i + 1 < argc) { cfg.batch_models = argv[++i]; }
        else if (arg == "--batch-out" && i + 1 < argc) { cfg.batch_out_prefix = argv[++i]; }
        else if ((arg == "--out" || arg == "-o") && i + 1 < argc) { cfg.out_prefix = argv[++i]; }
        else if (arg == "--out-ssf") { cfg.out_ssf = true; }
        else if (arg == "--out-gz") { cfg.out_gz = true; }
        else if (arg == "--manhattan-tsv") { cfg.manhattan_tsv = true; }
        else if (arg == "--qq-tsv") { cfg.qq_tsv = true; }
        else if ((arg == "--pheno" || arg == "-p") && i + 1 < argc) { cfg.pheno_file = argv[++i]; }
        else if (arg == "--input-file" && i + 1 < argc) { cfg.input_file = argv[++i]; }
        else if (arg == "--input-config" && i + 1 < argc) { cfg.input_config = argv[++i]; }
        else if (arg == "--write-template-config" && i + 1 < argc) { cfg.write_template_config = argv[++i]; }
        else if (arg == "--input-format" && i + 1 < argc) { cfg.input_format = argv[++i]; }
        else if (arg == "--input-model" && i + 1 < argc) { cfg.input_model = argv[++i]; }
        else if (arg == "--prepared-input-dir" && i + 1 < argc) { cfg.prepared_input_dir = argv[++i]; }
        else if (arg == "--prepare-input-only") { cfg.prepare_input_only = true; }
        else if (arg == "--id-col" && i + 1 < argc) { cfg.id_col = argv[++i]; }
        else if (arg == "--iid-col" && i + 1 < argc) { cfg.iid_col = argv[++i]; }
        else if (arg == "--time-col" && i + 1 < argc) { cfg.time_col = argv[++i]; }
        else if (arg == "--trait-cols" && i + 1 < argc) { cfg.trait_cols = argv[++i]; }
        else if (arg == "--input-trait" && i + 1 < argc) { cfg.input_trait = argv[++i]; }
        else if (arg == "--record-qcovars" && i + 1 < argc) { cfg.record_qcovars = argv[++i]; }
        else if (arg == "--record-ccovars" && i + 1 < argc) { cfg.record_ccovars = argv[++i]; }
        else if (arg == "--subject-qcovars" && i + 1 < argc) { cfg.subject_qcovars = argv[++i]; }
        else if (arg == "--subject-ccovars" && i + 1 < argc) { cfg.subject_ccovars = argv[++i]; }
        else if (arg == "--time-parser" && i + 1 < argc) { cfg.time_parser = argv[++i]; }
        else if (arg == "--time-unique-policy" && i + 1 < argc) { cfg.time_unique_policy = argv[++i]; }
        else if (arg == "--time-index-multiplier" && i + 1 < argc) { cfg.time_index_multiplier = argv[++i]; }
        else if (arg == "--preview-schema") { cfg.preview_schema = true; }
        else if (arg == "--infer-columns") { cfg.infer_columns = true; }
        else if (arg == "--wide-to-long") { cfg.wide_to_long = true; }
        else if (arg == "--long-pheno" && i + 1 < argc) { cfg.long_pheno_file = argv[++i]; cfg.rrm_model = true; }
        else if ((arg == "--mpheno" || arg == "-n") && i + 1 < argc) { cfg.mpheno = stoi(argv[++i]); }
        else if (arg == "--tpheno" && i + 1 < argc) { cfg.tpheno = stoi(argv[++i]); cfg.tpheno_set = true; }
        else if (arg == "--rrm-order" && i + 1 < argc) { cfg.rrm_order = stoi(argv[++i]); }
        else if (arg == "--rrm-basis" && i + 1 < argc) { cfg.rrm_basis = argv[++i]; }
        else if (arg == "--rrm-curve-grid" && i + 1 < argc) { cfg.rrm_curve_grid = stoi(argv[++i]); }
        else if (arg == "--rrm-top" && i + 1 < argc) { cfg.rrm_top = stoi(argv[++i]); }
        else if (arg == "--rrm-ai-iter" && i + 1 < argc) { cfg.rrm_ai_iter = stoi(argv[++i]); }
        else if (arg == "--rrm-sig-threshold" && i + 1 < argc) { cfg.rrm_sig_threshold = stod(argv[++i]); }
        else if (arg == "--rrm-unbalance") { cfg.rrm_unbalance = true; cfg.rrm_model = true; cfg.run_gwas = true; }
        else if (arg == "--rrm-gwas-method" && i + 1 < argc) { cfg.rrm_gwas_method = argv[++i]; }
        else if (arg == "--rrm-loco-global-null") { cfg.rrm_loco_global_null = true; }
        else if (arg == "--tpoint" && i + 1 < argc) { cfg.tpoint_col_name = argv[++i]; }
        else if (arg == "--tclass" && i + 1 < argc) { cfg.tclass_names = argv[++i]; }
        else if (arg == "--tcovar" && i + 1 < argc) { cfg.tcovar_names = argv[++i]; }
        else if (arg == "--rrm-compatible" && i + 1 < argc) { cfg.rrm_compatible = argv[++i]; cfg.rrm_model = true; }
        else if (arg == "--rrm-compatible-varcom" && i + 1 < argc) { cfg.rrm_compatible_varcom = argv[++i]; }
        else if (arg == "--mpheno-list" && i + 1 < argc) {
            string s = argv[++i];
            stringstream ss(s);
            string token;
            while(getline(ss, token, ',')) { cfg.mpheno_list.push_back(stoi(token)); }
        }
        else if (arg == "--multivariate") { cfg.multivariate = true; cfg.run_gwas = true; }
        else if (arg == "--emmvgwas" || arg == "--emmvgwas-full-gwas" || arg == "--emmvgwas-semi-gwas") {
            cfg.multivariate = true;
            cfg.run_gwas = true;
            cfg.multivariate_method = "emmv";
            cfg.multivariate_compatible = "emmvgwas";
            cfg.no_loco = true;
        }
        else if (arg == "--emmvgwas-vc-only" || arg == "--emmvgwas-vc-estimation") {
            cfg.multivariate = true;
            cfg.run_gwas = true;
            cfg.multivariate_method = "emmv";
            cfg.multivariate_compatible = "emmvgwas";
            cfg.multivariate_emmvgwas_output = true;
            cfg.multivariate_emmvgwas_vc_only = true;
            cfg.no_loco = true;
        }
        else if (arg == "--mv-method" && i + 1 < argc) { cfg.multivariate_method = argv[++i]; cfg.multivariate_method_explicit = true; }
        else if (arg == "--mv-compatible" && i + 1 < argc) { cfg.multivariate_compatible = argv[++i]; }
        else if (arg == "--mv-test" && i + 1 < argc) { cfg.multivariate_tests = argv[++i]; }
        else if (arg == "--mv-gemma-mode" && i + 1 < argc) { cfg.multivariate_gemma_mode = stoi(argv[++i]); }
        else if (arg == "--mv-pnr" && i + 1 < argc) { cfg.multivariate_pnr = stod(argv[++i]); }
        else if (arg == "--mv-emi" && i + 1 < argc) { cfg.multivariate_emi = stoi(argv[++i]); }
        else if (arg == "--mv-nri" && i + 1 < argc) { cfg.multivariate_nri = stoi(argv[++i]); }
        else if (arg == "--mv-emp" && i + 1 < argc) { cfg.multivariate_emp = stod(argv[++i]); }
        else if (arg == "--mv-nrp" && i + 1 < argc) { cfg.multivariate_nrp = stod(argv[++i]); }
        else if (arg == "--mv-crt") { cfg.multivariate_crt = true; }
        else if (arg == "--mv-refine-top" && i + 1 < argc) { cfg.multivariate_refine_top = stoi(argv[++i]); }
        else if (arg == "--mv-refine-p" && i + 1 < argc) { cfg.multivariate_refine_p = stod(argv[++i]); }
        else if (arg == "--mv-refine-discrepancy-top" && i + 1 < argc) { cfg.multivariate_refine_discrepancy_top = stoi(argv[++i]); }
        else if (arg == "--mv-refine-discrepancy" && i + 1 < argc) { cfg.multivariate_refine_discrepancy = stod(argv[++i]); }
        else if (arg == "--mv-refine-score-prepass" && i + 1 < argc) { cfg.multivariate_refine_score_prepass = stoi(argv[++i]); }
        else if (arg == "--mv-refit-p" && i + 1 < argc) { cfg.multivariate_refit_p = stod(argv[++i]); }
        else if (arg == "--mv-vce-solver" && i + 1 < argc) { cfg.multivariate_use_aireml = (string(argv[++i]) == "aireml"); }
        else if (arg == "--mv-em-min-iter" && i + 1 < argc) { cfg.multivariate_em_min_iter = stoi(argv[++i]); }
        else if (arg == "--mv-ai-max-iter" && i + 1 < argc) { cfg.multivariate_ai_max_iter = stoi(argv[++i]); }
        else if (arg == "--mv-ai-tol" && i + 1 < argc) { cfg.multivariate_ai_tol = stod(argv[++i]); }
        else if (arg == "--mv-output-vbeta") { cfg.multivariate_output_vbeta = true; }
        else if (arg == "--mv-emmvgwas-output" || arg == "--emmvgwas-output") { cfg.multivariate_emmvgwas_output = true; }
        else if ((arg == "--mv-emmvgwas-G" || arg == "--emmvgwas-G") && i + 1 < argc) { cfg.multivariate_emmvgwas_G_file = argv[++i]; }
        else if ((arg == "--mv-emmvgwas-R" || arg == "--emmvgwas-R") && i + 1 < argc) { cfg.multivariate_emmvgwas_R_file = argv[++i]; }
        else if (arg == "--repeatability-model") { cfg.repeatability_model = true; cfg.run_gwas = true; }
        else if (arg == "--maternal-model") { cfg.maternal_model = true; cfg.run_gwas = true; }
        else if (arg == "--maternal-id" && i + 1 < argc) { cfg.maternal_id_col = stoi(argv[++i]); }
        else if (arg == "--maternal-snp-test" && i + 1 < argc) {
            std::string v = argv[++i];
            if (v != "direct" && v != "maternal" && v != "joint") {
                throw std::runtime_error("Invalid --maternal-snp-test value '" + v + "' (expected direct|maternal|joint)");
            }
            cfg.maternal_snp_test = v;
        }
        else if (arg == "--common-env-model") { cfg.common_env_model = true; cfg.run_gwas = true; }
        else if (arg == "--common-env-file" && i + 1 < argc) { cfg.common_env_file = argv[++i]; }
        else if (arg == "--common-env-record-level") { cfg.common_env_record_level = true; }
        else if (arg == "--single-step") { cfg.single_step_model = true; cfg.run_gwas = true; }
        else if (arg == "--animal-model") { cfg.animal_model = true; cfg.run_gwas = true; }
        else if (arg == "--write-ped-kernel") { cfg.write_ped_kernel = true; }
        else if (arg == "--make-agh") { cfg.single_step_model = true; cfg.make_agh = true; cfg.write_agh = true; cfg.run_gwas = true; }
        else if (arg == "--write-agh") { cfg.write_agh = true; }
        else if (arg == "--ped" && i + 1 < argc) { cfg.ped_file = argv[++i]; }
        else if (arg == "--genotyped-list" && i + 1 < argc) { cfg.genotyped_list = argv[++i]; }
        else if (arg == "--h2-tau" && i + 1 < argc) { cfg.h_tau = stod(argv[++i]); }
        else if (arg == "--h2-omega" && i + 1 < argc) { cfg.h_omega = stod(argv[++i]); }
        else if (arg == "--blend-g" && i + 1 < argc) { cfg.blend_g = stod(argv[++i]); }
        else if (arg == "--single-step-method" && i + 1 < argc) {
            std::string v = argv[++i];
            if (v != "lmm" && v != "backsolve") {
                throw std::runtime_error("Invalid --single-step-method value '" + v + "' (expected lmm|backsolve)");
            }
            cfg.single_step_method = v;
            cfg.single_step_model = true;
            cfg.run_gwas = true;
        }
        else if (arg == "--ssgwas-window-snps" && i + 1 < argc) { cfg.ssgwas_window_snps = stoi(argv[++i]); }
        else if (arg == "--ssgwas-tune-g") { cfg.ssgwas_tune_g = true; }
        else if (arg == "--ssgwas-blend-a22") { cfg.ssgwas_blend_a22 = true; }
        else if (arg == "--ssgwas-solver" && i + 1 < argc) {
            std::string v = argv[++i];
            if (v != "ldlt" && v != "pcg") {
                throw std::runtime_error("Invalid --ssgwas-solver value '" + v + "' (expected ldlt|pcg)");
            }
            cfg.ssgwas_solver = v;
        }
        else if (arg == "--ssgwas-pcg-tol" && i + 1 < argc) {
            cfg.ssgwas_pcg_tol = std::stod(argv[++i]);
            if (!(cfg.ssgwas_pcg_tol > 0.0)) {
                throw std::runtime_error("--ssgwas-pcg-tol must be > 0");
            }
        }
        else if (arg == "--ssgwas-pcg-maxit" && i + 1 < argc) {
            cfg.ssgwas_pcg_maxit = std::stoi(argv[++i]);
            if (cfg.ssgwas_pcg_maxit <= 0) {
                throw std::runtime_error("--ssgwas-pcg-maxit must be > 0");
            }
        }
        else if (arg == "--ssgwas-pvalue-method" && i + 1 < argc) {
            std::string v = argv[++i];
            if (v != "off" && v != "exact" && v != "apy-exact" && v != "apy-approx") {
                throw std::runtime_error("Invalid --ssgwas-pvalue-method value '" + v +
                                         "' (expected off|exact|apy-exact|apy-approx)");
            }
            cfg.ssgwas_pvalue_method = v;
        }
        else if (arg == "--ssgwas-apy-exact-estimator") {
            if (i + 1 >= argc) {
                throw std::runtime_error("--ssgwas-apy-exact-estimator requires core-v2|legacy-hybrid-v1");
            }
            std::string v = argv[++i];
            if (v != "core-v2" && v != "legacy-hybrid-v1") {
                throw std::runtime_error("Invalid --ssgwas-apy-exact-estimator value '" + v +
                                         "' (expected core-v2|legacy-hybrid-v1)");
            }
            cfg.ssgwas_apy_exact_estimator = v;
            cfg.ssgwas_apy_exact_estimator_explicit = true;
        }
        else if (arg == "--ssgwas-apy-core-size" && i + 1 < argc) {
            std::string v = argv[++i];
            if (v == "auto") {
                cfg.ssgwas_apy_core_size = -1;  // sentinel for automatic C98/C99
            } else {
                int vi = stoi(v);
                if (vi < 1 || vi > 50000000) {
                    throw std::runtime_error("Invalid --ssgwas-apy-core-size value '" + v +
                                             "' (expected >= 1 or 'auto')");
                }
                cfg.ssgwas_apy_core_size = vi;
            }
        }
        else if (arg == "--ssgwas-apy-core-file" && i + 1 < argc) { cfg.ssgwas_apy_core_file = argv[++i]; }
        else if (arg == "--ssgwas-apy-seed" && i + 1 < argc) { cfg.ssgwas_apy_seed = stoi(argv[++i]); }
        else if (arg == "--ssgwas-apy-core-target-ratio" && i + 1 < argc) {
            cfg.ssgwas_apy_core_target_ratio = stod(argv[++i]);
        }
        else if (arg == "--ssgwas-apy-core-variance" && i + 1 < argc) {
            double v = stod(argv[++i]);
            if (v <= 0.0 || v > 1.0) {
                throw std::runtime_error("--ssgwas-apy-core-variance must be in (0, 1]");
            }
            cfg.ssgwas_apy_core_variance = v;
        }
        else if (arg == "--ssgwas-apy-core-eigensolver" && i + 1 < argc) {
            std::string v = argv[++i];
            if (v != "auto" && v != "exact" && v != "lanczos") {
                throw std::runtime_error("--ssgwas-apy-core-eigensolver must be auto|exact|lanczos");
            }
            cfg.ssgwas_apy_core_eigensolver = v;
        }
        else if (arg == "--ssgwas-grm-mode" && i + 1 < argc) {
            std::string v = argv[++i];
            if (v != "auto" && v != "dense" && v != "streamed") {
                throw std::runtime_error("--ssgwas-grm-mode must be auto|dense|streamed");
            }
            cfg.ssgwas_grm_mode = v;
        }
        else if (arg == "--ssgwas-gcc-solver" && i + 1 < argc) {
            std::string v = argv[++i];
            if (v != "dense" && v != "pcg") {
                throw std::runtime_error("--ssgwas-gcc-solver must be dense|pcg");
            }
            cfg.ssgwas_gcc_solver = v;
        }
        else if (arg == "--ssgwas-apy-inverse-mode" && i + 1 < argc) {
            std::string v = argv[++i];
            if (v != "dense" && v != "operator") {
                throw std::runtime_error("--ssgwas-apy-inverse-mode must be dense|operator");
            }
            cfg.ssgwas_apy_inverse_mode = v;
        }
        else if (arg == "--ssgwas-mme-apy-mode" && i + 1 < argc) {
            std::string v = argv[++i];
            if (v != "dense" && v != "operator") {
                throw std::runtime_error("--ssgwas-mme-apy-mode must be dense|operator");
            }
            cfg.ssgwas_mme_apy_mode = v;
        }
        else if (arg == "--ssgwas-pev-apy-mode" && i + 1 < argc) {
            std::string v = argv[++i];
            if (v != "dense" && v != "operator") {
                throw std::runtime_error("--ssgwas-pev-apy-mode must be dense|operator");
            }
            cfg.ssgwas_pev_apy_mode = v;
        }
        else if (arg == "--ssgwas-hinv-mode" && i + 1 < argc) {
            std::string v = argv[++i];
            if (v != "dense" && v != "operator") {
                throw std::runtime_error("--ssgwas-hinv-mode must be dense|operator");
            }
            cfg.ssgwas_hinv_mode = v;
        }
        else if (arg == "--ssgwas-pev-scan-mode" && i + 1 < argc) {
            std::string v = argv[++i];
            if (v != "scalar" && v != "batched") {
                throw std::runtime_error("--ssgwas-pev-scan-mode must be scalar|batched");
            }
            cfg.ssgwas_pev_scan_mode = v;
            cfg.ssgwas_pev_scan_mode_explicit = true;
        }
        else if (arg == "--ssgwas-pev-scan-block-size" && i + 1 < argc) {
            int v = stoi(argv[++i]);
            if (v <= 0) throw std::runtime_error("--ssgwas-pev-scan-block-size must be > 0");
            cfg.ssgwas_pev_scan_block_size = v;
        }
        else if (arg == "--ssgwas-pcg-preconditioner" && i + 1 < argc) {
            std::string v = argv[++i];
            if (v != "jacobi" && v != "lowrank") {
                throw std::runtime_error("--ssgwas-pcg-preconditioner must be jacobi|lowrank");
            }
            cfg.ssgwas_pcg_preconditioner = v;
        }
        else if (arg == "--ssgwas-pcg-preconditioner-rank" && i + 1 < argc) {
            int v = stoi(argv[++i]);
            if (v <= 0) throw std::runtime_error("--ssgwas-pcg-preconditioner-rank must be > 0");
            cfg.ssgwas_pcg_preconditioner_rank = v;
        }
        else if (arg == "--allow-large-exact-ssgwas-pvalue") { cfg.allow_large_exact_ssgwas_pvalue = true; }
        else if (arg == "--ssgwas-vce" && i + 1 < argc) { cfg.ssgwas_vce = argv[++i]; }
        else if (arg == "--write-ebv") { cfg.write_ebv = true; }
        else if (arg == "--threshold-model") { cfg.threshold_model = true; cfg.run_gwas = true; }
        else if (arg == "--threshold-dichotomize") { cfg.threshold_dichotomize = true; }
        else if (arg == "--link" && i + 1 < argc) { cfg.link_type = argv[++i]; }
        else if (arg == "--count-model") { cfg.count_model = true; cfg.run_gwas = true; }
        else if (arg == "--family" && i + 1 < argc) { cfg.family_type = argv[++i]; }
        else if (arg == "--offset" && i + 1 < argc) { cfg.offset_file = argv[++i]; }
        else if (arg == "--binomial-pheno" && i + 1 < argc) { cfg.binomial_file = argv[++i]; }
        else if (arg == "--ordinal-model") { cfg.ordinal_model = true; cfg.run_gwas = true; }
        else if (arg == "--ordinal-pheno" && i + 1 < argc) { cfg.ordinal_file = argv[++i]; }
        else if (arg == "--ordinal-mode" && i + 1 < argc) { cfg.ordinal_mode = argv[++i]; }
        else if (arg == "--ordinal-stacked-engine" && i + 1 < argc) { cfg.ordinal_stacked_engine = argv[++i]; }
        else if (arg == "--ordinal-vc-method" && i + 1 < argc) { cfg.ordinal_vc_method = argv[++i]; }
        else if (arg == "--ordinal-sigma2-g" && i + 1 < argc) { cfg.ordinal_sigma2_g = stod(argv[++i]); }
        else if (arg == "--ordinal-vc-grid-min" && i + 1 < argc) { cfg.ordinal_vc_grid_min = stod(argv[++i]); }
        else if (arg == "--ordinal-vc-grid-max" && i + 1 < argc) { cfg.ordinal_vc_grid_max = stod(argv[++i]); }
        else if (arg == "--ordinal-vc-grid-points" && i + 1 < argc) { cfg.ordinal_vc_grid_points = stoi(argv[++i]); }
        else if (arg == "--ordinal-allow-fallback") { cfg.ordinal_allow_fallback = true; }
        else if (arg == "--ordinal-test" && i + 1 < argc) { cfg.ordinal_test = argv[++i]; }
        else if (arg == "--psr-null-only") { cfg.psr_null_only = true; }
        else if (arg == "--psr-maxit" && i + 1 < argc) { cfg.psr_maxit = stoi(argv[++i]); }
        else if (arg == "--psr-tol" && i + 1 < argc) { cfg.psr_tol = stod(argv[++i]); }
        else if (arg == "--psr-mu-eps" && i + 1 < argc) { cfg.psr_mu_eps = stod(argv[++i]); }
        else if (arg == "--psr-min-weight" && i + 1 < argc) { cfg.psr_min_weight = stod(argv[++i]); }
        else if (arg == "--psr-max-weight" && i + 1 < argc) { cfg.psr_max_weight = stod(argv[++i]); }
        else if (arg == "--survival-model") { cfg.survival_model = true; cfg.run_gwas = true; }
        else if (arg == "--gxe-model") { cfg.gxe_model = true; cfg.run_gwas = true; }
        else if (arg == "--env-covar" && i + 1 < argc) { cfg.env_covar_col = stoi(argv[++i]); }
        else if (arg == "--covar" && i + 1 < argc) { cfg.covar_file = argv[++i]; }
        else if (arg == "--qcovar" && i + 1 < argc) { cfg.qcovar_file = argv[++i]; }
        else if (arg == "--make-grm") { cfg.make_grm = true; }
        else if (arg == "-gk") {
            cfg.make_grm = true;
            if (i + 1 < argc && argv[i + 1][0] != '-') cfg.grm_alg = stoi(argv[++i]);
        }
        else if (arg == "--gs" || arg == "-gs") {
            cfg.gemma_s_mode = 1;
            if (i + 1 < argc && argv[i + 1][0] != '-') cfg.gemma_s_mode = stoi(argv[++i]);
        }
        else if (arg == "--gq" || arg == "-gq") {
            cfg.gemma_q_mode = 1;
            if (i + 1 < argc && argv[i + 1][0] != '-') cfg.gemma_q_mode = stoi(argv[++i]);
        }
        else if (arg == "--gemma-eigen" || arg == "-eigen") { cfg.gemma_eigen = true; }
        else if (arg == "--calccor" || arg == "--gemma-calccor" || arg == "-calccor" || arg == "-gw") {
            cfg.gemma_calccor = true;
            if (i + 1 < argc && argv[i + 1][0] != '-') ++i;
        }
        else if ((arg == "--windowbp" || arg == "-windowbp") && i + 1 < argc) { cfg.calccor_window_bp = stoll(argv[++i]); }
        else if ((arg == "--windowcm" || arg == "-windowcm") && i + 1 < argc) { cfg.calccor_window_cm = stod(argv[++i]); }
        else if ((arg == "--windowns" || arg == "-windowns") && i + 1 < argc) { cfg.calccor_window_ns = stoi(argv[++i]); }
        else if ((arg == "--vc" || arg == "-vc") && i + 1 < argc) { cfg.vc_mode = stoi(argv[++i]); }
        else if (arg == "--reml" || arg == "--vce") { cfg.vc_mode = 2; cfg.standalone_reml = true; }
        else if (arg == "--reml-method" && i + 1 < argc) { cfg.reml_method = argv[++i]; }
        else if (arg == "--ci" || arg == "-ci") {
            cfg.ci_mode = 1;
            if (i + 1 < argc && argv[i + 1][0] != '-') cfg.ci_mode = stoi(argv[++i]);
        }
        else if ((arg == "--beta" || arg == "-beta") && i + 1 < argc) { cfg.beta_file = argv[++i]; }
        else if ((arg == "--cor" || arg == "-cor") && i + 1 < argc) { cfg.cor_file = argv[++i]; }
        else if ((arg == "--cat" || arg == "-cat") && i + 1 < argc) { cfg.category_file = argv[++i]; }
        else if ((arg == "--mcat" || arg == "-mcat") && i + 1 < argc) { cfg.multi_category_file = argv[++i]; }
        else if ((arg == "--wcat" || arg == "-wcat") && i + 1 < argc) { cfg.category_weight_file = argv[++i]; }
        else if ((arg == "--wsnp" || arg == "-wsnp") && i + 1 < argc) { cfg.snp_weight_summary_file = argv[++i]; }
        else if ((arg == "--study" || arg == "-study") && i + 1 < argc) { cfg.study_prefix = argv[++i]; }
        else if ((arg == "--ref" || arg == "-ref") && i + 1 < argc) { cfg.ref_prefix = argv[++i]; }
        else if ((arg == "--mstudy" || arg == "-mstudy") && i + 1 < argc) { cfg.multi_study_file = argv[++i]; }
        else if ((arg == "--mref" || arg == "-mref") && i + 1 < argc) { cfg.multi_ref_file = argv[++i]; }
        else if ((arg == "--summary-blocks" || arg == "-nblock") && i + 1 < argc) { cfg.summary_block_count = stoi(argv[++i]); }
        else if (arg == "--make-grm-alg" && i + 1 < argc) { cfg.grm_alg = stoi(argv[++i]); }
        else if (arg == "--make-grm-sparse") { cfg.make_grm_sparse = true; }
        else if (arg == "--sparse-cutoff" && i + 1 < argc) { cfg.sparse_cutoff = stod(argv[++i]); }
        else if (arg == "--gwas") { cfg.run_gwas = true; }
        else if (arg == "--scan-tier" && i + 1 < argc) {
            cfg.scan_tier = argv[++i];
            if (cfg.scan_tier != "approx" && cfg.scan_tier != "semi" && cfg.scan_tier != "exact") {
                throw std::runtime_error("Unknown --scan-tier: " + cfg.scan_tier + " (expected approx, semi, or exact)");
            }
        }
        else if (arg == "--fast-lmm") { cfg.fast_lmm = true; }
        else if (arg == "--fastlmm") { cfg.fastlmm = true; cfg.run_gwas = true; }
else if (arg == "--fastlmm-trace-threshold" && i + 1 < argc) { cfg.fastlmm_trace_threshold = stod(argv[++i]); }
    else if (arg == "--fastlmm-min-rank" && i + 1 < argc) { cfg.fastlmm_min_rank = stoi(argv[++i]); }
    else if (arg == "--fastlmm-max-rank" && i + 1 < argc) { cfg.fastlmm_max_rank = stoi(argv[++i]); }
    else if (arg == "--fastlmm-matrix-free") { cfg.fastlmm_matrix_free = true; }
    else if (arg == "--fastlmm-no-fallback") { cfg.fastlmm_no_fallback = true; }
        else if (arg == "--lmm-low-rank" && i + 1 < argc) { cfg.lmm_low_rank = stoi(argv[++i]); cfg.run_gwas = true; }
        else if (arg == "--lmm-low-rank-oversample" && i + 1 < argc) { cfg.lmm_low_rank_oversample = stoi(argv[++i]); }
        else if (arg == "--lmm-low-rank-seed" && i + 1 < argc) { cfg.lmm_low_rank_seed = stoll(argv[++i]); }
        else if (arg == "--lmm-low-rank-error-tol" && i + 1 < argc) { cfg.lmm_low_rank_error_tol = stod(argv[++i]); }
        else if ((arg == "--lmm-gemma-mode" || arg == "-lmm") && i + 1 < argc) {
            int mode = stoi(argv[++i]);
            if (mode == 5) {
                cfg.vc_mode = 2;
                cfg.run_gwas = false;
            } else if (mode == 9) {
                cfg.lmm_gemma_mode = 2;
                cfg.run_gwas = true;
            } else {
                cfg.lmm_gemma_mode = mode;
                cfg.run_gwas = true;
            }
        }
        else if (arg == "--lmm-gemma-refit") { cfg.lmm_gemma_refit = true; cfg.run_gwas = true; if (cfg.lmm_gemma_mode == 0) cfg.lmm_gemma_mode = 4; }
        else if (arg == "--lm") { cfg.linear_model = true; cfg.run_gwas = true; }
        else if (arg == "-lm") {
            cfg.linear_model = true;
            cfg.run_gwas = true;
            if (i + 1 < argc && argv[i + 1][0] != '-') cfg.linear_model_mode = stoi(argv[++i]);
        }
        else if (arg == "--lm-mode" && i + 1 < argc) { cfg.linear_model_mode = stoi(argv[++i]); }
        else if (arg == "--model") {
            if (i + 1 >= argc || argv[i + 1][0] == '-')
                throw std::runtime_error("--model requires a value; see --help models.");
            if (!selected_cli_model.empty()) throw std::runtime_error("Specify --model only once.");
            selected_cli_model = argv[++i];
            if (selected_cli_model == "lmm" || selected_cli_model == "fixed-lmm") cfg.eb_lmm = false;
            else if (selected_cli_model == "lmm-eb" || selected_cli_model == "eb-lmm") cfg.eb_lmm = true;
        }
        else if (arg == "--gwas-method" && i + 1 < argc) {
            std::string model = argv[++i];
            if (model == "lmm" || model == "fixed-lmm") cfg.eb_lmm = false;
            else if (model == "lmm-eb" || model == "eb-lmm") cfg.eb_lmm = true;
            else throw std::runtime_error("Unknown GWAS model/method: " + model + " (expected lmm or lmm-eb)");
        }
        else if (arg == "--test" && i + 1 < argc) {
            std::string test = argv[++i];
            if (test == "eb-wald") cfg.eb_lmm = true;
            else if (test == "wald" || test == "score" || test == "lrt") { /* accepted for CLI compatibility */ }
            else throw std::runtime_error("Unknown GWAS test: " + test + " (expected wald, score, lrt, or eb-wald)");
        }
        else if (arg == "--adjust" && i + 1 < argc) { ++i; }
        else if (arg == "--repeatability-mixed-precision") { cfg.repeatability_mixed_precision = true; }
        else if (arg == "--repeatability-grm-cache" && i + 1 < argc) { cfg.repeatability_grm_cache_prefix = argv[++i]; }
        else if (arg == "--repeatability-weight" && i + 1 < argc) { cfg.repeatability_weight_file = argv[++i]; }
        else if (arg == "--repeatability-refit-validation") { cfg.repeatability_refit_validation = true; }
        else if (arg == "--repeatability-refit-limit" && i + 1 < argc) { cfg.repeatability_refit_limit = stoi(argv[++i]); }
        else if (arg == "--repeatability-mode" && i + 1 < argc) { cfg.repeatability_mode = argv[++i]; }
        else if (arg == "--vce-mode" && i + 1 < argc) { cfg.vce_refine_mode = argv[++i]; }
        else if (arg == "--threads" && i + 1 < argc) { cfg.num_threads = stoi(argv[++i]); }
        else if (arg == "--omp" && i + 1 < argc) { cfg.num_threads = stoi(argv[++i]); }
        else if ((arg == "--memory-gib" || arg == "--memory-gb") && i + 1 < argc) {
            auto parsed = parse_memory_gib(argv[++i]);
            if (!parsed.ok) throw std::runtime_error("Invalid --memory-gib: " + parsed.error);
            cfg.memory_budget_specified = true;
            cfg.memory_budget_bytes = parsed.budget.bytes;
            cfg.memory_budget_source = to_string(parsed.budget.source);
        }
        else if (arg == "--spa") { cfg.use_spa = true; }
        else if (arg == "--skat") { cfg.run_skat = true; }
        else if (arg == "--burden") { cfg.run_burden = true; }
        else if (arg == "--set-file" && i + 1 < argc) { cfg.set_file = argv[++i]; }
        else if (arg == "--biobank-mode") { cfg.biobank_mode = true; }
        else if (arg == "--block-size" && i + 1 < argc) { cfg.block_size = stoi(argv[++i]); }
        else if (arg == "--regenie-model-snps" && i + 1 < argc) { cfg.regenie_model_snps = argv[++i]; }
        else if (arg == "--write-pred" && i + 1 < argc) { cfg.regenie_write_pred = argv[++i]; }
        else if (arg == "--use-pred" && i + 1 < argc) { cfg.regenie_use_pred = argv[++i]; }
        else if (arg == "--regenie-step" && i + 1 < argc) { cfg.regenie_step = stoi(argv[++i]); }
        else if ((arg == "--predict-effect" || arg == "-epm") && i + 1 < argc) { cfg.predict_effect_file = argv[++i]; }
        else if (arg == "-predict") {
            if (i + 1 < argc && argv[i + 1][0] != '-') ++i; // GEMMA mode number; CS infers prediction type from effect file.
        }
        else if (arg == "--predict-score" && i + 1 < argc) { cfg.predict_score_file = argv[++i]; }
        else if (arg == "--predict-snp-col" && i + 1 < argc) { cfg.predict_snp_col = argv[++i]; }
        else if (arg == "--predict-effect-col" && i + 1 < argc) { cfg.predict_effect_col = argv[++i]; }
        else if ((arg == "--bslmm" || arg == "-bslmm") && i + 1 < argc) { cfg.bslmm_mode = stoi(argv[++i]); }
        else if (arg == "--bslmm-mcmc") { cfg.bslmm_mode = 1; cfg.bslmm_mcmc = true; }
        else if (arg == "--bvsr") { cfg.bslmm_mode = 1; cfg.bvsr = true; }
        else if (arg == "--bslmm-top" && i + 1 < argc) { cfg.bslmm_top = stoi(argv[++i]); }
        else if (arg == "--bslmm-h2" && i + 1 < argc) { cfg.bslmm_h2 = stod(argv[++i]); }
        else if ((arg == "--bf" || arg == "-bf") && i + 1 < argc) { cfg.bslmm_bf_file = argv[++i]; }
        else if ((arg == "--hyp" || arg == "-hyp") && i + 1 < argc) { cfg.bslmm_hyp_file = argv[++i]; }
        else if (arg == "--bslmm-hgrid" && i + 1 < argc) { cfg.bslmm_hgrid = stoi(argv[++i]); }
        else if (arg == "--bslmm-rho-grid" && i + 1 < argc) { cfg.bslmm_rho_grid = stoi(argv[++i]); }
        else if (arg == "-hmin" && i + 1 < argc) { cfg.bslmm_h_min = stod(argv[++i]); }
        else if (arg == "-hmax" && i + 1 < argc) { cfg.bslmm_h_max = stod(argv[++i]); }
        else if (arg == "-rmin" && i + 1 < argc) { cfg.bslmm_rho_min = stod(argv[++i]); }
        else if (arg == "-rmax" && i + 1 < argc) { cfg.bslmm_rho_max = stod(argv[++i]); }
        else if (arg == "-pmin" && i + 1 < argc) { cfg.bslmm_p_min = stod(argv[++i]); cfg.bslmm_p_min_set = true; }
        else if (arg == "-pmax" && i + 1 < argc) { cfg.bslmm_p_max = stod(argv[++i]); cfg.bslmm_p_max_set = true; }
        else if (arg == "-gmean" && i + 1 < argc) { cfg.bslmm_gamma_mean = stod(argv[++i]); }
        else if (arg == "-hscale" && i + 1 < argc) { cfg.bslmm_h_scale = stod(argv[++i]); }
        else if (arg == "-rscale" && i + 1 < argc) { cfg.bslmm_rho_scale = stod(argv[++i]); }
        else if (arg == "-pscale" && i + 1 < argc) { cfg.bslmm_p_scale = stod(argv[++i]); }
        else if (arg == "-mh" && i + 1 < argc) { cfg.bslmm_mh_steps = stoi(argv[++i]); }
        else if (arg == "-smin" && i + 1 < argc) { cfg.bslmm_s_min = stoi(argv[++i]); }
        else if (arg == "-smax" && i + 1 < argc) { cfg.bslmm_top = stoi(argv[++i]); }
        else if (arg == "-rpace" && i + 1 < argc) { cfg.bslmm_rpace = stoi(argv[++i]); }
        else if (arg == "-wpace" && i + 1 < argc) { cfg.bslmm_wpace = stoi(argv[++i]); }
        else if (arg == "-w" && i + 1 < argc) { cfg.bslmm_burnin = stoi(argv[++i]); }
        else if (arg == "-s" && i + 1 < argc) { cfg.bslmm_sampling = stoi(argv[++i]); }
        else if (arg == "-seed" && i + 1 < argc) { cfg.bslmm_seed = stoll(argv[++i]); }
        else if (arg == "--mcmc-burnin" && i + 1 < argc) { cfg.bslmm_burnin = stoi(argv[++i]); }
        else if (arg == "--mcmc-samples" && i + 1 < argc) { cfg.bslmm_sampling = stoi(argv[++i]); }
        else if (arg == "--mcmc-thin" && i + 1 < argc) { cfg.bslmm_thin = stoi(argv[++i]); }
        else if (arg == "--seed" && i + 1 < argc) { cfg.bslmm_seed = stoll(argv[++i]); }
        else if (arg == "--chr" && i + 1 < argc) { cfg.scan_chr = argv[++i]; }
        else if (arg == "--from-bp" && i + 1 < argc) { cfg.scan_from_bp = stoll(argv[++i]); }
        else if (arg == "--to-bp" && i + 1 < argc) { cfg.scan_to_bp = stoll(argv[++i]); }
        else if (arg == "--no-covar-proj") { cfg.use_covariate_projection = false; }
        else if (arg == "--force-slow") { cfg.force_slow = true; }
        else if (arg == "--no-loco") { cfg.no_loco = true; }
        else if (arg == "--maf" && i + 1 < argc) { cfg.min_maf = stod(argv[++i]); }
        else if (arg == "--hwe" && i + 1 < argc) { cfg.min_hwe = stod(argv[++i]); }
        else if (arg == "--bgen-dom-max-missing" && i + 1 < argc) {
            cfg.bgen_dom_max_missing = stod(argv[++i]); bgen_dom_qc_requested = true;
        }
        else if (arg == "--bgen-dom-min-info" && i + 1 < argc) {
            cfg.bgen_dom_min_info = stod(argv[++i]); bgen_dom_qc_requested = true;
        }
        else if (arg == "--bgen-dom-entropy-warning" && i + 1 < argc) {
            cfg.bgen_dom_entropy_warning = stod(argv[++i]); bgen_dom_qc_requested = true;
        }
        else if (arg == "--extract" && i + 1 < argc) { cfg.extract_file = argv[++i]; }
        else if (arg == "--exclude" && i + 1 < argc) { cfg.exclude_file = argv[++i]; }
        else if (arg == "--lmm-null-solver" && i + 1 < argc) { cfg.lmm_null_solver = argv[++i]; }
        else if (arg == "--snp-weight" && i + 1 < argc) { cfg.snp_weight_file = argv[++i]; }
        else if (arg == "--pop-class" && i + 1 < argc) { cfg.pop_class_file = argv[++i]; }
        else if (arg == "--reml-maxit" && i + 1 < argc) { cfg.reml_maxit = stoi(argv[++i]); }
        else if (arg == "--reml-tol" && i + 1 < argc) { cfg.reml_tol = stod(argv[++i]); }
        else if (arg == "--reml-vg" && i + 1 < argc) { cfg.reml_vg = stod(argv[++i]); }
        else if (arg == "--reml-vpe" && i + 1 < argc) { cfg.reml_vpe = stod(argv[++i]); }
        else if (arg == "--reml-vc" && i + 1 < argc) { cfg.reml_vc = stod(argv[++i]); }
        else if (arg == "--reml-ve" && i + 1 < argc) { cfg.reml_ve = stod(argv[++i]); }
        else if (arg == "--help" || arg == "-h") { print_usage(argv[0]); return 0; }
        else if (arg == "--version") { cout << "CosmicScanner v" << COSMICSCANNER_VERSION << "\n"; return 0; }
        else { throw std::runtime_error("Unknown option: " + arg); }
    }

    applyCliModel(cfg, selected_cli_model, invocation);

    if (cfg.skat_o || cfg.acat_v) {
        if (int(cfg.skat_o) + int(cfg.acat_v) + int(cfg.run_skat) + int(cfg.run_burden) != 1)
            throw std::runtime_error("Choose one rare-variant test per run");
        cfg.rare_method = cfg.skat_o ? "skat-o" : "acat-v";
        cfg.run_skat = true;
    }
    if (!std::isfinite(cfg.rare_max_maf) || cfg.rare_max_maf <= 0.0 || cfg.rare_max_maf > 0.5)
        throw std::runtime_error("--rare-max-maf must be in (0, 0.5]");
    if (cfg.rare_min_mac < 1) throw std::runtime_error("--rare-min-mac must be positive");
    if (cfg.snp_model != "additive" && cfg.snp_model != "dominant" &&
        cfg.snp_model != "recessive" && cfg.snp_model != "heterozygote")
        throw std::runtime_error("Invalid --snp-model");
    if (cfg.snp_model != "additive" && !cfg.bgen_file.empty())
        throw std::runtime_error("Non-additive models require hardcalls; BGEN expected dosages are unsupported");
    if (cfg.gwas_model != "additive" && cfg.gwas_model != "additive_dominance" && cfg.gwas_model != "dominance")
        throw std::runtime_error("Invalid --gwas-model (expected additive, additive_dominance, or dominance)");
    if (!std::isfinite(cfg.bgen_dom_max_missing) || cfg.bgen_dom_max_missing < 0.0 ||
        cfg.bgen_dom_max_missing > 1.0)
        throw std::runtime_error("--bgen-dom-max-missing must be in [0,1]");
    if (!std::isfinite(cfg.bgen_dom_min_info) || cfg.bgen_dom_min_info < 0.0)
        throw std::runtime_error("--bgen-dom-min-info must be >= 0");
    if (!std::isfinite(cfg.bgen_dom_entropy_warning) || cfg.bgen_dom_entropy_warning < 0.0 ||
        cfg.bgen_dom_entropy_warning > 1.0)
        throw std::runtime_error("--bgen-dom-entropy-warning must be in [0,1]");
    if (bgen_dom_qc_requested && (cfg.bgen_file.empty() ||
        (cfg.gwas_model != "dominance" && !(cfg.make_grm && cfg.grm_alg == 4))))
        throw std::runtime_error("--bgen-dom-* QC options require BGEN dominance or BGEN --make-grm-alg 4");
    if (cfg.gwas_model == "additive_dominance") {
        if (!cfg.bgen_file.empty())
            throw std::runtime_error("--gwas-model additive_dominance is not supported for BGEN input: "
                                     "the BGEN dosage-only backend does not expose P(AB) or hardcall "
                                     "genotype classes, so it cannot represent the production "
                                     "dominance model. Convert to PLINK BED or PGEN hardcalls, or use "
                                     "the default additive model.");
        if (cfg.bfile_prefix.empty() && cfg.pfile_prefix.empty())
            throw std::runtime_error("--gwas-model additive_dominance requires PLINK BED (--bfile) or "
                                     "PGEN (--pfile) hardcall input");
        if (cfg.snp_model != "additive")
            throw std::runtime_error("--gwas-model additive_dominance currently requires additive --snp-model");
        if (cfg.multivariate || cfg.rrm_model || cfg.repeatability_model || cfg.biobank_mode)
            throw std::runtime_error("--gwas-model additive_dominance is only supported for standard single-trait LMM GWAS");
        if (cfg.run_skat || cfg.run_burden || cfg.rare_method.size())
            throw std::runtime_error("--gwas-model additive_dominance is not supported for rare-set tests");
        if (!cfg.use_covariate_projection)
            throw std::runtime_error("--gwas-model additive_dominance requires covariate projection "
                                     "(it is not compatible with --no-covar-proj)");
        if (cfg.out_gz)
            throw std::runtime_error("--gwas-model additive_dominance does not yet support --out-gz output");
        if (cfg.out_ssf)
            throw std::runtime_error("--gwas-model additive_dominance does not yet support --out-ssf (GWAS-SSF) output");
        if (cfg.manhattan_tsv || cfg.qq_tsv)
            throw std::runtime_error("--gwas-model additive_dominance does not yet support --manhattan-tsv/--qq-tsv output");
        if (cfg.fastlmm || cfg.lmm_low_rank > 0)
            throw std::runtime_error("--gwas-model additive_dominance always uses the exact fixed-null scan "
                                     "and does not support --fastlmm/--lmm-low-rank");
        cfg.run_gwas = true;
    }
    if (cfg.gwas_model == "dominance") {
        if (cfg.bfile_prefix.empty() && cfg.pfile_prefix.empty() && cfg.bgen_file.empty())
            throw std::runtime_error("--gwas-model dominance requires PLINK BED (--bfile), "
                                     "PGEN (--pfile) hardcall, or BGEN (--bgen) probability input");
        if (cfg.snp_model != "additive")
            throw std::runtime_error("--gwas-model dominance currently requires additive --snp-model");
        if (cfg.saige_style)
            throw std::runtime_error("--gwas-model dominance does not support --saige-style");
        if (cfg.multivariate || cfg.rrm_model || cfg.repeatability_model || cfg.biobank_mode)
            throw std::runtime_error("--gwas-model dominance is only supported for standard single-trait LMM GWAS");
        if (cfg.run_skat || cfg.run_burden || cfg.rare_method.size())
            throw std::runtime_error("--gwas-model dominance is not supported for rare-set tests");
        if (!cfg.use_covariate_projection)
            throw std::runtime_error("--gwas-model dominance requires covariate projection "
                                     "(it is not compatible with --no-covar-proj)");
        if (cfg.out_gz)
            throw std::runtime_error("--gwas-model dominance does not yet support --out-gz output");
        if (cfg.out_ssf)
            throw std::runtime_error("--gwas-model dominance does not yet support --out-ssf (GWAS-SSF) output");
        if (cfg.manhattan_tsv || cfg.qq_tsv)
            throw std::runtime_error("--gwas-model dominance does not yet support --manhattan-tsv/--qq-tsv output");
        if (cfg.fastlmm || cfg.lmm_low_rank > 0)
            throw std::runtime_error("--gwas-model dominance always uses the exact fixed-null scan "
                                     "and does not support --fastlmm/--lmm-low-rank");
        if (cfg.lmm_gemma_refit || cfg.lmm_gemma_mode != 0)
            throw std::runtime_error("--gwas-model dominance does not support --lmm-gemma-mode/--lmm-gemma-refit");
        if (!cfg.bgen_file.empty() && !cfg.pop_class_file.empty())
            throw std::runtime_error("BGEN probability dominance does not support --pop-class");
        if (!cfg.bgen_file.empty() && !cfg.snp_weight_file.empty())
            throw std::runtime_error("BGEN probability dominance does not support --snp-weight");
        cfg.run_gwas = true;
    }
    if (cfg.fastlmm && cfg.lmm_low_rank > 0)
        throw std::runtime_error("Choose either --fastlmm or --lmm-low-rank, not both");
    if (!(cfg.fastlmm_trace_threshold > 0.0 && cfg.fastlmm_trace_threshold <= 1.0) ||
        cfg.fastlmm_min_rank < 1 || cfg.fastlmm_max_rank < cfg.fastlmm_min_rank)
        throw std::runtime_error("FastLMM requires trace threshold in (0,1] and 1 <= min rank <= max rank");
    if (cfg.fastlmm_matrix_free && cfg.fastlmm_trace_threshold > 0.0 && !cfg.fastlmm_no_fallback)
        throw std::runtime_error("--fastlmm-matrix-free requires --fastlmm-no-fallback when trace threshold > 0 (use fixed rank mode)");
    if (cfg.lmm_low_rank < 0 || cfg.lmm_low_rank_oversample < 0 ||
        !(cfg.lmm_low_rank_error_tol > 0.0))
        throw std::runtime_error("Low-rank LMM requires rank/oversampling >= 0 and error tolerance > 0");
    if (!cfg.logistic_test.empty() && cfg.logistic_test != "score" && cfg.logistic_test != "wald")
        throw std::runtime_error("--logistic must be score or wald");
    const int multilocus_models = int(cfg.farmcpu) + int(cfg.blink) + int(cfg.firth_logistic) +
                                  int(!cfg.logistic_test.empty()) + int(cfg.logistic_glmm) +
                                  int(!cfg.condition_snp_file.empty());
    if (multilocus_models > 1)
        throw std::runtime_error("Choose one of --farmcpu, --blink, --firth-logistic, --logistic, --logistic-glmm, or --condition-snp");
    if (cfg.saige_style && !cfg.biobank_mode)
        throw std::runtime_error("--saige-style requires --biobank-mode");
    if (cfg.geno_orientation != "variants-by-samples" && cfg.geno_orientation != "samples-by-variants")
        throw std::runtime_error("--geno-orientation must be variants-by-samples or samples-by-variants");
    if (!cfg.geno_var_file.empty() && cfg.geno_matrix_file.empty())
        throw std::runtime_error("--geno-var requires --geno-matrix");
    if (!cfg.geno_sample_file.empty() && cfg.geno_matrix_file.empty())
        throw std::runtime_error("--geno-sample requires --geno-matrix");
    if (cfg.npy_orientation != "variants-by-samples" && cfg.npy_orientation != "samples-by-variants")
        throw std::runtime_error("--npy-orientation must be variants-by-samples or samples-by-variants");
    if (!cfg.npy_var_file.empty() && cfg.npy_file.empty())
        throw std::runtime_error("--npy-var requires --npy");
    if (!cfg.npy_sample_file.empty() && cfg.npy_file.empty())
        throw std::runtime_error("--npy-sample requires --npy");
    // --ped/--map act as a PLINK PED/MAP genotype source in standalone form;
    // in animal/pedigree models (--single-step, --animal-model, --make-agh)
    // --ped is instead the model pedigree and must not count as a genotype source.
    const bool ped_is_pedigree = cfg.single_step_model || cfg.animal_model || cfg.make_agh;
    if (!ped_is_pedigree && !cfg.ped_file.empty() && cfg.map_file.empty())
        throw std::runtime_error("--ped requires --map");
    if (!ped_is_pedigree && !cfg.map_file.empty() && cfg.ped_file.empty())
        throw std::runtime_error("--map requires --ped");
    const int genotype_sources = int(!cfg.bfile_prefix.empty()) + int(!cfg.bgen_file.empty()) +
        int(!cfg.pfile_prefix.empty() || !cfg.pgen_file.empty()) + int(!cfg.vcf_file.empty()) +
        int(!cfg.hmp_file.empty()) + int(!cfg.geno_matrix_file.empty()) + int(!cfg.npy_file.empty()) +
        (ped_is_pedigree ? 0 : int(!cfg.ped_file.empty()));
    if (genotype_sources > 1)
        throw std::runtime_error("Specify exactly one genotype source");
    if (!cfg.convert_from.empty() && cfg.convert_to.empty())
        throw std::runtime_error(
            "--convert " + cfg.convert_from + " requires --to <text|npy|ped|vcf>; without --to the "
            "request is a no-op and no conversion would be written");
    if (!cfg.vcf_file.empty()) {
        // R33-FIX (B2/B3): the multivariate and repeatability engines have no VCF scan path.
        // MultivariateModel only implements MvLmmEngine::scan() for Plink/Pgen/Bgen, and
        // RepeatabilityModel's STCG branch attaches the (empty under --vcf) PlinkReader as the
        // genotype-matrix operator. Routing VCF into either model would either fail deep in the
        // GRM builder or silently select the wrong genotype source, so the combination is
        // explicitly unsupported instead of being wired up optimistically.
        if (cfg.multivariate)
            throw std::runtime_error(
                "Multivariate GWAS does not support VCF input: supported genotype formats are "
                "PLINK BED (--bfile), PLINK PGEN (--pfile) and BGEN (--bgen). Convert the VCF to "
                "one of those formats before running --multivariate.");
        if (cfg.repeatability_model)
            throw std::runtime_error(
                "Repeatability GWAS does not support VCF input: supported genotype formats are "
                "PLINK BED (--bfile), PLINK PGEN (--pfile) and BGEN (--bgen). Convert the VCF to "
                "one of those formats before running --repeatability-model.");
    }
    if (cfg.conditional_joint && cfg.condition_snp_file.empty())
        throw std::runtime_error("--conditional-joint requires --condition-snp <file>");
    if (cfg.farmcpu_max_iterations < 1 || cfg.farmcpu_max_iterations > 100)
        throw std::runtime_error("--farmcpu-max-iterations must be between 1 and 100");
    if (cfg.farmcpu_pseudo_qtn < 1 || cfg.farmcpu_pseudo_qtn > 1000)
        throw std::runtime_error("--farmcpu-pseudo-qtn must be between 1 and 1000");
    if (multilocus_models && cfg.snp_model != "additive")
        throw std::runtime_error("FarmCPU, BLINK, logistic, and conditional GWAS currently require additive genotypes");

    // COJO validation
    const int cojo_models = int(cfg.cojo) + int(cfg.cojo_slct) + int(cfg.cojo_joint);
    if (cojo_models > 1)
        throw std::runtime_error("Choose one of --cojo, --cojo-slct, or --cojo-joint");
    if (cojo_models > 0) {
        if (cfg.cojo_summary.empty())
            throw std::runtime_error("--cojo/--cojo-slct/--cojo-joint requires --cojo-summary <file>");
        if (cfg.cojo_ld_ref.empty() && genotype_sources == 0)
            throw std::runtime_error("COJO requires either a genotype source or --cojo-ld-ref <file>");
        if (cfg.cojo_p <= 0.0 || cfg.cojo_p >= 1.0)
            throw std::runtime_error("--cojo-p must be in (0, 1)");
        if (cfg.cojo_r2 <= 0.0 || cfg.cojo_r2 >= 1.0)
            throw std::runtime_error("--cojo-r2 must be in (0, 1)");
        if (cfg.cojo_window <= 0)
            throw std::runtime_error("--cojo-wind must be positive");
    }
    if (cfg.fine_map) {
        if (cfg.finemap_summary.empty()) throw std::runtime_error("--fine-map requires --finemap-summary <file>");
        if (cfg.finemap_region.empty()) throw std::runtime_error("--fine-map requires --finemap-region chr:start-end");
        if (genotype_sources != 1 || (!cfg.hmp_file.empty() || !cfg.geno_matrix_file.empty() || !cfg.npy_file.empty()))
            throw std::runtime_error("--fine-map requires exactly one real LD source: --bfile, --pfile, --bgen, or --vcf");
        if (!std::isfinite(cfg.credible_level) || cfg.credible_level <= 0.0 || cfg.credible_level >= 1.0)
            throw std::runtime_error("--credible-level must be in (0, 1)");
        if (!std::isfinite(cfg.finemap_prior_var) || cfg.finemap_prior_var <= 0.0)
            throw std::runtime_error("--finemap-prior-var must be finite and positive");
        if (cfg.finemap_max_signals < 1) throw std::runtime_error("--finemap-max-signals must be positive");
        if (!std::isfinite(cfg.finemap_signal_p) || cfg.finemap_signal_p <= 0.0 || cfg.finemap_signal_p >= 1.0)
            throw std::runtime_error("--finemap-signal-p must be in (0, 1)");
    }

    // Meta-analysis validation
    if (cfg.meta_analysis) {
        if (cfg.meta_studies_file.empty() && cfg.meta_study_files.empty())
            throw std::runtime_error("--meta-analysis requires --meta-studies <file>");
        if (cfg.meta_method != "fixed" && cfg.meta_method != "random-dl" && cfg.meta_method != "random-reml")
            throw std::runtime_error("--meta-method must be fixed, random-dl, or random-reml");
        if (cfg.meta_het_p <= 0.0 || cfg.meta_het_p >= 1.0)
            throw std::runtime_error("--meta-het-p must be in (0, 1)");
        if (cfg.meta_min_studies < 1)
            throw std::runtime_error("--meta-min-studies must be >= 1");
        if (cfg.meta_gc_lambda <= 0.0)
            throw std::runtime_error("--meta-gc-lambda must be positive");
        if (cfg.meta_out_format != "metal" && cfg.meta_out_format != "simple")
            throw std::runtime_error("--meta-out-format must be metal or simple");
    }
    if (cfg.grm_alg < 0 || cfg.grm_alg > 6) throw std::runtime_error("--make-grm-alg/-gk must be between 0 and 6");
    if (!std::isfinite(cfg.sparse_cutoff) || cfg.sparse_cutoff < 0.0) throw std::runtime_error("--sparse-cutoff must be finite and non-negative");
    if (!std::isfinite(cfg.h_tau) || cfg.h_tau < 0.0 || !std::isfinite(cfg.h_omega) || cfg.h_omega < 0.0) {
        throw std::runtime_error("--h2-tau and --h2-omega must be finite and non-negative");
    }
    if (!std::isfinite(cfg.blend_g) || cfg.blend_g < 0.0 || cfg.blend_g > 1.0) {
        throw std::runtime_error("--blend-g must be between 0 and 1");
    }
    if (!std::isfinite(cfg.ssgwas_apy_core_target_ratio) || cfg.ssgwas_apy_core_target_ratio < 0.0) {
        throw std::runtime_error("--ssgwas-apy-core-target-ratio must be finite and non-negative (0 = unset)");
    }
    if (cfg.write_agh && !cfg.single_step_model) throw std::runtime_error("--write-agh requires --single-step or --make-agh");
    if (cfg.animal_model && cfg.ped_file.empty()) {
        throw std::runtime_error("--animal-model requires a pedigree: please provide --ped <file> (ID SIRE DAM, no header)");
    }
    if (cfg.single_step_model) {
        const bool is_backsolve = (cfg.single_step_method == "backsolve");
        if (is_backsolve && cfg.make_agh) {
            throw std::runtime_error("--make-agh (matrix build only) is a Method-1 workflow; it cannot be combined with --single-step-method backsolve");
        }
        if (is_backsolve && cfg.ped_file.empty()) {
            throw std::runtime_error("--single-step-method backsolve requires a pedigree: please provide --ped <file>");
        }
        if (is_backsolve && cfg.ssgwas_window_snps <= 0) {
            throw std::runtime_error("--ssgwas-window-snps must be a positive integer");
        }
        if (is_backsolve && cfg.grm_alg != 1) {
            throw std::runtime_error(
                "--single-step-method backsolve requires the VanRaden 2p(1-p) GRM (default --grm-alg 1); "
                "other GRM algorithms do not have a canonical back-solving denominator");
        }
        if (is_backsolve && cfg.reml_vg > 0.0 && cfg.reml_ve > 0.0) {
            if (!cfg.ssgwas_vce.empty()) {
                throw std::runtime_error("Provide either fixed --reml-vg/--reml-ve or --ssgwas-vce, not both");
            }
        } else if (is_backsolve && cfg.ssgwas_vce.empty()) {
            throw std::runtime_error(
                "--single-step-method backsolve requires variance components: "
                "provide fixed --reml-vg <f> --reml-ve <f> (canonical/deterministic path) or --ssgwas-vce <mode>");
        }
    }

    if (!cfg.write_template_config.empty()) {
        try {
            write_template_config(cfg, cfg.write_template_config);
            cout << "Input config template written to " << cfg.write_template_config << "\n";
            return 0;
        } catch (const exception& e) {
            cerr << "Error: " << e.what() << endl;
            return 1;
        }
    }

    // COJO and Meta-analysis are standalone summary-stats analyses that don't require genotype input
bool is_standalone_summary = cfg.cojo || cfg.cojo_slct || cfg.cojo_joint || cfg.meta_analysis || cfg.fine_map;
if (!is_standalone_summary && cfg.bfile_prefix.empty() && cfg.bgen_file.empty() && cfg.pfile_prefix.empty() && cfg.pgen_file.empty() && cfg.vcf_file.empty() && cfg.hmp_file.empty() && cfg.geno_matrix_file.empty() && cfg.npy_file.empty() && cfg.ped_file.empty()) {
    std::cerr << "Error: Specify a genotype source (e.g. --bfile PREFIX). See --help input.\n";
    printCliHelp("overview", argv[0]);
    return 1;
}

    GwasApp app(cfg);
    return app.run();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
