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
    cout << "Usage: " << prog_name << " --bfile <prefix> [options]\n";
    cout << "\n=== INPUT/OUTPUT ===\n";
    cout << "  --bfile <prefix>       PLINK 1.x binary prefix (.bed/.bim/.fam)\n";
    cout << "  --bgen <file>          BGEN file (Layout 2, zlib compressed)\n";
    cout << "  --sample <file>        Sample file for BGEN (Oxford format)\n";
    cout << "  --pfile <prefix>       PLINK 2.0 prefix (.pgen/.pvar/.psam)\n";
    cout << "  --pgen <file>          PLINK 2.0 genotype file\n";
    cout << "  --pvar <file>          PLINK 2.0 variant file\n";
    cout << "  --psam <file>          PLINK 2.0 sample file\n";
    cout << "  --out/-o <prefix>      Output prefix (default: cosmic_gwas)\n";
    cout << "  --out-ssf              Output in EBI GWAS-SSF format\n";
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
    cout << "  --gwas                 Run GWAS scan (required to produce association output;\n";
    cout << "                          all models need this flag explicitly, except --rrm-unbalance)\n";
    cout << "  --no-loco              Disable leave-one-chromosome-out (run global)\n";
    cout << "  --scan-tier <tier>     Per-SNP VC handling: approx, semi, exact\n";
    cout << "                          approx=fixed null; semi=fixed ratio/profile scale;\n";
    cout << "                          exact=full per-SNP variance-component refit\n";
    cout << "  --fast-lmm             Legacy alias for fixed-null single-trait scan\n";
    cout << "  --lmm-gemma-mode/-lmm <n> GEMMA LMM mode: 1 Wald, 2 LRT, 3 Score, 4 all\n";
    cout << "  --lmm-gemma-refit      Per-SNP alternative REML refit for GEMMA LMM Wald/LRT\n";
    cout << "  --lm/-lm <n>           GEMMA-style ordinary linear model GWAS\n";
    cout << "  --lm-mode <n>          GEMMA LM mode: 1 Wald, 2 LRT, 3 Score, 4 all\n";
    cout << "  --model lmm-eb         Empirical Bayes random-SNP LMM GWAS (experimental; report standard LMM too)\n";
    cout << "  --gwas-method eb-lmm   Alias for --model lmm-eb\n";
    cout << "  --test eb-wald         Alias enabling EB-LMM Wald test\n";
    cout << "  --adjust eb-mbonf      Experimental threshold in .eb.summary.txt; do not report alone\n";
    cout << "  --spa                  Use Saddle Point Approximation for binary traits\n";
    cout << "  --threads <n>          Number of OpenMP threads (default: max)\n";
    cout << "  --chr <chr>            Restrict scan to one chromosome\n";
    cout << "  --from-bp <pos>        Restrict scan start position\n";
    cout << "  --to-bp <pos>          Restrict scan end position\n";
    cout << "\n=== VCE OPTIONS (all models) ===\n";
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
    cout << "\n--- Random Regression Model (RRM) ---\n";
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
    cout << "\n--- Maternal Effects Model ---\n";
    cout << "  --maternal-model       Direct + maternal + PE multi-component model\n";
    cout << "  --pheno format: FID IID DAM_ID PHENO  (4+ columns)\n";
    cout << "  --maternal-id <col>    Maternal ID column (default: 3, 1-indexed)\n";
    cout << "  VCE: direct additive (ad), maternal additive (am),\n";
    cout << "       direct-maternal covariance (adm), PE (c), residual (e)\n";
    cout << "  Scan uses effective VCs: Vg=ad+am+(c) with LOCO GRM.\n";
    cout << "\n--- Single-Step GWAS (ssGWAS) Model ---\n";
    cout << "  --single-step          ssGWAS with A22-corrected pedigree + genomic H matrix\n";
    cout << "  --make-agh             Build A-inverse, G, and H-inverse matrices and exit\n";
    cout << "  --write-agh            Write A-inverse, G, and H-inverse during ssGWAS\n";
    cout << "  --ped <file>           Pedigree: ID SIRE DAM (no header)\n";
    cout << "  --genotyped-list <f>   File listing genotyped individual IDs\n";
    cout << "  --write-ebv            Reserved GEBV output switch (not stable in v1.0)\n";
    cout << "  --h2-tau <f>          Tau scaling for A22 blend (default: 1.0)\n";
    cout << "  --h2-omega <f>        Omega scaling for A22 blend (default: 1.0)\n";
    cout << "  --blend-g <f>         G blend weight with I (default: 0.05)\n";
    cout << "  Scan uses H⁻¹-estimated VCs with LOCO G.\n";
    cout << "\n--- Single-Step Method 2: ssGBLUP back-solving (Phase 6, experimental) ---\n";
    cout << "  --single-step-method lmm|backsolve\n";
    cout << "                        lmm       = Method 1 fixed-SNP association (default, backward compatible)\n";
    cout << "                        backsolve = Method 2: ssGBLUP genetic evaluation -> GEBV -> SNP\n";
    cout << "                                    back-solving -> per-SNP/window genetic variance\n";
    cout << "  --reml-vg <f> --reml-ve <f>\n";
    cout << "                        Fixed variance components for Method 2 (canonical/deterministic path)\n";
    cout << "  --ssgwas-vce <mode>   VCE estimation mode for Method 2 (ai|emai|em|hi|exact|...).\n";
    cout << "                        If omitted, --reml-vg/--reml-ve are required (fixed-VC path)\n";
    cout << "  --ssgwas-window-snps <n>\n";
    cout << "                        Sliding SNP-count window for %var aggregation (default: 10)\n";
    cout << "  --ssgwas-tune-g       Align G to A22 (G = b*G + a, preGSf90-style) before inversion\n";
    cout << "  --ssgwas-blend-a22    Blend G* with A22 instead of I (HIBLUP/preGSf90 reference; backsolve only)\n";
    cout << "  --ssgwas-solver ldlt|pcg\n";
    cout << "                        MME solver (default: ldlt reference; pcg = matrix-free PCG)\n";
    cout << "  --ssgwas-pvalue-method off|exact|apy-exact|apy-approx\n";
    cout << "                        Frequentist P-values for back-solved SNP effects (Phases C1/C2/D).\n";
    cout << "                        Default off. exact = C1 reference (selected MME inverse via LDLT,\n";
    cout << "                        O(n_genotyped^2) memory). apy-exact = APY route (Leite et al.\n";
    cout << "                        2024 GSE 56:59): APY G^-1 in the MME + exact core-slice selected\n";
    cout << "                        inverse (requires --ssgwas-solver ldlt; O(n_core^2) uncertainty\n";
    cout << "                        block; validated scalable APY uncertainty path).\n";
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
    cout << "  --ssgwas-apy-core-file <file>\n";
    cout << "                        Core-animal ID file (one ID per line, pedigree IDs); overrides\n";
    cout << "                        the seeded random core selection\n";
    cout << "  --ssgwas-apy-seed <s>  Seed of the deterministic random core selection (default 12345)\n";
    cout << "  --ssgwas-apy-core-target-ratio <x>\n";
    cout << "                        Harness-assigned label core_size / K98 (recorded in meta only; 0 = unset)\n";
    cout << "  --allow-large-exact-ssgwas-pvalue\n";
    cout << "                        Bypass the 8 GiB selected-inverse memory guard of\n";
    cout << "                        --ssgwas-pvalue-method exact (and apy-exact core slice)\n";
    cout << "  Method 2 outputs: <out>.ss_backsolve.ebv.txt / .snp.txt / .window.txt / .meta.txt\n";
    cout << "  With --ssgwas-pvalue-method exact additionally:\n";
    cout << "    <out>.ss_backsolve.pvalue.txt (EFFECT_ALLELE / SE_ALLELE / Z / P / LOG10P)\n";
    cout << "\n--- Pedigree Animal Model (Simple Animal Model) ---\n";
    cout << "  --animal-model        Pedigree-based Simple Animal Model GWAS: y = Xb + x*beta + Za + e,\n";
    cout << "                        a ~ N(0, A*sigma2_a) with full-pedigree A22 for the analysis sample set.\n";
    cout << "                        Requires --ped. A is built from the pedigree only (no genomic GRM);\n";
    cout << "                        LOCO does not apply (pedigree kernel is the same for all chromosomes).\n";
    cout << "                        Outputs *.assoc.txt (standard format) and *.animal_model.null.txt.\n";
    cout << "  --write-ped-kernel    Write the analysis A22 kernel (*.animal_model.k.txt + .k.id) for\n";
    cout << "                        external benchmarking (e.g. GEMMA -k).\n";

#ifdef COSMICSCANNER_ENABLE_EXPERIMENTAL
    cout << "\n--- Threshold/Probit Model (EXPERIMENTAL) ---\n";
    cout << "  --threshold-model      Binary trait GWAS (PQL + probit/logit)\n";
    cout << "  --link probit|logit    Link function (default: probit)\n";
    cout << "  PQL with fixed Ve=1 (probit) or Ve=pi^2/3 (logit).\n";
    cout << "  Warning: use --spa for improved calibration.\n";
    cout << "\n--- Count/Binomial PSR-GWAS Model (EXPERIMENTAL) ---\n";
    cout << "  --count-model          Count/binomial trait GWAS (PSR/PQL + weighted LMM)\n";
    cout << "  --family poisson|negbin|binomial Distribution (default: poisson)\n";
    cout << "  --offset <file>        One numeric log-offset/exposure value per sample for count traits\n";
    cout << "  --binomial-pheno <file> Binomial FID IID EVENT TRIAL file for event/trial traits\n";
    cout << "\n--- Ordinal PSR-GWAS Model (EXPERIMENTAL) ---\n";
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
    cout << "\n--- Survival (Cox) Model (EXPERIMENTAL) ---\n";
    cout << "  --survival-model       Time-to-event GWAS (Cox + martingale residual)\n";
    cout << "  --pheno format: FID IID TIME STATUS  (2 columns after FID/IID)\n";
    cout << "  Cox PH model fitted first, martingale residuals used in LMM.\n";
    cout << "\n--- GxE Reaction Norm Model (EXPERIMENTAL) ---\n";
    cout << "  --gxe-model            SNP-by-environment interaction GWAS\n";
    cout << "  --env-covar <col>      Environment variable column in --pheno\n";
    cout << "  VCE: Vg0 (intercept) + Vg1 (slope on env) + residual.\n";
    cout << "  Scan: 2-df Wald test for beta_g + beta_gxe jointly.\n";
    cout << "  Output: BETA_G, SE_G, BETA_GxE, SE_GxE, JOINT_CHISQ, JOINT_P\n";
#endif

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

#ifdef COSMICSCANNER_ENABLE_EXPERIMENTAL
    cout << "  --biobank-mode         REGENIE-style 2-step for UKB-scale data (EXPERIMENTAL)\n";
    cout << "  --regenie-model-snps <f> Model SNP list for Step 1\n";
    cout << "  --write-pred <file>    Write step 1 LOCO predictions\n";
    cout << "  --use-pred <file>      Use existing LOCO predictions\n";
    cout << "  --regenie-step <n>     0=both, 1=predictions only, 2=association only\n";
    cout << "  --predict-effect/-epm <f> GEMMA-style prediction from SNP effect table\n";
    cout << "  --predict-score <f>    Output prediction file (default: <out>.predict.txt)\n";
    cout << "  --predict-snp-col <c>  SNP ID column in effect table (default: SNP/rs/id)\n";
    cout << "  --predict-effect-col <c> Effect column (default: BETA/beta/alpha/effect)\n";
    cout << "  --bslmm/-bslmm <n>     GEMMA-style multi-locus model: 1 BSLMM, 2 ridge/GBLUP, 3 probit, 4/5 DAP (EXPERIMENTAL)\n";
    cout << "  --bslmm-top <n>        Sparse top-SNP count for --bslmm 1/3 approximation (default: 100)\n";
    cout << "  --bslmm-h2 <f>         Ridge/GBLUP shrinkage h2 in (0,1), default 0.5\n";
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
    cout << "  --skat                 Run SKAT rare-variant test (--set-file)\n";
    cout << "  --burden               Run Burden rare-variant test (--set-file)\n";
#else
    cout << "  --skat                 Run SKAT rare-variant test (--set-file)\n";
    cout << "  --burden               Run Burden rare-variant test (--set-file)\n";
#endif
    cout << "  --set-file <file>      Variant set definition for SKAT/Burden\n";
    cout << "\nModel notes:\n";
    cout << "  Models are mutually exclusive. Use exactly one --*-model flag.\n";
    cout << "  Default (no flag): standard univariate LMM (single-trait).\n";
    cout << "  Report bugs at: https://github.com/anomalyco/PanOmics/issues\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    GwasConfig cfg;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--input-config" && i + 1 < argc) {
            apply_input_config(cfg, argv[++i]);
        }
    }

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--bfile" && i + 1 < argc) { cfg.bfile_prefix = argv[++i]; }
        else if (arg == "--bgen" && i + 1 < argc) { cfg.bgen_file = argv[++i]; }
        else if (arg == "--sample" && i + 1 < argc) { cfg.sample_file = argv[++i]; }
        else if (arg == "--pfile" && i + 1 < argc) { cfg.pfile_prefix = argv[++i]; }
        else if (arg == "--pgen" && i + 1 < argc) { cfg.pgen_file = argv[++i]; }
        else if (arg == "--pvar" && i + 1 < argc) { cfg.pvar_file = argv[++i]; }
        else if (arg == "--psam" && i + 1 < argc) { cfg.psam_file = argv[++i]; }
        else if ((arg == "--out" || arg == "-o") && i + 1 < argc) { cfg.out_prefix = argv[++i]; }
        else if (arg == "--out-ssf") { cfg.out_ssf = true; }
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
        else if (arg == "--ssgwas-pvalue-method" && i + 1 < argc) {
            std::string v = argv[++i];
            if (v != "off" && v != "exact" && v != "apy-exact" && v != "apy-approx") {
                throw std::runtime_error("Invalid --ssgwas-pvalue-method value '" + v +
                                         "' (expected off|exact|apy-exact|apy-approx)");
            }
            cfg.ssgwas_pvalue_method = v;
        }
        else if (arg == "--ssgwas-apy-core-size" && i + 1 < argc) {
            int v = stoi(argv[++i]);
            if (v < 1 || v > 50000000) {
                throw std::runtime_error("Invalid --ssgwas-apy-core-size value '" + std::to_string(v) +
                                         "' (expected >= 1)");
            }
            cfg.ssgwas_apy_core_size = v;
        }
        else if (arg == "--ssgwas-apy-core-file" && i + 1 < argc) { cfg.ssgwas_apy_core_file = argv[++i]; }
        else if (arg == "--ssgwas-apy-seed" && i + 1 < argc) { cfg.ssgwas_apy_seed = stoi(argv[++i]); }
        else if (arg == "--ssgwas-apy-core-target-ratio" && i + 1 < argc) {
            cfg.ssgwas_apy_core_target_ratio = stod(argv[++i]);
        }
        else if (arg == "--allow-large-exact-ssgwas-pvalue") { cfg.allow_large_exact_ssgwas_pvalue = true; }
        else if (arg == "--ssgwas-vce" && i + 1 < argc) { cfg.ssgwas_vce = argv[++i]; }
        else if (arg == "--write-ebv") { cfg.write_ebv = true; }

#ifdef COSMICSCANNER_ENABLE_EXPERIMENTAL
        else if (arg == "--threshold-model") { cfg.threshold_model = true; cfg.run_gwas = true; }
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
#endif
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
        else if ((arg == "--model" || arg == "--gwas-method") && i + 1 < argc) {
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
        else if (arg == "--spa") { cfg.use_spa = true; }
        else if (arg == "--skat") { cfg.run_skat = true; }
        else if (arg == "--burden") { cfg.run_burden = true; }
        else if (arg == "--set-file" && i + 1 < argc) { cfg.set_file = argv[++i]; }

#ifdef COSMICSCANNER_ENABLE_EXPERIMENTAL
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
#endif
        else if (arg == "-rpace" && i + 1 < argc) { cfg.bslmm_rpace = stoi(argv[++i]); }
        else if (arg == "-wpace" && i + 1 < argc) { cfg.bslmm_wpace = stoi(argv[++i]); }
        else if (arg == "-w" && i + 1 < argc) { cfg.bslmm_burnin = stoi(argv[++i]); }
        else if (arg == "-s" && i + 1 < argc) { cfg.bslmm_sampling = stoi(argv[++i]); }
        else if (arg == "-seed" && i + 1 < argc) { cfg.bslmm_seed = stoll(argv[++i]); }
        else if (arg == "--chr" && i + 1 < argc) { cfg.scan_chr = argv[++i]; }
        else if (arg == "--from-bp" && i + 1 < argc) { cfg.scan_from_bp = stoll(argv[++i]); }
        else if (arg == "--to-bp" && i + 1 < argc) { cfg.scan_to_bp = stoll(argv[++i]); }
        else if (arg == "--no-covar-proj") { cfg.use_covariate_projection = false; }
        else if (arg == "--force-slow") { cfg.force_slow = true; }
        else if (arg == "--no-loco") { cfg.no_loco = true; }
        else if (arg == "--maf" && i + 1 < argc) { cfg.min_maf = stod(argv[++i]); }
        else if (arg == "--hwe" && i + 1 < argc) { cfg.min_hwe = stod(argv[++i]); }
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
        else { throw std::runtime_error("Unknown option: " + arg); }
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

    if (cfg.bfile_prefix.empty() && cfg.bgen_file.empty() && cfg.pfile_prefix.empty() && cfg.pgen_file.empty()) {
        print_usage(argv[0]);
        return 1;
    }

    GwasApp app(cfg);
    return app.run();
}
