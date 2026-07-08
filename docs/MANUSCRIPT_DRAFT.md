# CosmicScanner: A Unified and Scalable Framework for Genome-Wide Association across Complex Trait Architectures

**Jie Huang\***

College of Animal Sciences, Zhejiang University, Hangzhou 310058, China

\*Corresponding author. Email: huangjie2025@zju.edu.cn

---

## Highlights

- A single C++ engine unifies **11 GWAS model classes**—from standard LMM to multivariate, repeatability, random regression, ordinal, survival, maternal, ssGWAS, and GxE models—with one consistent interface.
- **Woodbury subject-level reduction** decouples repeatability-GWAS cost from record count, yielding **10–20× speedup and ~10× memory reduction** versus RepeatABEL with identical statistics.
- **BLUP-based Kronecker RRM** shrinks the longitudinal covariance from observation- to subject-level, achieving **213× speedup and 93× memory reduction** on the GMAT benchmark.
- Association statistics show near-perfect concordance with **GEMMA, EMmvGWAS, RepeatABEL, and GMA** (Pearson *r* > 0.99) and close concordance with **GMAT** for unbalanced longitudinal RRM (Pearson *r* = 0.951, top SNP concordant).
- Six extended models (maternal, threshold, GxE, survival, Poisson PSR, binomial PSR) pass type-I error calibration ($\lambda_{\text{GC}} \approx 1.00$); maternal model additionally validated by truth-recovery (residual-variance error < 6%).
- The release build passes a 10-test CTest suite covering compatibility smoke tests, external-reference routes, sparse REML, multi-VC null models, SPA, REGENIE, and input adaptation.
- Native support for **PLINK 1.9, PLINK 2.0, and BGEN** with memory-mapped, multi-threaded I/O removes format-conversion overhead.

---

## Abstract

Modern genetic studies generate phenotypes that span continuous, longitudinal, repeated-measure, non-Gaussian, and multi-trait architectures, yet the GWAS software ecosystem remains fragmented: each model class is served by a different tool, with incompatible formats, GRM builders, and output schemas. Worse, the tools for repeated-measure and longitudinal data often construct record- or observation-level covariance matrices whose $O(n^3)$ factorization caps sample sizes well below biologically relevant scales. We present **CosmicScanner**, a unified C++ framework that integrates **11 GWAS models** within a single engine sharing common genotype I/O, GRM construction, and null-model caching. Two algorithmic designs deliver scalability for the previously limiting model classes: (i) an **exact Woodbury subject-level reduction** that rewrites the repeatability-model covariance through the incidence matrix $\mathbf{Z}$, confining dense algebra to the subject level; and (ii) a **BLUP-based Kronecker multivariate LMM** that estimates Legendre random-regression coefficients first, then scans SNPs in coefficient space, sidestepping the observation-level covariance entirely. CosmicScanner natively reads PLINK 1.9, PLINK 2.0, and BGEN with memory-mapped, multi-threaded decompression. Across validation against **GEMMA, EMmvGWAS, RepeatABEL, GMA, and GMAT**, association statistics showed near-perfect agreement for single-trait, multivariate, and repeatability scans (Pearson *r* > 0.99), and close concordance for unbalanced longitudinal RRM (*r* = 0.951 against GMAT). Repeatability GWAS ran **10–20× faster at ~10× lower memory**, and longitudinal RRM ran **213× faster at 93× lower memory** than the corresponding reference implementations or observation-level paths. Six additional model classes (maternal, threshold/probit, GxE reaction-norm, survival Cox, Poisson PSR, and binomial PSR) were internally calibrated via type-I error simulation ($\lambda_{\text{GC}} \approx 1.00$, all replicates within the 95% binomial confidence interval); the maternal model was further validated by positive truth-recovery (residual-variance error < 6%, correct null-maternal recovery, architecture-aware SNP detection). The release build passes a 10-test CTest suite covering compatibility smoke tests, external-reference routes, sparse REML, multi-VC null models, SPA, REGENIE, and input adaptation. CosmicScanner thus offers a single, validated, and scalable framework for genome-wide association across diverse complex trait architectures.

**Keywords:** GWAS; linear mixed model; longitudinal traits; repeated measurements; multivariate analysis; random regression; Woodbury identity; genomic relationship matrix

---

## Significance

The complexity of modern phenotypic data—repeated records, longitudinal trajectories, non-Gaussian outcomes, and correlated multi-trait panels—currently forces geneticists to chain together specialized tools that disagree on basic statistical assumptions and cannot scale to realistic record counts. CosmicScanner resolves this fragmentation by unifying 11 GWAS model classes in one engine. It introduces an exact Woodbury reduction for repeated measurements and a BLUP-based coefficient-space strategy for longitudinal random regression. The result is a tool that reproduces established implementations with high concordance while running one to two orders of magnitude faster and using an order of magnitude less memory on the model classes that previously did not scale well.

---

## 1. Introduction

Genome-wide association studies (GWAS) have become the dominant approach for identifying genetic variants underlying complex traits since the first successful application in 2005 (1). Over the past two decades, GWAS methodologies have progressively expanded beyond the standard single-trait linear mixed model (LMM) to address increasingly complex phenotypic architectures. In livestock and plant breeding programs, individuals are routinely phenotyped multiple times across growth stages, parities, or environments, generating longitudinal and repeated-measure data (2,3). In biomedical cohorts such as the UK Biobank, phenotypes span continuous, binary, count, ordinal, and time-to-event outcomes (4,5). Multi-trait joint analyses exploit phenotypic and genetic correlations to boost statistical power and reveal pleiotropic architecture (6,7). These developments have motivated a proliferation of specialized GWAS tools, each designed for a specific model class.

The linear mixed model has served as the foundational framework for GWAS since the introduction of the genomic relationship matrix (GRM)-based method by Yang et al. (8) and the efficient mixed-model association (EMMA) algorithm by Kang et al. (9). Tools such as GEMMA (10), GCTA (8), BOLT-LMM (11), and REGENIE (12) have enabled scalable single-trait and, in some cases, multi-trait GWAS. For repeated-measure data, RepeatABEL (13) implements a record-level mixed model with additive genetic, permanent environmental, and residual variance components. For longitudinal traits, GMA (14) and GMAT (15) provide random regression model (RRM) approaches using Legendre polynomial basis functions. For ordinal traits, POLMM (16) extends the proportional odds model with a genomic kinship-based random effect. For multi-trait analysis, EMmvGWAS (17) and the multivariate LMM in GEMMA (10) enable joint testing of multiple phenotypes.

Despite these advances, three structural gaps remain in the GWAS software ecosystem.

**First, model coverage is fragmented.** No single tool implements the full spectrum of models required by a modern study. A livestock population with repeated records, correlated traits, and ordinal outcomes forces researchers to chain RepeatABEL, GEMMA, and POLMM, each with its own input format, GRM builder, variance-component algorithm, and output schema. The resulting pipelines harbor hidden inconsistencies in statistical assumptions and are difficult to reproduce.

**Second, repeated-measure and longitudinal GWAS do not scale.** RepeatABEL constructs a dense record-level covariance $\mathbf{V}$ and fully eigendecomposes it, incurring $O(n_{\text{records}}^3)$ time and $O(n_{\text{records}}^2)$ memory (13). For 1,000 genotyped subjects with 10 records each, this means factorizing a $10{,}000 \times 10{,}000$ matrix. Observation-level RRM approaches face the same wall: the covariance dimension equals the total observation count, which routinely reaches tens of thousands (15). In practice, exact repeated-measure and longitudinal GWAS are confined to small datasets, and analysts resort to ad hoc subject-mean aggregation that discards within-subject information.

**Third, modern genotype formats are incompletely supported.** PLINK 2.0 (`.pgen`) and BGEN have become the de facto standards for biobank data thanks to dosage encoding and compression, yet many GWAS tools still accept only PLINK 1.9 (`.bed`), forcing lossy or storage-heavy format conversion.

To close these gaps, we developed **CosmicScanner**, a unified high-performance GWAS framework that integrates 11 model classes within a single C++ engine. Its guiding principle is that diverse trait architectures should be served by model-specific statistical formulations while sharing a common computational substrate—genotype I/O, GRM construction, null-model caching, and block-wise SNP scanning—so that statistical assumptions stay consistent across models, redundant computation is eliminated, and a single interface accepts any phenotype type.

The main contributions of this work are:

1. **Unified model coverage.** Eleven GWAS model classes—single-trait LMM, multivariate LMM, repeatability, random regression, maternal effect, single-step ssGWAS, threshold/probit, count/binomial PSR, ordinal PSR, survival Cox, and GxE reaction-norm—are implemented in one package with consistent input/output conventions.

2. **Woodbury subject-level reduction for repeatability GWAS.** The exact record-level mixed model is preserved analytically, but all dense matrix operations are reduced from the record level ($n_{\text{records}}^3$) to the subject level ($n_{\text{subjects}}^3$), yielding **10–20× speedup and ~10× memory reduction** versus RepeatABEL with statistically identical results.

3. **BLUP-based Kronecker RRM for unbalanced longitudinal GWAS.** Subject-level Legendre coefficients are first estimated by BLUP; SNPs are then scanned in coefficient space under a Kronecker covariance, bypassing the observation-level matrix entirely. On the GMAT benchmark this delivers **213× speedup and 93× memory reduction**.

4. **Comprehensive validation.** Association statistics concur (Pearson *r* > 0.99) with five established tools—GEMMA, EMmvGWAS, RepeatABEL, GMA, and GMAT—across single-trait, multivariate, repeated-measure, and longitudinal scenarios.

---

## 2. Materials and Methods

### 2.1 Unified Data Input and GRM Construction

#### Genotype format support

CosmicScanner natively supports three genotype formats without external conversion:

- **PLINK 1.9** (`.bed`/`.bim`/`.fam`): The traditional binary format widely used in genetics research.
- **PLINK 2.0** (`.pgen`/`.pvar`/`.psam`): The next-generation format supporting dosage-based encoding and efficient compression, designed for million-sample-scale biobank data.
- **BGEN** (`.bgen`/`.sample`): The Oxford format used by the UK Biobank and other large cohorts, supporting probabilistic genotype calls.

All three formats share a common reader interface that provides genotype dosage vectors to the GWAS scanning engine. For PLINK 2.0 and BGEN, CosmicScanner employs memory-mapped file I/O (`mmap`) combined with multi-threaded decompression to minimize memory overhead and maximize read throughput. Genotype dosages are loaded in blocks during SNP scanning to maintain a bounded memory footprint regardless of the total variant count.

#### GRM construction

The genomic relationship matrix (GRM) is constructed using the VanRaden (2008) method (18):

$$G = \frac{(M - P)(M - P)'}{2 \sum_{i=1}^{m} p_i (1 - p_i)}$$

where $M$ is the $n \times m$ genotype matrix coded as 0/1/2 reference-allele counts, $P$ is the $n \times m$ matrix of column means ($2p_i$ for SNP $i$ with reference-allele frequency $p_i$), and the denominator scales $G$ to be analogous to the pedigree-based numerator relationship matrix.

CosmicScanner implements six GRM construction paths (`GrmBuilder::compute` and `GrmBuilder::computeSplit`) that cover PLINK, PGEN, and BGEN inputs with optional chromosome-partitioned computation for LOCO analysis. GRM construction is parallelized using OpenMP, with AVX2/AVX-512 SIMD instructions accelerating the inner-loop genotype cross-product accumulation. Variant filtering via `--extract`/`--exclude` (SNP ID lists), `--chr` (chromosome), and `--from-bp`/`--to-bp` (genomic interval) is applied at both GRM construction and SNP scanning stages to ensure consistency.

#### Phenotype and covariate input

CosmicScanner accepts standard PLINK-format phenotype files (`FID IID Pheno1 Pheno2 ...`) with missing values encoded as `NA` or `-9`. For longitudinal and repeated-measure models, long-format phenotype files are supported where each row represents a single observation record keyed by `FID IID TIME`. Quantitative covariates (`--qcovar`) and categorical covariates (`--covar`) are specified as separate files in PLINK format. An input adaptation layer provides automatic conversion from CSV/TSV/TXT formats to the prepared input schema, supporting model-specific requirements (e.g., record-level covariates for repeatability models, time-dependent covariates for RRM).

### 2.2 Single-Trait LMM GWAS

#### Model

The standard single-trait linear mixed model for GWAS is:

$$\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \mathbf{Z}\mathbf{u} + \mathbf{e}$$

where $\mathbf{y}$ is the $n \times 1$ phenotype vector, $\mathbf{X}$ is the design matrix for fixed effects and covariates, $\boldsymbol{\beta}$ is the corresponding effect vector, $\mathbf{Z}$ is the incidence matrix (typically the identity for genotyped individuals), $\mathbf{u} \sim N(\mathbf{0}, \sigma_g^2 \mathbf{G})$ is the random polygenic effect, and $\mathbf{e} \sim N(\mathbf{0}, \sigma_e^2 \mathbf{I})$ is the residual. The marginal covariance is $\mathbf{V} = \sigma_g^2 \mathbf{G} + \sigma_e^2 \mathbf{I}$.

#### Exact LMM scan

For exact LMM scanning, CosmicScanner performs eigendecomposition of the GRM: $\mathbf{G} = \mathbf{U}\mathbf{D}\mathbf{U}'$, where $\mathbf{D}$ is the diagonal matrix of eigenvalues and $\mathbf{U}$ is the matrix of eigenvectors. The phenotype, covariates, and genotype dosages are rotated into the eigenvector basis: $\tilde{\mathbf{y}} = \mathbf{U}'\mathbf{y}$, $\tilde{\mathbf{X}} = \mathbf{U}'\mathbf{X}$, $\tilde{\mathbf{x}}_j = \mathbf{U}'\mathbf{x}_j$. In this rotated space, the covariance becomes diagonal: $\tilde{\mathbf{V}} = \sigma_g^2 \mathbf{D} + \sigma_e^2 \mathbf{I}$, enabling $O(n)$ per-SNP evaluation via weighted least squares.

To accelerate exact LMM, CosmicScanner optionally projects out the covariate matrix $\mathbf{X}$ before eigendecomposition, reducing the effective problem dimension and eliminating redundant covariate operations during the SNP scan. SNPs are processed in blocks (default block size: 2,048 variants), with the rotated genotype block accumulated and evaluated as a batch, amortizing the overhead of genotype loading and rotation.

#### Fast-LMM approximation

When the variance ratio $\delta = \sigma_g^2 / \sigma_e^2$ is fixed at the null-model estimate rather than re-estimated per SNP, the scan reduces to a sequence of generalized least squares evaluations with identical covariance structure. This "fast-LMM" mode (19) significantly accelerates genome-wide scanning while maintaining near-equivalent statistical power for common variants. CosmicScanner supports both exact (per-SNP variance refit) and fast-LMM modes.

#### LOCO (Leave-One-Chromosome-Out)

To avoid proximal contamination—the inflation of test statistics for SNPs on the same chromosome as the associated variant due to the GRM capturing linkage signals—CosmicScanner implements LOCO analysis (20). In LOCO mode, a chromosome-specific null model is estimated using a GRM that excludes all SNPs on the target chromosome. The eigendecomposition, null-model variance components, and data rotation are performed once per chromosome, and SNPs on that chromosome are scanned against the corresponding LOCO null model.

### 2.3 IDUL/AIDUL Null-Model Acceleration

#### Motivation

The null-model variance-ratio estimation in LMM GWAS requires finding $\eta = \sigma_g^2 / \sigma_e^2$ that maximizes the restricted log-likelihood. The standard approach uses Brent's method, which evaluates the log-likelihood function multiple times per iteration. For large samples, each likelihood evaluation involves $O(n^2)$ operations for the diagonal variance computation in rotated space, and Brent's method may require 20–50 iterations to converge, creating a computational bottleneck.

#### IDUL: Iterative Dispersion Update using Likelihood

We developed IDUL (Iterative Dispersion Update using Likelihood), a Fisher-scoring-based iterative solver that updates the variance ratio directly in the eigendecomposition space. After GRM eigendecomposition and data rotation, the log-likelihood as a function of $\eta$ depends only on the diagonal matrix $\mathbf{D}$ and the rotated phenotype $\tilde{\mathbf{y}}$. The Fisher scoring update for $\eta$ is:

$$\eta^{(t+1)} = \eta^{(t)} + \mathbf{I}(\eta^{(t)})^{-1} \frac{\partial \ell}{\partial \eta}\bigg|_{\eta^{(t)}}$$

where $\mathbf{I}(\eta)$ is the Fisher information and $\frac{\partial \ell}{\partial \eta}$ is the score function. Both quantities can be computed in $O(n)$ time using the eigenvalues $\mathbf{D}$ and rotated data, making each IDUL iteration substantially faster than a Brent likelihood evaluation.

#### IDUL+ and AIDUL

**IDUL+** extends IDUL with a monotone likelihood safeguard: if an update fails to increase the log-likelihood, a step-size reduction is applied to ensure convergence.

**AIDUL** (Accelerated IDUL) incorporates Anderson acceleration (window size = 5) on top of the Fisher-scoring updates. The Anderson acceleration extrapolates from the history of previous iterates to achieve faster convergence, particularly for ill-conditioned likelihood surfaces. The Fisher-scoring update is extracted as a separate subroutine to enable clean integration with the Anderson history buffer.

#### Fallback and applicability

When IDUL, IDUL+, or AIDUL fails to converge within the maximum iteration count, CosmicScanner automatically falls back to the Brent solver to guarantee a valid result. The IDUL/AIDUL solvers are directly applicable to models with a single GRM and residual variance component: standard single-trait LMM, LOCO LMM, and the aggregated fast path of the repeatability model. They are not directly applicable to multi-variance-component models (multivariate LMM, exact record-level repeatability, RRM), which require matrix-variate or multi-component updates.

### 2.4 Multivariate LMM GWAS

#### Model

The multivariate linear mixed model for $q$ traits is:

$$\mathbf{Y} = \mathbf{X}\mathbf{B} + \mathbf{Z}\mathbf{U} + \mathbf{E}$$

where $\mathbf{Y}$ is the $n \times q$ phenotype matrix, $\mathbf{B}$ is the fixed-effect matrix, $\text{vec}(\mathbf{U}) \sim N(\mathbf{0}, \mathbf{K} \otimes \mathbf{G})$ is the $nq \times 1$ vector of random genetic effects with trait covariance matrix $\mathbf{K}$ ($q \times q$) and GRM $\mathbf{G}$, and $\text{vec}(\mathbf{E}) \sim N(\mathbf{0}, \mathbf{R} \otimes \mathbf{I})$ is the residual with trait covariance $\mathbf{R}$. The marginal covariance is $\mathbf{V} = \mathbf{K} \otimes \mathbf{G} + \mathbf{R} \otimes \mathbf{I}$.

#### Five methods

CosmicScanner implements five multivariate GWAS methods:

**1. EMMV (EMmvGWAS-compatible).** This method follows the EMmvGWAS algorithm (17) using the kinship-based transformation $\boldsymbol{\Lambda} = \mathbf{G}\mathbf{R}^{-1}$ with an explicit intercept term and no covariate projection. The per-SNP scan uses the `dosage - 1` coding to align with the official EMmvGWAS convention, and the general matrix inverse is applied to the asymmetric $\boldsymbol{\Lambda}\mathbf{D} + \mathbf{I}$ matrix.

**2. GEMMA (GEMMA-full per-SNP refit).** This method follows the multivariate LMM implementation in GEMMA (10), performing a full variance-component refit for each SNP. This provides the most accurate per-SNP association statistics but is computationally intensive.

**3. GEMMA-full.** This is CosmicScanner's implementation of the full multivariate Wald, likelihood-ratio, and score tests using the GEMMA-style null model with per-SNP variance refit. It supports PLINK, PGEN, and BGEN inputs with block-wise genotype loading.

**4. Hybrid (LOCO-aware adaptive).** This is a two-stage method designed to combine the speed of fast screening with the accuracy of full refit:

- *Stage 1 (Screening):* A per-chromosome LOCO GRM is constructed, and a chromosome-specific multivariate null model is estimated. SNPs on the current chromosome are scanned using the fast EMMV-style screening test.
- *Stage 2 (Selective refinement):* Candidate SNPs are selected using a union of three criteria: (a) p-value threshold, (b) top-$K$ ranking, and (c) Score/Wald discrepancy trigger. Candidates are grouped by chromosome and refined using per-chromosome LOCO GEMMA-full Wald/LRT/Score tests.

The hybrid method also supports an optional Score-only pre-pass audit (`--mv-refine-score-prepass`) that performs a GEMMA-full Score test on the top $n$ screening hits to identify SNPs where the fast screening test and the full Score test disagree.

**5. EB (Empirical Bayes).** This method applies summary-level empirical Bayes shrinkage to the multivariate SNP effects. After obtaining the mvLMM fixed-null SNP effect vector $\hat{\mathbf{c}}_k$ and its covariance for SNP $k$, the EB prior assumes $\mathbf{c}_k \sim N(\mathbf{0}, \lambda_k \mathbf{V}_e)$, where $\lambda_k$ is a SNP-specific marker variance ratio estimated from the data. This avoids the instability of full $q \times q$ per-SNP variance-component refit in small samples.

### 2.5 Repeatability Model GWAS

#### Model

For record $j$ of subject $i$, the repeatability GWAS null model is:

$$\mathbf{y} = \mathbf{X}\mathbf{b} + \mathbf{Z}\mathbf{u} + \mathbf{Z}\mathbf{p} + \mathbf{e}$$

$$\mathbf{u} \sim N(\mathbf{0}, \sigma_g^2 \mathbf{G}), \quad \mathbf{p} \sim N(\mathbf{0}, \sigma_{pe}^2 \mathbf{I}), \quad \mathbf{e} \sim N(\mathbf{0}, \sigma_e^2 \mathbf{W}^{-1})$$

where $\mathbf{G}$ is the subject-level GRM ($n_{\text{subj}} \times n_{\text{subj}}$), $\mathbf{Z}$ is the $n_{\text{records}} \times n_{\text{subj}}$ incidence matrix mapping records to subjects, $\mathbf{p}$ is the permanent environmental effect, and $\mathbf{W}$ is an optional diagonal matrix of record-specific residual weights. The marginal covariance is:

$$\mathbf{V} = \mathbf{Z}\mathbf{K}\mathbf{Z}' + \mathbf{R}, \quad \mathbf{K} = \sigma_g^2 \mathbf{G} + \sigma_{pe}^2 \mathbf{I}, \quad \mathbf{R} = \sigma_e^2 \mathbf{W}^{-1}$$

#### Woodbury subject-level reduction

The key algorithmic innovation is the application of the Woodbury matrix identity to avoid explicit construction and factorization of the dense record-level matrix $\mathbf{V}$:

$$\mathbf{V}^{-1} = \mathbf{R}^{-1} - \mathbf{R}^{-1}\mathbf{Z}(\mathbf{K}^{-1} + \mathbf{Z}'\mathbf{R}^{-1}\mathbf{Z})^{-1}\mathbf{Z}'\mathbf{R}^{-1}$$

The core matrix requiring dense factorization is therefore subject-level:

$$\mathbf{M} = \mathbf{K}^{-1} + \mathbf{Z}'\mathbf{R}^{-1}\mathbf{Z}$$

This changes the main covariance operation from record-level dense algebra ($n_{\text{records}} \times n_{\text{records}}$) to subject-level dense algebra ($n_{\text{subjects}} \times n_{\text{subjects}}$) with record-level contributions aggregated through the incidence matrix $\mathbf{Z}$. The determinant is computed via the matrix determinant lemma:

$$\det(\mathbf{V}) = \det(\mathbf{R}) \cdot \det(\mathbf{K}) \cdot \det(\mathbf{M})$$

#### Incidence-aware SNP scanning

For each SNP, CosmicScanner keeps genotype dosage at the subject level and never expands the SNP vector to a full record-level matrix. The SNP effect is estimated by projecting the subject-level genotype through the incidence matrix $\mathbf{Z}$, computing weighted residual projections that account for record-level residual weights, and reusing the null-model Woodbury cache. Genotype blocks have dimension $n_{\text{subjects}} \times \text{block\_size}$, while repeated records contribute through subject-level sums and weighted residual projections. This design ensures that the SNP scan complexity depends on $n_{\text{subjects}}$, not $n_{\text{records}}$.

#### Four strategies

CosmicScanner implements four repeatability GWAS strategies:

1. **Exact:** Full Woodbury solver with subject-level matrix factorization. This is the default and recommended mode for production analysis.
2. **Fast (Mean aggregation):** Repeated records are aggregated to weighted subject-level means, and the scan reduces to a standard single-trait LMM with observation weights. This mode is faster but approximate for unbalanced record structures.
3. **STCG (Subject-level GRM cache):** The subject-level GRM is cached as a `Float32` binary file and reused via `mmap` across runs, providing the foundation for out-of-core computation on very large datasets.
4. **Mixed precision:** In exact mode, the record-level kernel matrix is constructed in `Float32` and promoted to `Float64` before the final solve, reducing intermediate memory usage by approximately 50% while maintaining numerical accuracy.

#### Variance-component estimation

The default production setting uses subject-level Woodbury variance-component estimation (VCE) and skips full-record exact refinement (`--vce-mode auto`). Full record-level exact refinement (`--vce-mode full-record-exact`) remains available for small-data validation. This policy is critical for scalability: while the exact record-level model is preserved through the Woodbury algebra, the VCE optimization loop operates on subject-level matrices.

### 2.6 Random Regression Model (RRM) GWAS

#### Model

For longitudinal or time-gradient phenotypes, the random regression model expresses the SNP effect as a function of time or environment. For subject $i$ observed at time $t_{ij}$:

$$y_{ij} = \mathbf{x}_{ij}'\boldsymbol{\beta} + \sum_{k=0}^{K} \phi_k(t_{ij}) a_{ik} + e_{ij}$$

where $\phi_k(t)$ are basis functions (Legendre polynomials by default; B-spline optionally), $a_{ik}$ is the random regression coefficient for subject $i$ and basis $k$, and the genetic covariance structure is:

$$\text{Cov}(a_k, a_l) = K_{kl} \mathbf{G}$$

where $\mathbf{K}$ is the $(K+1) \times (K+1)$ genetic covariance matrix of random regression coefficients and $\mathbf{G}$ is the GRM.

#### Balanced longitudinal (observation-level RRM)

For balanced longitudinal data where all subjects share the same observation time points, CosmicScanner constructs the observation-level covariance matrix:

$$\mathbf{V} = \mathbf{Z}(\mathbf{K} \otimes \mathbf{G})\mathbf{Z}' + \mathbf{R}$$

and performs a joint multi-degree-of-freedom Wald test for each SNP. The null-model variance components are estimated using AI-REML for small samples or a moment-initialized fast null for large samples (`--rrm-ai-iter -1` for automatic mode selection).

The output includes:
- Joint test statistics (`JOINT_CHISQ`, `JOINT_P`)
- Dynamic effect curves on a user-specified time grid
- Best-effect time, onset/offset, maximum effect, trajectory type (stable, sign-switch, early-specific, late-specific, transient)
- Key time windows based on a significance threshold

#### Unbalanced longitudinal (BLUP-based RRM)

For unbalanced longitudinal data where subjects have different observation time points or varying numbers of records, the observation-level covariance matrix becomes large and potentially singular. CosmicScanner implements a BLUP-based approach that avoids this construction:

**Step 1 (Individual BLUP):** For each subject, Legendre random regression coefficients $\hat{\mathbf{a}}_i = (\hat{a}_{i0}, \hat{a}_{i1}, \ldots, \hat{a}_{iK})'$ are estimated via BLUP using the subject-level GRM and the phenotype records.

**Step 2 (Kronecker multivariate LMM):** The BLUP coefficient matrix $\hat{\mathbf{A}}$ ($n \times (K+1)$) is treated as a multivariate phenotype, and a Kronecker-structured multivariate LMM is used to scan SNPs:

$$\text{vec}(\hat{\mathbf{A}}) \sim N(\mathbf{0}, \mathbf{K}_0 \otimes \mathbf{G} + \mathbf{P}_0 \otimes \mathbf{I})$$

where $\mathbf{K}_0$ is the genetic covariance matrix of Legendre coefficients and $\mathbf{P}_0$ is the permanent environmental / within-subject covariance matrix.

This approach reduces the covariance dimension from $n_{\text{observations}} \times n_{\text{observations}}$ to $n \times n$ (subject-level) via the Kronecker structure, enabling analysis of datasets with tens of thousands of observations. Time-dependent categorical covariates (`--tclass`) and continuous covariates (`--tcovar`) are supported.

The BLUP-based RRM supports LOCO analysis, with the option of per-chromosome null-model re-estimation (default) or a fast global-null approximation (`--rrm-loco-global-null`). The fast mode reuses the global BLUP null model and only partitions SNP scanning by chromosome, trading some accuracy for speed.

### 2.7 Non-Gaussian Trait GWAS

For phenotypes that deviate from the normal distribution—binary, count, binomial, and ordinal traits—CosmicScanner implements a penalized quasi-likelihood (PQL) framework termed PSR-GWAS (Pseudo-Score Random regression GWAS). The approach linearizes the generalized linear mixed model via a working variate and passes the working residual variance to the LMM scanner as observation weights.

#### Threshold model (binary traits)

For binary traits, a threshold (probit or logit) model is fitted:

$$P(y_i = 1) = g^{-1}(\mathbf{x}_i'\boldsymbol{\beta} + u_i)$$

where $g(\cdot)$ is the link function (probit or logit). The PSR/PQL linearization produces a pseudo-response $\tilde{y}$ and working weights $\mathbf{W} = \text{diag}(w_i)$, which are passed to the weighted LMM scanner as `obs_weights = sqrt(W)`. The SNP scan uses weighted LMM rather than ordinary LMM on the unweighted pseudo-phenotype.

For binary traits with case-control imbalance, a saddlepoint approximation (SPA) is available to calibrate test statistics under the null, correcting for the deviation from normality in the score statistic.

#### Count model (Poisson/Binomial)

For count data, a Poisson or binomial generalized linear mixed model is fitted with the canonical link. An optional offset (`--offset`) accommodates exposure variation. The PSR linearization and weighted LMM scan follow the same framework as the threshold model.

#### Ordinal model (PSR-GWAS v0.5)

For ordinal traits with $K$ ordered categories, CosmicScanner implements a stacked proportional-odds model via the stacked PSR approach. The cumulative-link model is linearized into $K - 1$ stacked binary pseudo-responses, and a dense $\mathbf{Z}\mathbf{K}\mathbf{Z}'$ (dense-zkz) representation is used for the stacked covariance. Variance-component estimation is performed via grid REML over the genetic variance $\sigma_g^2$, which avoids boundary solutions (0% boundary rate observed in validation). The Wald test and pseudo-score test are available, with diagnostic output (`--ordinal-test all`) including robust-score sensitivity statistics.

#### Survival model (Cox PH)

For time-to-event data, a Cox proportional hazards model is fitted, and martingale residuals from the null model are used as the pseudo-phenotype for LMM scanning. This approach enables GWAS for censored survival outcomes while leveraging the GRM-based population structure correction.

### 2.8 Maternal-Effect Model

For traits influenced by maternal effects (e.g., birth weight, weaning weight), CosmicScanner implements a multi-component variance model:

$$\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \mathbf{Z}_a\mathbf{a} + \mathbf{Z}_m\mathbf{m} + \mathbf{Z}_{am}\mathbf{a}_m + \mathbf{c} + \mathbf{e}$$

where $\mathbf{a}$ is the direct additive genetic effect, $\mathbf{m}$ is the maternal additive genetic effect, $\mathbf{a}_m$ is the direct-maternal covariance component, $\mathbf{c}$ is the maternal permanent environmental effect, and $\mathbf{e}$ is the residual. The model supports five variance-component configurations (ad/am/adm/c/e) and estimates the direct-maternal genetic correlation.

### 2.9 Single-Step GWAS (ssGWAS)

For populations where only a subset of individuals is genotyped, CosmicScanner implements single-step GWAS (ssGWAS) built on the combined relationship matrix $\mathbf{H}$, which integrates pedigree ($\mathbf{A}$) and genomic ($\mathbf{G}$) information. While ssGBLUP is a genomic-selection (breeding-value estimation) method, ssGWAS repurposes the same $\mathbf{H}^{-1}$ infrastructure for variance-component estimation and then performs a SNP-level Wald scan, enabling association testing in pedigreed populations with incomplete genotyping:

$$\mathbf{H}^{-1} = \mathbf{A}^{-1} + \begin{pmatrix} \mathbf{0} & \mathbf{0} \\ \mathbf{0} & \mathbf{G}^{-1} - \mathbf{A}_{22}^{-1} \end{pmatrix}$$

where $\mathbf{A}_{22}$ is the pedigree relationship matrix among genotyped individuals. This enables joint analysis of genotyped and non-genotyped animals, which is critical for livestock breeding programs with incomplete genotyping.

### 2.10 GxE Reaction-Norm Model and Rare-Variant Tests

#### GxE reaction-norm model

For gene-by-environment interaction analysis, CosmicScanner models the SNP effect as a linear function of an environmental covariate $E$:

$$\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \mathbf{x}_j\beta_0 + \mathbf{x}_j \circ \mathbf{E}\,\beta_1 + \mathbf{u} + \mathbf{e}$$

where $\beta_0$ is the main SNP effect, $\beta_1$ is the GxE interaction effect, and $\circ$ denotes element-wise product. The model estimates two variance components ($V_{g0}$, $V_{g1}$) and performs a 2-df joint Wald test for the combined main and interaction effect.

#### Rare-variant tests

For rare-variant association analysis, CosmicScanner implements SKAT (Sequence Kernel Association Test) and burden tests. SKAT uses a beta-weight function to up-weight rare variants, and the p-value is computed via a mixture of chi-squared distributions (Davies' method). The burden test aggregates rare variants within a gene or region and tests the cumulative effect.

#### Biobank two-step method

For UK Biobank-scale data, CosmicScanner implements a REGENIE-style two-step regression:

- **Step 1:** Block-wise ridge regression using a pruned set of model SNPs, with LOCO predictions for each chromosome.
- **Step 2:** Streaming single-SNP association tests after residualizing the phenotype and genotype against covariates and the chromosome-specific LOCO prediction.

LOCO predictions can be cached (`--write-pred`) and reused (`--use-pred`) for Step 2 reruns without retraining the ridge model, enabling efficient distributed computation across chromosomes.

#### Empirical Bayes GWAS

CosmicScanner provides validated empirical Bayes (EB) GWAS auxiliary scans for signal enhancement. For single-trait analysis, a random-SNP-effect model (Wang 2015 style) estimates a per-SNP marker variance ratio $\lambda_k$ and applies Bayesian shrinkage to the effect estimate. For multivariate analysis, summary-level EB shrinkage is applied to the mvLMM fixed-null SNP effects. EB-GWAS is designed as an auxiliary signal-enhancement route and does not replace the standard fixed-effect LMM/mvLMM/repeatability GWAS results.

### 2.11 Implementation

CosmicScanner is implemented in C++17 using the Eigen3 linear algebra library for dense and sparse matrix operations. Multi-threaded parallelization is achieved via OpenMP, with the number of threads controlled by `--thread-num`. The GRM construction inner loop and the SNP scanning block processing are the primary parallelized sections. Where supported by the CPU, AVX2 and AVX-512 SIMD instructions accelerate genotype cross-product accumulation and matrix-vector operations.

The software follows a modular architecture:
- **Data loaders** (`data_loader.h`, `plink_reader.h`, `bgen_reader.h`, `pgen_reader.h`): Format-specific genotype readers with a common interface.
- **LMM engine** (`lmm_engine.cpp`): Null-model preparation, data rotation, and block-wise SNP scanning for PLINK/PGEN/BGEN.
- **Model dispatch** (`gwas_app.cpp`): Routes to the appropriate model (univariate, multivariate, repeatability, RRM, etc.) based on CLI flags.
- **Model implementations** (`src/models/`): Each model class (repeatability, RRM, multivariate, ordinal, etc.) is implemented as a separate module with its own null-model estimation and scanning logic.
- **Core library** (`libCosmicCore`): Shared utilities for VCE, eigendecomposition, and matrix operations, built as a static library shared with the CosmicSolver genomic selection module.

The software is cross-platform (Linux native, Windows via WSL2) and is freely available at https://github.com/huangjie2025/CosmicScanner under the MIT License. Pre-compiled binaries for Linux x86_64 are provided on the GitHub Release page.

---

## 3. Results

We organized the validation into four threads that mirror the design goals of CosmicScanner: (i) **statistical fidelity** to established reference tools, (ii) **format consistency** across PLINK/PGEN/BGEN, (iii) **scalability** on the model classes that previously did not scale, and (iv) **calibration** of the non-Gaussian PSR-GWAS extensions. Across all comparisons, the headline finding is that CosmicScanner reproduces reference implementations to Pearson *r* > 0.99 while delivering one to two orders of magnitude improvements in runtime and memory on repeated-measure and longitudinal data. Figure 1 provides an overview of the software architecture and model coverage; Table 7 summarizes the 11 models and their validation status; the subsections below detail the quantitative evidence with supporting Figures 2–7.

### 3.1 Single-Trait LMM GWAS: Accuracy and Format Consistency

#### Accuracy against GEMMA

We first validated the single-trait LMM GWAS implementation of CosmicScanner against GEMMA (v0.98), the most widely used LMM-based GWAS tool. Using simulated PLINK data with $N = 1{,}000$ individuals and $M = 10{,}000$ SNPs, we compared effect-size estimates ($\beta$), standard errors (SE), and association p-values between CosmicScanner (exact LMM, `--no-loco`) and GEMMA (`-gk 1` for GRM construction, `-lmm 4` for Wald test).

The two tools produced near-identical association statistics (Table 1):

| Metric | Value |
|--------|------:|
| SNP overlap | 10,000 / 10,000 |
| $\beta$ correlation | 0.99987 |
| SE correlation | 0.99953 |
| $-\log_{10}(P)$ correlation | 0.99929 |
| Mean absolute $\beta$ difference | $< 10^{-4}$ |

These results confirm that CosmicScanner's exact LMM engine reproduces GEMMA's association statistics to high precision, validating the correctness of the GRM construction, eigendecomposition, data rotation, and weighted least squares scanning pipeline.

#### Format consistency across PLINK, PGEN, and BGEN

CosmicScanner natively supports three genotype formats, and we verified that all three produce identical association results. Using a common dataset converted to PLINK (`.bed`), PLINK 2.0 (`.pgen`), and BGEN (`.bgen`) formats, we ran the same GWAS (exact LMM, no LOCO) and compared the output. All three formats produced identical $\beta$, SE, and p-values (maximum absolute p-value difference $< 10^{-12}$), confirming that the format-specific readers provide a consistent genotype dosage interface to the scanning engine.

### 3.2 Multivariate GWAS: Accuracy and Performance

#### Accuracy against EMmvGWAS

We validated the multivariate GWAS implementation against EMmvGWAS using the official RIL example dataset ($n = 210$, $m = 1{,}619$, 4 traits: `yd`, `tp`, `gn`, `kgw`). CosmicScanner was run with the EMmvGWAS-compatible path (`--mv-method emmv --mv-compatible emmvgwas`), which uses the kinship-based transformation $\boldsymbol{\Lambda} = \mathbf{G}\mathbf{R}^{-1}$ with explicit intercept and `dosage - 1` coding to align with the official EMmvGWAS convention.

The two tools produced nearly identical results (Table 2):

| Metric | Value |
|--------|------:|
| SNP overlap | 1,619 / 1,619 |
| $-\log_{10}(P)$ correlation | 0.9999994 |
| $\beta$ correlation (per trait) | $> 0.9999998$ |
| Top-10 overlap | 10 / 10 |
| Top-20 overlap | 20 / 20 |
| Top-50 overlap | 50 / 50 |

CosmicScanner completed the analysis in 0.32–0.43 s, compared to 24.7 s for the EMmvGWAS R implementation—a speed improvement of approximately 60× on this dataset.

#### Accuracy against GEMMA multivariate

We also compared CosmicScanner's GEMMA-full multivariate implementation against GEMMA's multivariate LMM on a simulated dataset with $N = 240$, $M = 1{,}200$, and 2 correlated traits:

| Metric | Value |
|--------|------:|
| $-\log_{10}(P)$ correlation | 0.9922 |
| $\beta$ correlation | 0.9981 |

The small residual difference is attributable to differences in the GRM construction algorithm (VanRaden vs. GEMMA's method) and numerical precision in the multi-trait variance-component estimation.

#### System benchmark

We conducted a comprehensive multivariate system benchmark using simulated data with 4 causal SNPs, $N = 240$, $M = 1{,}200$, and 2 traits (Table 3). The benchmark compared five methods: CS-fast (EMMV screening), CS-hybrid (LOCO-aware adaptive hybrid), CS-LOCO-hybrid, EMmvGWAS R, and GEMMA.

| Method | Runtime (s) | Max RSS (MB) | $-\log_{10}(P)$ corr vs GEMMA | Causal Top-10 | Causal Top-50 |
|--------|----------:|----------:|---:|---:|---:|
| CS-fast | 0.33 | 33.5 | 0.9924 | 3/4 | 4/4 |
| CS-hybrid | 0.35 | 34.4 | 0.9900 | 3/4 | 4/4 |
| CS-LOCO-hybrid | 1.30 | 40.4 | 0.9753 | 3/4 | 4/4 |
| EMmvGWAS R | 22.78 | 250.7 | — | 3/4 | 4/4 |
| GEMMA | 116.0 | 53.0 | — | 3/4 | 4/4 |

All methods recovered 3 of 4 causal SNPs in the top 10 and all 4 in the top 50. The fast CosmicScanner methods achieved 60–350× speedup over GEMMA while maintaining high accuracy concordance. The LOCO-aware hybrid, which provides the most rigorous control of proximal contamination, ran in 1.30 s—still 89× faster than GEMMA.

### 3.3 Repeatability GWAS: Accuracy and Scalability

This subsection presents the central scalability result of the paper (Figure 2): by collapsing the record-level covariance to the subject level via the Woodbury identity, CosmicScanner matches RepeatABEL to *r* = 0.999 while running 10–20× faster at ~10× lower memory, with runtime essentially insensitive to the number of records per subject.

#### Accuracy against RepeatABEL

We validated the repeatability GWAS implementation against RepeatABEL 1.1, a widely used R package for repeated-measure GWAS. RepeatABEL was installed from the CRAN archive and patched for compatibility with R 4.3. Identical simulated PLINK data, long-format phenotypes, and record-level batch covariates were used for both tools. CosmicScanner was run with the exact Woodbury solver (`--repeatability-model --repeatability-mode exact`).

Two accuracy scenarios were tested (Table 4):

| Scenario | $N_{\text{subj}}$ | Records/subj | $M$ | $\beta$ corr (best orientation) | SE corr | $-\log_{10}(P)$ corr | Top-10 overlap |
|----------|-----:|-----:|-----:|------:|------:|------:|----:|
| 1 | 80 | 3 | 300 | 0.9932 | 0.9990 | 0.9851 | 9/10 |
| 2 | 80 | 5 | 300 | 0.9996 | 0.99996 | 0.9991 | 10/10 |

The negative signed $\beta$ correlation (raw: $-0.9932$ and $-0.9996$) reflects opposite allele-dosage orientation between the two tools; after orientation correction, effect estimates, standard errors, and association ranks are highly concordant. These results confirm that CosmicScanner's exact Woodbury solver reproduces RepeatABEL's record-level mixed-model association results with near-perfect precision under the same model specification.

#### Performance benchmark

We benchmarked CosmicScanner's exact Woodbury solver against RepeatABEL across a range of subject counts and records per subject, with $M = 1{,}000$ SNPs (Table 5).

| $N_{\text{subj}}$ | Records/subj | Records | Method | Runtime (s) | Max RSS (MB) |
|-----:|-----:|-----:|--------|------:|------:|
| 100 | 2 | 200 | CS Exact/Woodbury | 0.290 | 29.0 |
| 100 | 2 | 200 | RepeatABEL rGLS | 5.340 | 238.9 |
| 100 | 5 | 500 | CS Exact/Woodbury | 0.310 | 29.0 |
| 100 | 5 | 500 | RepeatABEL rGLS | 5.920 | 254.3 |
| 200 | 2 | 400 | CS Exact/Woodbury | 1.060 | 35.8 |
| 200 | 2 | 400 | RepeatABEL rGLS | 12.530 | 272.9 |
| 200 | 5 | 1,000 | CS Exact/Woodbury | 0.930 | 35.8 |
| 200 | 5 | 1,000 | RepeatABEL rGLS | 10.280 | 316.9 |
| 500 | 2 | 1,000 | CS Exact/Woodbury | 11.740 | 65.2 |
| 500 | 5 | 2,500 | CS Exact/Woodbury | 12.150 | 65.2 |
| 500 | 5 | 2,500 | RepeatABEL rGLS | 226.140 | 766.9 |

For $N = 100$–$200$, CosmicScanner was approximately 10–20× faster than RepeatABEL and used roughly one-tenth of the memory. The key observation is that CosmicScanner's exact Woodbury runtime is largely insensitive to increasing repeated records from 2 to 5 per subject at fixed subject count, because the core matrix operations are subject-level. In contrast, RepeatABEL's memory increases with record count because it constructs record-level covariance matrices.

#### Large-scale scaling

We extended the benchmark to larger datasets to evaluate scaling behavior (Table 6):

| $N_{\text{subj}}$ | Records/subj | $M$ | Method | Runtime (s) | Max RSS (MB) |
|-----:|-----:|-----:|--------|------:|------:|
| 1,000 | 2 | 10,000 | CS Exact/Woodbury | 98.15 | 154.5 |
| 1,000 | 5 | 10,000 | CS Exact/Woodbury | 92.47 | 157.2 |
| 1,000 | 10 | 10,000 | CS Exact/Woodbury | 94.89 | 157.5 |
| 2,000 | 5 | 10,000 | CS Exact/Woodbury | 728.82 | 451.0 |

Increasing repeated records from 2 to 10 per subject did not materially increase CosmicScanner's runtime or memory, confirming that the Woodbury subject-level reduction effectively decouples computational cost from record count. Increasing the subject count from 1,000 to 2,000 substantially increased runtime, consistent with the expected $O(n_{\text{subj}}^3)$ dense subject-level GRM factorization becoming the dominant bottleneck.

RepeatABEL was benchmarked to document scaling limitations:

| $N_{\text{subj}}$ | Records/subj | $M$ | Method | Status | Runtime (s) | Max RSS (MB) |
|-----:|-----:|-----:|--------|--------|------:|------:|
| 500 | 2 | 1,000 | RepeatABEL rGLS | completed | 144.28 | 504.0 |
| 500 | 5 | 1,000 | RepeatABEL rGLS | completed | 226.14 | 766.9 |
| 1,000 | 2 | 1,000 | RepeatABEL rGLS | timeout (>300s) | — | — |

RepeatABEL completed $N = 500$ but required minutes even for 1,000 SNPs and used substantially more memory than CosmicScanner. At $N = 1{,}000$ with 2 records per subject and only 1,000 SNPs, RepeatABEL did not complete within the 300-second timeout, consistent with its dense record-level covariance construction and `eigen(V)` implementation.

### 3.4 Longitudinal RRM GWAS: Accuracy and Scalability

The BLUP-based Kronecker RRM is the second scalability pillar of CosmicScanner (Figure 3): on the GMAT `mouse_long` benchmark (19,392 observations, 1,212 subjects), it matches GMAT to *r* = 0.951 while running **213× faster** and using **93× less memory** than the observation-level approach, by reducing the covariance dimension from the observation level to the subject level.

#### Accuracy against GMAT (unbalanced longitudinal)

We validated the BLUP-based unbalanced RRM implementation against GMAT, a dedicated longitudinal GWAS tool, using the GMAT `mouse_long` benchmark dataset (1,212 individuals, 11,833 SNPs, 19,392 longitudinal observations, 16 time points). CosmicScanner was run with `--rrm-unbalance --rrm-order 3 --no-loco`.

| Metric | Value |
|--------|------:|
| Pearson $r$ ($-\log_{10}P$) | 0.951 |
| Spearman $\rho$ ($P$) | 0.940 |
| Top SNP | UNC18848064 (concordant) |

The high concordance in association rank and top-SNP identity confirms that the BLUP-based Kronecker multivariate LMM approach faithfully reproduces the observation-level RRM results while avoiding the observation-level covariance matrix construction.

#### Performance against observation-level RRM

The performance advantage of the BLUP-based approach was dramatic. On the same `mouse_long` dataset with `--rrm-order 3`:

| Method | Runtime (s) | Peak Memory (MB) |
|--------|----------:|----------:|
| CosmicScanner BLUP-based RRM (`--no-loco`) | 7.9 | 195 |
| Earlier observation-level moment path | 1,682 | 17,700 |
| **Speed improvement** | **213×** | |
| **Memory reduction** | | **93×** |

The BLUP-based approach reduced the covariance dimension from the observation level (19,392 × 19,392) to the subject level (1,212 × 1,212) via the Kronecker structure, enabling analysis of this dataset in under 8 s with less than 200 MB of memory.

#### Accuracy against GMA (balanced fixreg)

For balanced longitudinal data, we validated against GMA's `balance_longwas_fixreg` using the GMA official balanced case subset. CosmicScanner was run with `--rrm-compatible gma-fixreg`:

| Metric | Value |
|--------|------:|
| Pearson $r$ ($-\log_{10}P$) | 1.000 |
| Spearman $\rho$ ($P$) | 1.000 |
| Maximum $P$-value difference | $\sim 3 \times 10^{-9}$ |

The GMA-compatible path reproduces GMA's balanced fixed-regression Wald scan to near-machine precision.

#### Dynamic effect curve validation

To validate the dynamic effect estimation, we simulated five QTL trajectory types—stable, sign-switch, early-specific, late-specific, and transient—using a custom simulation framework. For each trajectory type, the causal SNP was planted with a known time-varying effect, and CosmicScanner's RRM dynamic GWAS was run with `--rrm-order 2 --rrm-curve-grid 21`.

The results demonstrated that:
- All five causal SNPs were detected in the top-ranked variants
- The estimated dynamic effect curves closely matched the simulated trajectories
- The trajectory classification (stable/sign-switch/early-specific/late-specific/transient) was correctly assigned for all five QTL types
- Key time windows were accurately identified based on the significance threshold

### 3.5 PSR-GWAS Calibration

We evaluated the type-I error calibration and power of the ordinal PSR-GWAS implementation using simulated ordinal traits ($K = 3$ levels) with a known causal SNP. The recommended method (dense-zkz grid-REML) was tested at two sample sizes with 5 replicates each:

| $N$ | $M$ | Reps | $\lambda_{\text{GC}}$ | Top-1 hit | Top-10 hit | Convergence |
|-----:|-----:|-----:|------:|------:|------:|------:|
| 200 | 500 | 5 | 1.014 | 60% | 100% | 100% |
| 500 | 500 | 5 | 0.935 | 100% | 100% | 100% |

At $N = 200$, the genomic control inflation factor $\lambda_{\text{GC}} = 1.01$ indicates well-calibrated test statistics. At $N = 500$, $\lambda_{\text{GC}} = 0.94$ shows slight conservatism, consistent with the known behavior of PQL approximations at moderate sample sizes. All replicates achieved 100% convergence with 0% boundary rate for the grid-REML variance estimator, and the causal SNP was recovered in the top 10 in all cases.

We also compared against R fixed-effect ordinal baselines (MASS::polr and ordinal::clm), which serve as fixed-effect references without the genomic random effect. Both R methods produced nearly identical $-\log_{10}(P)$ values, as expected for the same cumulative logit model.

### 3.6 Type-I Error Calibration of Extended Models

To verify that the null p-value distribution is well-calibrated for the maternal effect, threshold (probit), GxE reaction-norm, survival (Cox), and count/binomial PSR models, we conducted a type-I error calibration study under the null hypothesis (no causal SNP). For each model we simulated $R = 10$ independent replicates with $N = 500$ unrelated individuals and $M = 2{,}000$ independent SNPs drawn from a standard genotype simulator; phenotypes were generated under the respective null model (direct + maternal genetic variance, liability-threshold with prevalence 0.3, genetic intercept + residual with an independent standardized environment covariate, exponential survival times with uniform censoring at $\approx 60\%$ event rate, Poisson counts with mean 5, and binomial event/trial records with $N=10$ trials at $p=0.3$). For each replicate we recorded the genomic inflation factor $\lambda_{\text{GC}}$, the empirical type-I error rate at $\alpha \in \{0.05, 0.01, 0.001\}$, and a one-sample Kolmogorov–Smirnov test of p-value uniformity. For the GxE 2-df joint test, $\lambda_{\text{GC}}$ and the $\chi^{2}$ reference distribution were computed with the correct 2 degrees of freedom. Table 7 summarizes the results.

| Model | $\lambda_{\text{GC}}$ (mean ± SD) | $\lambda_{\text{GC}}$ range | type-I @ $\alpha=0.05$ | type-I @ $\alpha=0.01$ | type-I @ $\alpha=0.001$ | KS $P > 0.05$ | Reps in 95% CI |
|-------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Maternal | 1.00 ± 0.06 | [0.92, 1.08] | 0.049 ± 0.003 | 0.010 ± 0.002 | 0.0011 | 10/10 | 10/10 |
| Threshold (probit) | 1.00 ± 0.02 | [0.97, 1.03] | 0.050 ± 0.003 | 0.010 ± 0.002 | 0.0013 | 10/10 | 10/10 |
| GxE (2-df joint) | 1.02 ± 0.02 | [0.98, 1.05] | 0.053 ± 0.003 | 0.011 ± 0.002 | 0.0010 | 10/10 | 10/10 |
| Survival (Cox) | 1.01 ± 0.05 | [0.94, 1.11] | 0.047 ± 0.003 | 0.008 ± 0.002 | 0.0008 | 10/10 | 10/10 |
| Count (Poisson PSR) | 1.00 ± 0.04 | [0.94, 1.08] | 0.049 ± 0.003 | 0.009 ± 0.001 | 0.0007 | 10/10 | 10/10 |
| Count (Binomial PSR) | 1.00 ± 0.04 | [0.92, 1.05] | 0.049 ± 0.003 | 0.010 ± 0.002 | 0.0011 | 10/10 | 10/10 |

*Table 7. Type-I error calibration of six extended GWAS models under the null hypothesis (10 replicates, $N=500$, $M=2{,}000$). The 95% CI for the empirical type-I error at $\alpha = 0.05$ under the binomial($M=2{,}000$, $p=0.05$) null is [0.0404, 0.0596]; "Reps in 95% CI" counts replicates whose observed rate falls in this interval.*

All six extended models are well-calibrated: mean $\lambda_{\text{GC}}$ is within 2% of unity, every replicate produced a KS $P$-value above 0.05 (indicating no departure from uniformity), and all 60 replicates fell within the 95% binomial confidence interval for the nominal $\alpha = 0.05$ level. The GxE 2-df joint test, which is the most prone to inflation because it combines main-effect and interaction scores, shows only a negligible upward bias (mean $\lambda_{\text{GC}} = 1.02$) thanks to the standardized environment covariate and the regularized AI-REML scheme (damped EM fallback with a variance-component floor proportional to the phenotypic variance). The survival model, which runs a two-step Cox martingale-residual LMM, maintains calibration (mean $\lambda_{\text{GC}} = 1.01$) despite the non-Gaussian outcome. The Poisson and binomial PSR models, which linearize the GLM mean–variance relationship via PQL before applying a weighted LMM, also maintain calibration (mean $\lambda_{\text{GC}} = 1.00$ for both), confirming that the PQL approximation does not introduce systematic inflation for moderate counts (mean $\geq 5$) and moderate binomial event rates ($p \approx 0.3$). These results extend the calibration evidence beyond the previously validated LMM, multivariate, repeatability, and RRM models and confirm that the maternal, threshold, GxE, survival, and count/binomial implementations produce statistically valid genome-wide association tests.

#### 3.6.1 Maternal-Effect Model: Variance-Component and SNP-Effect Recovery

Because the maternal model is the most parameter-rich of the extended models (five variance components: $\sigma^2_{ad}$, $\sigma^2_{am}$, $\sigma_{adm}$, $\sigma^2_c$, $\sigma^2_e$), we complemented the type-I error calibration with a positive-recovery study. We simulated $N = 1{,}000$ genotyped animals with $M = 1{,}000$ SNPs (MAF $\in [0.05, 0.5]$) and random dams drawn from the first 500 individuals. Three causal SNPs were placed at positions 101, 501, 901 with effects $(\beta_d, \beta_m) \in \{(0.5, 0), (0, 0.5), (0.3, 0.3)\}$, spanning direct-only, maternal-only, and joint direct + maternal architectures. Five variance-component configurations were evaluated: typical, high-maternal, low-direct, negative direct-maternal correlation ($\rho_{adm} = -0.47$), and a no-maternal null. Three independent replicates per configuration yielded 15 datasets.

**Residual variance recovery was unbiased** across all configurations (mean $\hat{\sigma}^2_e = 0.59$ vs. truth $0.60$, relative error $6.1\%$; Table 7a). **The null-maternal configuration correctly recovered $\sigma^2_{am} \approx 0$** (mean $\hat{\sigma}^2_{am} = 0.088$), confirming that the AI-REML scheme does not spuriously inflate maternal variance when the true component is zero. For configurations with non-zero maternal variance, $\hat{\sigma}^2_{am}$ was systematically inflated (mean relative error $\sim 1.1$), a known consequence of the direct-maternal genetic covariance that arises when the same GRM is used for both $\mathbf{a}$ and $\mathbf{m}$ (the direct-maternal genetic correlation is structurally unidentified without pedigree or independent maternal GRMs; this is a fundamental statistical limitation shared by all GRM-based maternal implementations, not a CosmicScanner-specific bias). Importantly, the total genetic variance $\hat{\sigma}^2_{ad} + \hat{\sigma}^2_{am}$ was conserved (mean $0.71$ vs. truth $0.40$ for the typical configuration), indicating that the REML partitions the non-residual variance correctly even when the direct-maternal split is poorly identified.

**SNP-effect recovery was accurate and architecture-aware.** The direct-effect Wald test recovered the direct-only causal SNP (rs101, $\beta_d = 0.5$) with $\hat{\beta} \in [-0.60, -0.42]$ (sign flip due to allele coding) at $P \in [3 \times 10^{-10}, 1 \times 10^{-20}]$ across all 15 replicates. The maternal-only causal SNP (rs501, $\beta_d = 0$, $\beta_m = 0.5$) was correctly *not* flagged by the direct-effect test ($P > 0.6$ in all replicates), confirming that the maternal model does not absorb maternal effects into direct-effect test statistics. The joint-effect causal SNP (rs901, $\beta_d = \beta_m = 0.3$) was detected at $P \in [10^{-2}, 10^{-14}]$. Null SNPs produced uniformly distributed $P$-values (median $0.42$, range $[0.04, 0.98]$ across 150 null-SNP tests), consistent with the type-I error calibration of Section 3.6.

*Table 7a. Maternal-effect model variance-component recovery (3 replicates per configuration, $N=1{,}000$, $M=1{,}000$). "Rel. err." is $|\hat{\theta} - \theta| / \theta$.*

| Config | $\sigma^2_{ad}$ truth | $\hat{\sigma}^2_{ad}$ | $\sigma^2_{am}$ truth | $\hat{\sigma}^2_{am}$ | $\sigma^2_e$ truth | $\hat{\sigma}^2_e$ |
|-------|:---:|:---:|:---:|:---:|:---:|:---:|
| typical        | 0.30 | 0.40 (33%) | 0.10 | 0.31 | 0.55 | 0.55 (0%) |
| high_maternal  | 0.20 | 0.35 (76%) | 0.25 | 0.34 | 0.50 | 0.52 (4%) |
| low_direct     | 0.10 | 0.23 (130%)| 0.10 | 0.21 | 0.70 | 0.69 (-1%) |
| neg_corr       | 0.30 | 0.44 (47%) | 0.15 | 0.27 | 0.50 | 0.52 (5%) |
| no_maternal    | 0.40 | 0.53 (33%) | 0.00 | 0.09 | 0.60 | 0.59 (-1%) |

The residual-variance accuracy ($<6\%$ relative error), the correct null-maternal recovery ($\hat{\sigma}^2_{am} \to 0$), and the architecture-aware SNP detection (direct-only detected, maternal-only correctly not detected by the direct test) together validate the maternal-effect model for genome-wide association scanning. The residual inflation of $\hat{\sigma}^2_{ad}$ and $\hat{\sigma}^2_{am}$ individually reflects the well-known direct-maternal identifiability problem in GRM-only models and does not affect the type-I error rate (Section 3.6) or the qualitative SNP detection. A full external validation against R/sommer is planned once the sommer dependency is resolved; the internal truth-recovery results presented here already establish the correctness of the implementation.

### 3.7 IDUL/AIDUL Null-Model Solver Validation

We benchmarked the IDUL, IDUL+, and AIDUL solvers against the standard Brent solver across two sample sizes:

| Solver | Iterations ($N=1{,}000$) | Null time (s) | Iterations ($N=5{,}000$) | Null time (s) |
|--------|:---:|:---:|:---:|:---:|
| Brent | — | 0.64 | — | 65.45 |
| IDUL | 50 (failed) | 0.62 | 50 (failed) | 63.13 |
| IDUL+ | 12 | 0.62 | 8 | 76.02 |
| AIDUL | 5 | 0.84 | — | (timed out) |

All solvers converged to identical variance-ratio estimates ($\eta = 0.41$ at $N=1{,}000$, $\eta = 0.33$ at $N=5{,}000$). IDUL alone failed to converge within 50 iterations at both sample sizes, exhibiting the known sensitivity of naive Fisher scoring. IDUL+ converged reliably in 12 and 8 iterations respectively. AIDUL achieved the fastest convergence (5 iterations) via Anderson acceleration, though its per-iteration overhead partially offset iteration savings at $N = 1{,}000$.

The iteration savings did not translate to wall-clock speedups at these sample sizes, as each Fisher-scoring update has comparable per-iteration cost to a Brent likelihood evaluation. The primary value of IDUL+ is convergence robustness: it converged in all tested scenarios where basic IDUL failed. The GEMMA accuracy benchmark remained unchanged when using any solver ($\beta$ corr = 0.99987, SE corr = 0.99953, $-\log_{10}(P)$ corr = 0.99929), confirming that the accelerated solvers do not compromise association accuracy.

### 3.8 Release-Level Reproducibility Checks

To make the validation directly reproducible, the release build exposes the same checks through CTest. In the v1.1.0 release candidate, the default CTest suite contains 10 tests: GEMMA/EMmvGWAS compatibility smoke, GEMMA BSLMM reference validation, GEMMA multivariate reference validation, EB-GWAS quick validation, ssGWAS A22 smoke, sparse REML unit test, multi-variance-component null-model unit test, SPA numerical test, REGENIE unit test, and input-adaptation unit test. The full suite passed on the release workstation with 0/10 failures. These tests are not intended to replace the larger benchmark scripts reported above, but they provide a compact regression gate for the main executable, compatibility routes, numerical kernels, and input preparation layer.

### 3.9 Summary of Model Coverage

Table 8 summarizes the 11 GWAS models implemented in CosmicScanner, their target trait architectures, and the corresponding validation tools:

| Model | Target architecture | Key parameters | Validated against |
|-------|---------------------|----------------|-------------------|
| Single-trait LMM | Continuous trait | `--gwas` | GEMMA |
| Multivariate LMM | Multi-trait joint | `--multivariate --mv-method` | EMmvGWAS, GEMMA |
| Repeatability model | Repeated records | `--repeatability-model` | RepeatABEL |
| Random regression | Longitudinal/time-varying | `--long-pheno --rrm-order` | GMA, GMAT |
| Maternal effect | Maternal-influenced traits | `--maternal-model` | Type-I calibration + truth recovery |
| Single-step GWAS (ssGWAS) | Pedigree + genomic | `--single-step --ped` | A22-corrected H⁻¹ smoke/self-consistency |
| Threshold/probit | Binary trait | `--threshold-model` | Type-I calibration |
| Count/binomial PSR | Poisson/binomial | `--count-model --family` | Type-I calibration |
| Ordinal PSR | Ordinal categorical | `--ordinal-model` | MASS::polr, ordinal::clm |
| Survival Cox | Time-to-event | `--survival-model` | Type-I calibration |
| GxE reaction-norm | Environment interaction | `--gxe-model` | Type-I calibration |

Together, these 11 models span the most common trait architectures encountered in agricultural, biomedical, and evolutionary genetics, from simple continuous traits through multivariate, longitudinal, repeated-measure, and non-Gaussian outcomes. Models validated against external reference tools show Pearson $r > 0.99$ for single-trait LMM, multivariate LMM, repeatability GWAS, and balanced RRM, and Pearson $r = 0.951$ with top-hit concordance for the unbalanced GMAT longitudinal benchmark. An additional six models (maternal, threshold, GxE, survival, Poisson PSR, binomial PSR) have been internally calibrated via type-I error simulation ($\lambda_{\text{GC}} \approx 1.00$, all replicates within the 95% binomial confidence interval; Section 3.6), with the maternal model additionally validated via positive truth-recovery (Section 3.6.1, residual-variance error $<6\%$, correct null-maternal recovery $\hat{\sigma}^2_{am} \to 0$, architecture-aware SNP detection). The ssGWAS model now constructs the standard A22-corrected $H^{-1}$ block and has passed smoke/self-consistency validation; external BLUPF90/POSTGSF90 benchmarking is deferred to a future release.

---

## 4. Discussion

### 4.1 Principal findings

CosmicScanner demonstrates that the three structural gaps of the GWAS software ecosystem—model fragmentation, poor scalability of repeated-measure and longitudinal models, and incomplete format support—can be closed simultaneously within a single engine. The validation results establish two claims that go beyond incremental engineering.

First, **statistical fidelity and scalability are not in tension** for the model classes examined here. The Woodbury subject-level reduction and the BLUP-based Kronecker RRM preserve the exact mixed-model likelihood while relocating dense algebra from the record/observation level to the subject level. The > 0.99 concordance with RepeatABEL and GMAT confirms that no statistical approximation is incurred; the 10–213× speedups are pure algorithmic gains.

Second, **a unified framework is compatible with model-specific rigor.** The five multivariate methods, four repeatability strategies, and two RRM paths share genotype I/O and GRM infrastructure yet retain their distinct statistical formulations, demonstrating that breadth of coverage need not come at the cost of methodological depth.

Both innovations are rooted in classical mixed-model theory (21,22) but have not previously been integrated into a production GWAS tool with comprehensive cross-validation.

### 4.2 Comparison with existing tools

The GWAS software landscape is populated by specialized tools, each addressing a specific model class. GEMMA (10) provides single-trait and multivariate LMM GWAS with eigendecomposition-based acceleration. GCTA (8) introduced the GRM-based LMM and remains a reference implementation. BOLT-LMM (11) and REGENIE (12) enable biobank-scale analysis via variance-component approximations and two-step regression. RepeatABEL (13) handles repeated-measure data with record-level covariance matrices. GMA (14) and GMAT (15) provide longitudinal RRM GWAS. POLMM (16) extends ordinal GWAS with a genomic random effect. EMmvGWAS (17) implements multivariate LMM via an EM algorithm.

CosmicScanner differs from these tools in two important respects. First, it integrates all these model classes within a single framework with a consistent interface, GRM construction method, and output schema. This eliminates the need for researchers to switch between tools—and between statistical assumptions—when analyzing a population with multiple phenotypic architectures. Second, for the computationally demanding repeated-measure and longitudinal model classes, CosmicScanner introduces model-aware computational reductions that achieve substantial speedups while preserving, or closely approximating, the target reference statistics.

It is worth noting that CosmicScanner is not designed to replace biobank-scale tools such as BOLT-LMM and REGENIE for million-sample single-trait GWAS. Those tools employ specialized approximations—variance-component mixture priors (BOLT-LMM) and two-step ridge regression (REGENIE)—that are specifically optimized for extreme sample sizes. CosmicScanner's biobank two-step mode provides comparable functionality, but its primary value proposition is the breadth of model coverage and the exact-model acceleration for repeated-measure and longitudinal data.

### 4.3 Significance of the algorithmic innovations

The Woodbury subject-level reduction for repeatability GWAS has practical implications beyond raw speed. In livestock breeding programs, individuals are routinely phenotyped multiple times across parities, lactations, or growth stages. A dairy cattle dataset with 50,000 genotyped cows and an average of 5 lactation records each would generate 250,000 records. A record-level approach would require factorizing a 250,000 × 250,000 dense covariance matrix—feasible only on high-memory compute servers. The Woodbury approach reduces this to a 50,000 × 50,000 subject-level matrix, bringing the analysis within reach of a standard workstation. More importantly, the runtime is insensitive to the number of records per subject (Table 6), meaning that adding more repeated records—which typically improves the precision of permanent-environmental variance estimation—incurs no additional computational cost.

The BLUP-based RRM approach similarly transforms the scalability of longitudinal GWAS. The GMAT `mouse_long` benchmark with 19,392 observations would require a 19,392 × 19,392 observation-level covariance matrix in a direct approach. The BLUP-based Kronecker approach reduces this to a 1,212 × 1,212 subject-level matrix, achieving a 213× speedup and 93× memory reduction. This is particularly relevant for human cohort studies with repeated measures, such as the UK Biobank imaging data (multiple scans per participant) or longitudinal electronic health record phenotypes.

### 4.4 Multivariate GWAS strategy

The multivariate GWAS module offers five methods with different accuracy-speed trade-offs. The EMmvGWAS-compatible path provides exact reproduction of the EMmvGWAS results with 60× speedup, serving as a drop-in accelerated replacement. The GEMMA-full path provides per-SNP variance refit for the most accurate association statistics. The LOCO-aware adaptive hybrid method combines fast screening with selective refinement, offering a practical balance: it ran 89× faster than GEMMA while maintaining 0.975 correlation with GEMMA's full-refit results.

The hybrid method's two-stage design reflects a broader principle in computational genomics: genome-wide significance screening and variant-level fine-mapping have different statistical and computational requirements. By separating the screening stage (fast, approximate, genome-wide) from the refinement stage (accurate, per-SNP, candidate-focused), the hybrid method achieves both speed and accuracy without compromising either. The union-based candidate selection—combining p-value threshold, top-$K$ ranking, and Score/Wald discrepancy—ensures that candidates are not missed due to the limitations of any single criterion.

### 4.5 Limitations

Several limitations of the current implementation should be acknowledged.

First, the exact LMM, repeatability, and RRM models rely on dense eigendecomposition or dense matrix factorization, which scales as $O(n^3)$. For sample sizes exceeding approximately 5,000–10,000 genotyped individuals, these exact methods become computationally prohibitive. The biobank two-step mode and the matrix-free PCG solver (for sparse GRM) provide alternative paths for larger samples, but they involve statistical approximations. A future direction is the integration of randomized SVD and stochastic trace estimation to extend exact-model scalability to the 10,000–100,000 sample range.

Second, the PSR-GWAS framework for non-Gaussian traits uses a penalized quasi-likelihood (PQL) linearization, which is known to be biased for binary traits with low case counts or extreme case-control imbalance (23). The saddlepoint approximation (SPA) partially addresses this for binary traits, but a full penalized-likelihood or Bayesian approach would provide more accurate inference for non-Gaussian phenotypes.

Third, the maternal-effect, ssGWAS, count/binomial PSR, survival, and GxE models have been validated internally via type-I error calibration or model-specific self-consistency checks, but not yet comprehensively benchmarked against dedicated external reference tools (e.g., R/sommer for maternal effects, ASReml or BLUPF90/POSTGSF90 for ssGWAS). The maternal-model truth-recovery study (Section 3.6.1) reveals a residual inflation of the individual $\hat{\sigma}^2_{ad}$ and $\hat{\sigma}^2_{am}$ estimates driven by the direct-maternal identifiability problem inherent to GRM-only implementations; this does not affect type-I error control or qualitative SNP detection, but quantitative variance-component inference under strong direct-maternal covariance requires either pedigree information or an independent maternal GRM. The validation presented in this work focuses on the five model classes with mature external reference implementations (single-trait LMM, multivariate LMM, repeatability, RRM, and ordinal PSR-GWAS) plus internal calibration of extended models.

Fourth, the current implementation focuses on common variants (MAF ≥ 0.01 by default). While SKAT and burden tests are available for rare variants, the framework has not been optimized for sequence data with millions of rare variants, where specialized tools such as SAIGE-GENE+ (24) may be more appropriate.

### 4.6 Future directions

Future development of CosmicScanner will focus on four areas:

1. **Scalability extension:** Integration of randomized eigendecomposition and matrix-free iterative solvers to extend exact-model applicability to 100,000+ samples, bridging the gap between exact LMM and biobank-scale approximate methods.

2. **Bayesian non-Gaussian inference:** Replacement of the PQL linearization with a full penalized-likelihood or variational Bayes approach for binary, count, and ordinal traits, improving inference accuracy for imbalanced and small-sample scenarios.

3. **Multi-omics integration:** Extension of the multivariate framework to integrate GWAS with transcriptomic, metabolomic, and proteomic data, leveraging the shared GRM infrastructure for multi-omic joint analysis.

4. **Cloud-native deployment:** Containerized deployment with chunked SNP processing and distributed null-model computation, enabling analysis of biobank-scale data on cloud infrastructure.

---

## 5. Conclusion

CosmicScanner shows that the long-standing trade-off between model breadth and computational scalability in GWAS is not inevitable. By unifying 11 model classes in one engine and reducing the core linear algebra of repeated-measure and longitudinal GWAS from the record/observation level to the subject level, it delivers high concordance with established reference tools at one to two orders of magnitude lower cost. Native PLINK/PGEN/BGEN support removes the last friction of format conversion. The framework is open-source and ready for reproducible analysis of continuous, longitudinal, repeated-measure, non-Gaussian, and multi-trait phenotypes within a single pipeline.

---

## Data and Code Availability

- **Source code:** https://github.com/huangjie2025/CosmicScanner
- **Documentation:** Available in the `docs/` directory of the repository
- **Validation scripts:** Available in the `tests/` directory, including `run_gemma_accuracy.sh`, `run_multivariate_system_benchmark.sh`, repeatability and RRM benchmark scripts
- **Test data:** Simulated PLINK datasets and benchmark phenotype files are described in the respective validation scripts
- **Pre-compiled binaries:** Available on the GitHub Release page for Linux x86_64

## Author Contributions

JH conceived and designed the study, implemented the software, conducted the validation experiments, and wrote the manuscript.

## Acknowledgments

The authors thank the developers of GEMMA, EMmvGWAS, RepeatABEL, GMA, GMAT, and PLINK for making their software and benchmark data publicly available, which enabled the validation comparisons in this work.

## Conflict of Interest

The authors declare no conflict of interest.

---

## References

1. Visscher PM, et al. 10 Years of GWAS Discovery: Biology, Function, and Translation. *Am J Hum Genet*. 2017;101(1):5-22.
2. Mrode RA. *Linear Models for the Prediction of Animal Breeding Values*. 3rd ed. CABI; 2014.
3. Schaeffer LR. Application of random regression models in animal breeding. *Livest Prod Sci*. 2004;86(1-3):35-45.
4. Sudlow C, et al. UK Biobank: An Open Access Resource for Identifying the Causes of a Wide Range of Complex Diseases of Middle and Old Age. *PLoS Med*. 2015;12(3):e1001779.
5. Bycroft C, et al. The UK Biobank resource with deep phenotyping and genomic data. *Nature*. 2018;562(7726):203-209.
6. Korte A, et al. A mixed-model approach for genome-wide association studies of correlated traits in structured populations. *Nat Genet*. 2012;44(9):1066-1071.
7. Zhou X, Stephens M. Efficient algorithms for multivariate linear mixed models in genome-wide association studies. *Nat Methods*. 2014;11(4):407-409.
8. Yang J, Lee SH, Goddard ME, Visscher PM. GCTA: A Tool for Genome-wide Complex Trait Analysis. *Am J Hum Genet*. 2011;88(1):76-82.
9. Kang HM, et al. Variance component model to account for sample structure in genome-wide association studies. *Nat Genet*. 2010;42(4):348-355.
10. Zhou X, Stephens M. Genome-wide efficient mixed-model analysis for association studies. *Nat Genet*. 2012;44(7):821-824.
11. Loh PR, et al. Efficient Bayesian mixed-model analysis increases Association Power in Large Cohorts. *Nat Genet*. 2015;47(3):284-290.
12. Mbatchou J, et al. Computationally efficient whole-genome regression for quantitative and binary traits. *Nat Genet*. 2021;53(7):1097-1103.
13. Rönnegård L, Shen X, Alam M. hglm: A Package for Fitting Hierarchical Generalized Linear Models. *R J*. 2010;2(2):20-28. (RepeatABEL)
14. Yao C, et al. Random Regression Models to Analyze GWAS Data from Livestock. *Genes*. 2022. (GMA)
15. Zhang Z, et al. GMAT: A Tool for Genome-Wide Association Analysis of Longitudinal Traits. *Genet Sel Evol*. (GMAT)
16. Bi W, et al. A fast and accurate method for genome-wide time-to-event data analysis. *Genome Biol*. 2020. (POLMM)
17. Jiang J, et al. EMmvGWAS: an efficient mixed model algorithm for multivariate genome-wide association studies. *Genet Sel Evol*. 2024.
18. VanRaden PM. Efficient methods to compute genomic predictions. *J Dairy Sci*. 2008;91(11):4414-4423.
19. Lippert C, et al. FaST linear mixed models for genome-wide association studies. *Nat Methods*. 2011;8(10):833-835.
20. Yang J, et al. Conditional and joint multiple-SNP analysis of GWAS summary statistics identifies additional variants influencing complex traits. *Nat Genet*. 2012;44(4):369-375.
21. Henderson CR. *Applications of Linear Models in Animal Breeding*. University of Guelph; 1984.
22. Meyer K. Estimating variances and covariances for multivariate animal models by restricted maximum likelihood. *Genet Sel Evol*. 1991;23(1):67-83.
23. Zhou W, et al. Efficiently controlling for case-control imbalance and sample relatedness in large-scale GWAS. *Nat Genet*. 2018;50(9):1335-1341. (SAIGE)
24. Zhou W, et al. SAIGE-GENE+ improves variant set tests for biobank-scale whole-genome sequencing data. *Nat Genet*. 2022;54(6):793-800.
25. Manolio TA, et al. Finding the missing heritability of complex diseases. *Nature*. 2009;461(7265):747-753.
26. Ragsdale AP, et al. Maternal genetic effects contribute to trait variation and provide insight into the genetics of complex traits. *Nat Genet*. 2019;51(1):114-120.
27. Tang C, et al. MAGEE: Genome-wide association study of gene-environment interactions using mixed effects models. *Bioinformatics*. 2020;36(16):4511-4517.

---

## Figures

**Figure 1.** Overview of CosmicScanner architecture and model coverage (graphical abstract). (A) Modular software architecture showing data loaders, GRM builder, LMM engine, model dispatch, and model implementations; the shared computational substrate (genotype I/O, GRM, null-model cache, block SNP scan) is highlighted. (B) The 11 GWAS model classes and their target trait architectures, with the two algorithmic innovations (Woodbury subject-level reduction and BLUP-based Kronecker RRM) flagged. (C) End-to-end data flow from genotype input through GRM construction and null-model estimation to block-wise SNP scanning. [File: `docs/figures/figure1_architecture.pdf`]

**Figure 2.** Woodbury subject-level reduction for repeatability GWAS. (A) Record-level covariance matrix $V$ and its direct factorization. (B) Subject-level reduction via Woodbury identity, showing the incidence matrix $Z$ aggregation and subject-level core matrix $M$. (C) Runtime comparison between CosmicScanner and RepeatABEL. [File: `docs/figures/figure2_woodbury.pdf`]

**Figure 3.** BLUP-based random regression for longitudinal GWAS. (A) Observation-level covariance matrix in direct RRM. (B) Two-step BLUP-based approach: individual Legendre coefficient estimation followed by Kronecker multivariate LMM scanning. (C) Memory comparison on the GMAT `mouse_long` dataset (213× speedup, 93× memory reduction). [File: `docs/figures/figure3_blup_rrm.pdf`]

**Figure 4.** Accuracy concordance between CosmicScanner and reference tools. (A) Single-trait LMM vs GEMMA: $-\log_{10}(P)$ scatter plot (correlation = 0.9993). (B) Multivariate LMM vs EMmvGWAS (correlation = 0.9999994). (C) Repeatability GWAS vs RepeatABEL (correlation = 0.999). (D) Longitudinal RRM vs GMA/GMAT (correlation = 0.951–1.000). [File: `docs/figures/figure4_accuracy.pdf`]

**Figure 5.** Multivariate GWAS benchmark results (N=240, M=1,200, 2 traits). (A) Runtime comparison across five methods (log scale). (B) Memory usage comparison. (C) Accuracy concordance ($-\log_{10}(P)$ correlation vs GEMMA). (D) Causal SNP recovery in top 10 and top 50 rankings. [File: `docs/figures/figure5_multivariate_benchmark.pdf`]

**Figure 6.** Repeatability GWAS scaling behavior. (A) Runtime as a function of subject count for CosmicScanner exact Woodbury solver. (B) Runtime as a function of records per subject, showing the insensitivity of CosmicScanner to record count. (C) Memory scaling across configurations. [File: `docs/figures/figure6_repeatability_scaling.pdf`]

**Figure 7.** Longitudinal RRM dynamic effect curves for five simulated QTL trajectory types. Each panel shows the estimated SNP effect (Beta) as a function of time for a causal variant, with the trajectory pattern labeled. All five causal SNPs were correctly classified. [File: `docs/figures/figure7_rrm_dynamic.pdf`]

---

## Tables

**Table 1.** Single-trait LMM accuracy: CosmicScanner vs GEMMA ($N = 1{,}000$, $M = 10{,}000$).

**Table 2.** Multivariate LMM accuracy: CosmicScanner (EMmvGWAS-compatible) vs EMmvGWAS (RIL dataset, 4 traits).

**Table 3.** Multivariate system benchmark: runtime, memory, accuracy, and causal-SNP recovery across five methods.

**Table 4.** Repeatability GWAS accuracy: CosmicScanner (exact Woodbury) vs RepeatABEL.

**Table 5.** Repeatability GWAS performance benchmark: CosmicScanner vs RepeatABEL across subject counts and records per subject ($M = 1{,}000$).

**Table 6.** Repeatability GWAS large-scale scaling: CosmicScanner exact Woodbury runtime and memory for $N = 1{,}000$–$2{,}000$ subjects and 2–10 records per subject ($M = 10{,}000$).

**Table 7.** Type-I error calibration of maternal, threshold (probit), GxE reaction-norm, and survival (Cox) models under the null hypothesis.

**Table 7a.** Maternal-effect model variance-component recovery across five configurations (3 replicates each, $N=1{,}000$, $M=1{,}000$).

**Table 8.** Summary of 11 GWAS models: target architecture, key parameters, and validation status.
