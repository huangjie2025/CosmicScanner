#!/usr/bin/env Rscript
# PSR-GWAS external method availability and ordinal R baseline.

options(warn=1)
pkg_available <- function(pkg) requireNamespace(pkg, quietly=TRUE)

args <- commandArgs(trailingOnly=TRUE)
mode <- if (length(args) >= 1) args[[1]] else "availability"

if (mode == "availability") {
  out <- if (length(args) >= 2) args[[2]] else "tests/results/psr_gwas_smoke/psr_method_availability.tsv"
  rows <- list(
    data.frame(trait_type="normal", method="lm", method_available=TRUE, error_message=""),
    data.frame(trait_type="binary", method="glm_binomial", method_available=TRUE, error_message=""),
    data.frame(trait_type="binary", method="logistf", method_available=pkg_available("logistf"), error_message=ifelse(pkg_available("logistf"), "", "R package logistf is not available")),
    data.frame(trait_type="binary", method="GMMAT", method_available=pkg_available("GMMAT"), error_message=ifelse(pkg_available("GMMAT"), "", "R package GMMAT is not available")),
    data.frame(trait_type="binomial", method="glm_cbind", method_available=TRUE, error_message=""),
    data.frame(trait_type="binomial", method="glmer", method_available=pkg_available("lme4"), error_message=ifelse(pkg_available("lme4"), "", "R package lme4 is not available")),
    data.frame(trait_type="binomial", method="glmmTMB", method_available=pkg_available("glmmTMB"), error_message=ifelse(pkg_available("glmmTMB"), "", "R package glmmTMB is not available")),
    data.frame(trait_type="poisson", method="glm_poisson", method_available=TRUE, error_message=""),
    data.frame(trait_type="poisson", method="glmer", method_available=pkg_available("lme4"), error_message=ifelse(pkg_available("lme4"), "", "R package lme4 is not available")),
    data.frame(trait_type="poisson", method="glmmTMB", method_available=pkg_available("glmmTMB"), error_message=ifelse(pkg_available("glmmTMB"), "", "R package glmmTMB is not available")),
    data.frame(trait_type="ordinal", method="MASS::polr", method_available=pkg_available("MASS"), error_message=ifelse(pkg_available("MASS"), "", "R package MASS is not available")),
    data.frame(trait_type="ordinal", method="ordinal::clm", method_available=pkg_available("ordinal"), error_message=ifelse(pkg_available("ordinal"), "", "R package ordinal is not available")),
    data.frame(trait_type="ordinal", method="POLMM", method_available=pkg_available("POLMM"), error_message=ifelse(pkg_available("POLMM"), "", "R package POLMM is not available"))
  )
  tab <- do.call(rbind, rows)
  dir.create(dirname(out), recursive=TRUE, showWarnings=FALSE)
  write.table(tab, out, sep="\t", quote=FALSE, row.names=FALSE)
  cat("[PSR-GWAS] External benchmark method availability written to", out, "\n")
} else if (mode == "ordinal_gwas") {
  method <- if (length(args) >= 2) args[[2]] else stop("method required")
  plink_prefix <- if (length(args) >= 3) args[[3]] else stop("plink_prefix required")
  pheno_file <- if (length(args) >= 4) args[[4]] else stop("pheno_file required")
  out_prefix <- if (length(args) >= 5) args[[5]] else stop("out_prefix required")

  pkg_name <- sub("::.*", "", method)
  if (!pkg_available(pkg_name)) {
    cat(sprintf("R package %s is not available\n", pkg_name))
    quit(save="no", status=0)
  }

  if (grepl("^MASS", method)) suppressPackageStartupMessages(library(MASS))
  else if (grepl("^ordinal", method)) suppressPackageStartupMessages(library(ordinal))

  # Read phenotype
  pheno <- read.table(pheno_file, header=FALSE, stringsAsFactors=FALSE)
  n <- nrow(pheno)
  y_raw <- pheno[,3]
  if (length(unique(y_raw)) < 3) {
    cat(sprintf("[PSR-GWAS] R baseline %s: response has %d levels (<3), skipping\n", method, length(unique(y_raw))))
    summary_line <- data.frame(
      method = method,
      method_available = TRUE,
      runtime_sec = 0.0,
      n_snps = 0,
      valid_snp_count = 0,
      top_snp = NA_character_,
      max_log10p = NA_real_,
      error_message = sprintf("response has %d levels (<3 required)", length(unique(y_raw)))
    )
    out_summary <- paste0(out_prefix, ".r_baseline_summary.tsv")
    write.table(summary_line, out_summary, sep="\t", quote=FALSE, row.names=FALSE)
    cat("[PSR-GWAS] R baseline ordinal GWAS: skipped (insufficient levels)\n")
    quit(save="no", status=0)
  }
  y <- factor(y_raw)

  # Read BIM for SNP metadata
  bim_file <- paste0(plink_prefix, ".bim")
  if (!file.exists(bim_file)) {
    cat(sprintf("BIM file not found: %s\n", bim_file))
    quit(save="no", status=1)
  }
  bim <- read.table(bim_file, header=FALSE, stringsAsFactors=FALSE)
  m <- nrow(bim)
  snp_ids <- bim[,2]

  # Read BED (PLINK binary)
  bed_file <- paste0(plink_prefix, ".bed")
  if (!file.exists(bed_file)) {
    cat(sprintf("BED file not found: %s\n", bed_file))
    quit(save="no", status=1)
  }
  con <- file(bed_file, "rb")
  magic1 <- readBin(con, "raw", 1)
  magic2 <- readBin(con, "raw", 1)
  mode3 <- readBin(con, "raw", 1)
  bytes_per_snp <- ceiling(n / 4)
  G <- matrix(NA_integer_, nrow=n, ncol=m)
  lookup <- c("0"=2L, "1"=NA_integer_, "2"=1L, "3"=0L)
  for (j in seq_len(m)) {
    raw_bytes <- readBin(con, "raw", bytes_per_snp)
    geno <- integer(n)
    for (i in seq_len(n)) {
      b <- as.integer(raw_bytes[(i-1)%/%4 + 1])
      shift <- ((i-1) %% 4) * 2
      code <- (b %/% (2^shift)) %% 4
      geno[i] <- lookup[as.character(code)]
    }
    G[,j] <- geno
  }
  close(con)

  results <- data.frame(
    snp = snp_ids,
    chr = bim[,1],
    bp = bim[,4],
    a1 = bim[,5],
    a2 = bim[,6],
    n_complete = 0L,
    converged = 0L,
    p_value = NA_real_,
    log10p = NA_real_,
    stringsAsFactors=FALSE
  )

  runtime_start <- proc.time()[3]
  complete <- complete.cases(y)
  for (j in seq_len(m)) {
    x <- G[,j]
    dat <- data.frame(y=y, x=x)
    dat <- dat[complete.cases(dat$x, dat$y), , drop=FALSE]
    results$n_complete[j] <- nrow(dat)
    if (nrow(dat) < 10 || length(unique(dat$y)) < 2) next
    pheno_counts <- table(dat$y)
    if (any(pheno_counts < 2)) next
    if (var(as.numeric(dat$x), na.rm=TRUE) == 0) next
    tryCatch({
      if (method == "MASS::polr") {
        fit <- polr(y ~ x, data=dat, method="logistic", Hess=TRUE)
      } else if (method == "ordinal::clm") {
        fit <- clm(y ~ x, data=dat, link="logit")
      }
      cf <- summary(fit)$coefficients
      x_row <- grep("^x$", rownames(cf))
      if (length(x_row) == 0) next
      cn <- colnames(cf)
      se_col <- grep("^[Ss]td", cn, value=TRUE)[1]
      est_col <- grep("^[Ee]st|^[Vv]alue", cn, value=TRUE)[1]
      if (length(se_col) == 0 || length(est_col) == 0) next
      se <- cf[x_row, se_col]
      est <- cf[x_row, est_col]
      if (is.na(se) || se == 0) next
      z <- est / se
      p <- 2 * pnorm(-abs(z))
      results$converged[j] <- 1L
      results$p_value[j] <- p
      results$log10p[j] <- -log10(p)
    }, error=function(e) NULL)
  }
  runtime_sec <- proc.time()[3] - runtime_start

  out_gwas <- paste0(out_prefix, ".r_gwas.tsv")
  write.table(results, out_gwas, sep="\t", quote=FALSE, row.names=FALSE)

  max_log10p <- suppressWarnings(max(results$log10p, na.rm=TRUE))
  if (is.infinite(max_log10p) || is.na(max_log10p)) max_log10p <- NA_real_
  valid_count <- sum(results$converged)
  top_snp <- if (valid_count > 0) results$snp[which.max(results$log10p)] else NA_character_
  summary_line <- data.frame(
    method = method,
    method_available = TRUE,
    runtime_sec = runtime_sec,
    n_snps = m,
    valid_snp_count = valid_count,
    top_snp = top_snp,
    max_log10p = max_log10p,
    error_message = ""
  )
  out_summary <- paste0(out_prefix, ".r_baseline_summary.tsv")
  write.table(summary_line, out_summary, sep="\t", quote=FALSE, row.names=FALSE)
  cat(sprintf("[PSR-GWAS] R baseline ordinal GWAS written for %s: %d/%d valid SNPs\n", method, valid_count, m))
} else if (mode == "ordinal_availability") {
  out <- if (length(args) >= 2) args[[2]] else "tests/results/psr_v05_ordinal_benchmark/v05_ordinal_r_baseline_availability.tsv"
  rows <- list(
    data.frame(
      method = "MASS::polr",
      method_available = pkg_available("MASS"),
      error_message = ifelse(pkg_available("MASS"), "", "R package MASS is not available"),
      supported_trait_types = "ordinal",
      model_type = "cumulative logit (fixed effects)"
    ),
    data.frame(
      method = "ordinal::clm",
      method_available = pkg_available("ordinal"),
      error_message = ifelse(pkg_available("ordinal"), "", "R package ordinal is not available"),
      supported_trait_types = "ordinal",
      model_type = "cumulative logit (fixed effects)"
    ),
    data.frame(
      method = "POLMM",
      method_available = pkg_available("POLMM"),
      error_message = ifelse(pkg_available("POLMM"), "", "R package POLMM is not available"),
      supported_trait_types = "ordinal",
      model_type = "POLMM ordinal mixed model"
    )
  )
  tab <- do.call(rbind, rows)
  dir.create(dirname(out), recursive=TRUE, showWarnings=FALSE)
  write.table(tab, out, sep="\t", quote=FALSE, row.names=FALSE)
  cat("[PSR-GWAS] Ordinal R baseline availability written to", out, "\n")
}
