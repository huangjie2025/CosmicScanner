#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: psr_candidate_refine.R <data_prefix> <assoc> <family> <out> <top_n> [method]")
}

data_prefix <- args[[1]]
assoc_path <- args[[2]]
trait_type <- tolower(args[[3]])
out_path <- args[[4]]
top_n <- as.integer(args[[5]])
method <- if (length(args) >= 6) tolower(args[[6]]) else "glm"

supported_methods <- c("glm", "firth", "glmer", "glmmtmb")
if (!method %in% supported_methods) stop("Unsupported refinement method: ", method)
if (method == "firth" && trait_type != "binary") stop("Firth refinement is currently supported for binary traits only")
if (!trait_type %in% c("binary", "poisson", "binomial")) stop("Unsupported refinement trait type: ", trait_type)

read_table <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
empty_row <- function(snp, elapsed, error_message) {
  data.frame(
    trait_type = trait_type,
    method = method,
    snp = snp,
    beta = NA_real_,
    se = NA_real_,
    wald_or_stat = NA_real_,
    p_value = NA_real_,
    runtime = elapsed,
    converged = FALSE,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}
result_row <- function(snp, beta, se, stat, p, elapsed, converged, error_message = "") {
  data.frame(
    trait_type = trait_type,
    method = method,
    snp = snp,
    beta = beta,
    se = se,
    wald_or_stat = stat,
    p_value = p,
    runtime = elapsed,
    converged = converged,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

assoc <- read_table(assoc_path)
if (!"SNP" %in% names(assoc)) stop("Association file must contain SNP column: ", assoc_path)
p_col <- if ("P" %in% names(assoc)) "P" else if ("P_VALUE" %in% names(assoc)) "P_VALUE" else stop("Association file must contain P or P_VALUE")
assoc <- assoc[is.finite(assoc[[p_col]]), , drop = FALSE]
assoc <- assoc[order(assoc[[p_col]]), , drop = FALSE]
assoc <- head(assoc, top_n)

geno <- read_table(paste0(data_prefix, ".geno.tsv"))
key <- paste(geno$FID, geno$IID, sep = ":")

if (trait_type == "poisson") {
  pheno <- read.table(paste0(data_prefix, ".poisson.pheno"), stringsAsFactors = FALSE)
  names(pheno) <- c("FID", "IID", "y")
  offset <- scan(paste0(data_prefix, ".offset"), quiet = TRUE)
  pheno$offset <- offset
  model_family <- poisson()
  response <- pheno$y
  model_offset <- pheno$offset
  weights <- NULL
} else if (trait_type == "binary") {
  pheno <- read.table(paste0(data_prefix, ".binary.pheno"), stringsAsFactors = FALSE)
  names(pheno) <- c("FID", "IID", "y")
  model_family <- binomial()
  response <- pheno$y
  model_offset <- NULL
  weights <- NULL
} else {
  pheno <- read.table(paste0(data_prefix, ".binomial.pheno"), stringsAsFactors = FALSE)
  names(pheno) <- c("FID", "IID", "event", "trial")
  model_family <- binomial()
  response <- pheno$event / pheno$trial
  model_offset <- NULL
  weights <- pheno$trial
}

pheno_key <- paste(pheno$FID, pheno$IID, sep = ":")
idx <- match(pheno_key, key)
if (anyNA(idx)) stop("Genotype matrix missing phenotype samples")
geno <- geno[idx, , drop = FALSE]

missing_pkg_message <- NULL
if (method == "firth" && !requireNamespace("logistf", quietly = TRUE)) missing_pkg_message <- "R package 'logistf' is required for method=firth"
if (method == "glmer" && !requireNamespace("lme4", quietly = TRUE)) missing_pkg_message <- "R package 'lme4' is required for method=glmer"
if (method == "glmmtmb" && !requireNamespace("glmmTMB", quietly = TRUE)) missing_pkg_message <- "R package 'glmmTMB' is required for method=glmmTMB"

rows <- list()
for (snp in assoc$SNP) {
  if (!snp %in% names(geno)) next
  elapsed <- 0
  if (!is.null(missing_pkg_message)) {
    elapsed <- system.time({})[["elapsed"]]
    rows[[length(rows) + 1]] <- empty_row(snp, elapsed, missing_pkg_message)
    next
  }
  g <- as.numeric(geno[[snp]])
  dat <- data.frame(y = response, g = g, fid = pheno$FID)
  fit <- NULL
  error_message <- ""
  elapsed <- system.time({
    fit <- tryCatch({
      if (method == "firth") {
        logistf::logistf(y ~ g, data = dat)
      } else if (method == "glmer") {
        if (trait_type == "poisson") {
          if (is.null(model_offset)) lme4::glmer(y ~ g + (1 | fid), family = poisson(), data = dat) else lme4::glmer(y ~ g + offset(model_offset) + (1 | fid), family = poisson(), data = dat)
        } else {
          lme4::glmer(y ~ g + (1 | fid), family = binomial(), data = dat, weights = weights)
        }
      } else if (method == "glmmtmb") {
        if (trait_type == "poisson") {
          if (is.null(model_offset)) glmmTMB::glmmTMB(y ~ g + (1 | fid), family = poisson(), data = dat) else glmmTMB::glmmTMB(y ~ g + offset(model_offset) + (1 | fid), family = poisson(), data = dat)
        } else {
          glmmTMB::glmmTMB(y ~ g + (1 | fid), family = binomial(), data = dat, weights = weights)
        }
      } else if (is.null(model_offset)) {
        glm(y ~ g, family = model_family, data = dat, weights = weights)
      } else {
        glm(y ~ g + offset(model_offset), family = model_family, data = dat, weights = weights)
      }
    }, error = function(e) {
      error_message <<- conditionMessage(e)
      NULL
    })
  })[["elapsed"]]
  if (is.null(fit)) {
    rows[[length(rows) + 1]] <- empty_row(snp, elapsed, error_message)
    next
  }

  if (method == "firth") {
    beta <- unname(fit$coefficients["g"])
    se <- unname(sqrt(diag(fit$var))["g"])
    p <- unname(fit$prob["g"])
    stat <- (beta / se)^2
    converged <- is.null(fit$conv) || isTRUE(fit$conv)
  } else if (method == "glmer") {
    co <- summary(fit)$coefficients
    if (!"g" %in% rownames(co)) next
    beta <- unname(co["g", "Estimate"])
    se <- unname(co["g", "Std. Error"])
    stat <- (beta / se)^2
    p <- pchisq(stat, df = 1, lower.tail = FALSE)
    converged <- is.null(lme4::isSingular(fit)) || !lme4::isSingular(fit)
  } else if (method == "glmmtmb") {
    co <- summary(fit)$coefficients$cond
    if (!"g" %in% rownames(co)) next
    beta <- unname(co["g", "Estimate"])
    se <- unname(co["g", "Std. Error"])
    stat <- (beta / se)^2
    p <- pchisq(stat, df = 1, lower.tail = FALSE)
    converged <- isTRUE(fit$sdr$pdHess)
  } else {
    co <- summary(fit)$coefficients
    if (!"g" %in% rownames(co)) next
    beta <- unname(co["g", "Estimate"])
    se <- unname(co["g", "Std. Error"])
    stat <- (beta / se)^2
    p <- pchisq(stat, df = 1, lower.tail = FALSE)
    converged <- isTRUE(fit$converged)
  }
  rows[[length(rows) + 1]] <- result_row(snp, beta, se, stat, p, elapsed, converged)
}

out <- if (length(rows)) do.call(rbind, rows) else data.frame(
  trait_type = character(), method = character(), snp = character(), beta = numeric(), se = numeric(),
  wald_or_stat = numeric(), p_value = numeric(), runtime = numeric(), converged = logical(), error_message = character()
)
write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("candidate_refinement_output\t", out_path, "\n", sep = "")
cat("candidate_refinement_rows\t", nrow(out), "\n", sep = "")
