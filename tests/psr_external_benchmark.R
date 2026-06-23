#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
out <- if (length(args) >= 1) args[[1]] else "tests/results/psr_gwas_smoke/psr_method_availability.tsv"

pkg_available <- function(pkg) requireNamespace(pkg, quietly=TRUE)
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
