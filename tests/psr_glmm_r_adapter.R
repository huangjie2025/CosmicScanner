#!/usr/bin/env Rscript

# Adapter hook for the official PSR-GLMM/R package once its source/demo scripts
# are available locally. Until then, this script provides a deterministic R GLM
# baseline on the validation data exported by run_psr_gwas_validation.sh.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: psr_glmm_r_adapter.R <data_prefix> <result_dir>")
data_prefix <- args[[1]]
result_dir <- args[[2]]

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (is.na(this_file) || !nzchar(this_file)) this_file <- "psr_glmm_r_adapter.R"
script_dir <- dirname(normalizePath(this_file, mustWork = FALSE))

source_path <- Sys.getenv("PSR_GLMM_R_SOURCE", unset = "")
if (nzchar(source_path)) {
  if (!file.exists(source_path)) stop("PSR_GLMM_R_SOURCE not found: ", source_path)
  source(source_path)
  message("Sourced official PSR-GLMM/R code from: ", source_path)
  message("Set PSR_GLMM_R_FUNCTION and extend this adapter once the official API is confirmed.")
}

refine_script <- file.path(script_dir, "psr_candidate_refine.R")

run_refine <- function(label, family) {
  assoc <- file.path(result_dir, paste0("cosmic_", label, ".assoc.txt"))
  out <- file.path(result_dir, paste0("r_baseline_", label, ".tsv"))
  rscript <- file.path(R.home("bin"), "Rscript")
  cmd <- sprintf("%s %s %s %s %s %s 240",
                 shQuote(rscript), shQuote(refine_script), shQuote(data_prefix), shQuote(assoc), shQuote(family), shQuote(out))
  start <- proc.time()[["elapsed"]]
  status <- system(cmd)
  runtime <- proc.time()[["elapsed"]] - start
  if (status != 0) stop("R baseline failed for ", label)
  data.frame(trait = label, output = out, runtime_seconds = runtime)
}

timing <- rbind(run_refine("binary", "binary"), run_refine("poisson", "poisson"), run_refine("binomial", "binomial"))
write.table(timing, file.path(result_dir, "psr_glmm_r_adapter.timing.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
