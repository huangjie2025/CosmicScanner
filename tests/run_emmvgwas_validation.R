args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: run_emmvgwas_validation.R <P.csv> <Z.csv> <out_prefix>")
}

P <- as.matrix(read.csv(args[1], check.names = FALSE))
Zdf <- read.csv(args[2], check.names = FALSE)
snp_ids <- Zdf[[1]]
Z <- as.matrix(Zdf[, -1, drop = FALSE])
rownames(Z) <- snp_ids

K <- t(Z) %*% Z
K <- K / mean(diag(K))

library(EMmvGWAS)
res <- gwas_all(K, P, Z, name = args[3], num_cores = 1)
res$SNP_ID <- snp_ids[as.integer(res$snp)]
write.csv(res, paste0(args[3], "_association_with_ids.csv"), row.names = FALSE)
