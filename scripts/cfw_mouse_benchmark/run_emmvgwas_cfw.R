args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: run_emmvgwas_cfw.R <pheno.tsv> <Z.csv> <out_prefix>")
}

pheno <- read.delim(args[1], check.names = FALSE, stringsAsFactors = FALSE)
Zdf <- read.csv(args[2], check.names = FALSE, stringsAsFactors = FALSE)
snp_ids <- Zdf[[1]]
Zall <- as.matrix(Zdf[, -1, drop = FALSE])
rownames(Zall) <- snp_ids

ids <- as.character(pheno$IID)
traits <- setdiff(colnames(pheno), c("FID", "IID"))
Pall <- as.matrix(pheno[, traits, drop = FALSE])
storage.mode(Pall) <- "double"
complete <- complete.cases(Pall)
ids <- ids[complete]
P <- Pall[complete, , drop = FALSE]

available <- intersect(ids, colnames(Zall))
if (length(available) < 10) {
  stop("Too few complete samples shared by phenotype and Z matrix")
}
P <- P[match(available, ids), , drop = FALSE]
Z <- Zall[, available, drop = FALSE]

K <- t(Z) %*% Z
K <- K / mean(diag(K))

library(EMmvGWAS)
res <- gwas_all(K, P, Z, name = args[3], num_cores = 1)
res$SNP_ID <- snp_ids[as.integer(res$snp)]
write.csv(res, paste0(args[3], "_association_with_ids.csv"), row.names = FALSE)
write.table(data.frame(samples = nrow(P), traits = ncol(P), snps = nrow(Z)),
            paste0(args[3], "_run_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
