
# docker run --rm -it -v "%cd%:/app/214-NANOSTRING-18-ANALYSIS" -v "%cd%\..\211-ONCOGENES:/app/211-ONCOGENES" -w /app/214-NANOSTRING-18-ANALYSIS 214-nanostring-18-analysis


# our data: log2tpm; results/002_tpm/nanostring-18patients-log2tpm.csv

# tcga: log2tpm: ../211-ONCOGENES/tcga_data/log2tpm_tcga.csv

library(sva)

nano <- read.csv("results/002_tpm/nanostring-18patients-log2tpm.csv", row.names = 1, check.names = FALSE)
tcga <- read.csv("../211-ONCOGENES/tcga_data/log2tpm_tcga.csv", row.names = 1, check.names = FALSE)

common_genes <- intersect(rownames(nano), rownames(tcga))

nano <- nano[common_genes, ]
tcga <- tcga[common_genes, ]

expr <- as.matrix(cbind(nano, tcga))
batch <- rep(c("nanostring", "tcga"), times = c(ncol(nano), ncol(tcga)))

expr_combat <- ComBat(dat = expr, batch = batch)

write.csv(expr_combat, "results/004_batch_merging/nanostring_tcga_combat_log2tpm.csv")

