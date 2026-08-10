# 003_cms_prediction.R
# Self-contained CMS prediction pipeline, mirroring 211-ONCOGENES/005_cms_prediction.R.
# Input : results/002_tpm/nanostring-18patients-log2tpm.csv
#         Rows = gene symbols (NanoString panel, ~769 genes), columns = samples.
# Output: results/003_cms/CMS_labels_nanostring.csv
#
# CAVEAT: DeepCC's functional spectra score MSigDBv7 gene sets using genome-wide
# expression. This NanoString panel covers only ~769 genes, so most gene sets will
# have sparse/no coverage and predictions here are expected to be noisier/less
# reliable than the RNA-seq-based calls in 211-ONCOGENES. Proceeding per request.
#
# HOW TO RUN WITH DOCKER IN THE CURRENT SETTINGS (Windows)
# docker run --rm -it -v "%cd%:/app/214-NANOSTRING-18-ANALYSIS" -v "%cd%\..\174-CMS-CLASSIFIER:/app/174-CMS-CLASSIFIER" -w /app/214-NANOSTRING-18-ANALYSIS 214-nanostring-18-analysis

library(org.Hs.eg.db)
library(reticulate)
use_condaenv("r-tf23", required = TRUE)
library(DeepCC)
library(tensorflow)
library(data.table)

out_dir <- "results/003_cms"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load expression matrix ────────────────────────────────────────────────
mat <- read.csv(
  "results/002_tpm/nanostring-18patients-log2tpm.csv",
  row.names = 1,
  check.names = FALSE
)
# mat: genes x samples (symbols as row names)
cat("Loaded matrix:", nrow(mat), "genes x", ncol(mat), "samples\n")

# ── 2. Map gene symbols -> Entrez IDs ────────────────────────────────────────
sym_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = rownames(mat),
  columns = "ENTREZID",
  keytype = "SYMBOL"
)

# Keep one symbol -> one Entrez; drop unmapped
sym_map <- sym_map[!is.na(sym_map$ENTREZID) & !duplicated(sym_map$SYMBOL), ]
mat     <- mat[sym_map$SYMBOL, ]
rownames(mat) <- sym_map$ENTREZID

# Drop duplicate Entrez IDs (safety check)
mat <- mat[!duplicated(rownames(mat)), ]
cat("After mapping:", nrow(mat), "genes retained\n")

# ── 3. Transpose -> samples x genes (DeepCC format) ──────────────────────────
expr <- t(mat) |> as.data.frame()

# ── 4. Functional spectra ────────────────────────────────────────────────────
fs <- getFunctionalSpectra(expr, geneSets = "MSigDBv7")

# ── 5. CMS prediction ────────────────────────────────────────────────────────
cms_model <- load_DeepCC_model("../174-CMS-CLASSIFIER/deepcc_model/CRC_TCGA")
labels    <- get_DeepCC_label(cms_model, fs, cutoff = 0.5)
probs     <- get_DeepCC_label(cms_model, fs, cutoff = 0.5, prob_mode = TRUE, prob_raw = TRUE)

# ── 6. Export ────────────────────────────────────────────────────────────────
out <- data.frame(sample = rownames(fs), CMS = labels, probs, row.names = NULL)
write.csv(out, file = file.path(out_dir, "CMS_labels_nanostring.csv"), row.names = FALSE)
cat("Done. Results saved to", file.path(out_dir, "CMS_labels_nanostring.csv"), "\n")
print(out)
