library(data.table)
library(biomaRt)

out_dir <- "results/002_tpm"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load raw NanoString counts ──────────────────────────────────────────
# Rows = probes (genes + QC controls), columns = 23 samples (18 patients).
raw <- fread("data/Nanostring rawdata for 18 patients.csv")
setnames(raw, "Probe Name", "Gene")

# Drop NanoString QC probes (positive/negative controls) — not real genes,
# not part of the expression signal used for TPM/UMAP.
qc_probes <- grepl("^(POS|NEG)_[A-H]$", raw$Gene)
cat("Dropping", sum(qc_probes), "QC control probes:", paste(raw$Gene[qc_probes], collapse = ", "), "\n")
raw <- raw[!qc_probes]

genes <- raw$Gene
counts <- as.matrix(raw[, -1])
rownames(counts) <- genes
cat("Panel:", nrow(counts), "genes x", ncol(counts), "samples\n")

# ── 2. Effective gene length (needed to turn counts into TPM) ─────────────
# NanoString counts are not sequencing reads, but the reference pipeline's
# TPM values are length-normalized, so we length-normalize here too for a
# comparable unit. Effective length = length of the longest known transcript
# per HGNC symbol (Ensembl biomaRt), cached locally so re-runs don't hit the
# network.
length_cache_path <- file.path(out_dir, "gene_lengths_biomart.csv")

if (file.exists(length_cache_path)) {
  cat("Loading cached gene lengths from", length_cache_path, "\n")
  gene_lengths <- fread(length_cache_path)
} else {
  cat("Querying Ensembl biomaRt for transcript lengths (requires internet)...\n")
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  bm <- getBM(
    attributes = c("hgnc_symbol", "transcript_length"),
    filters = "hgnc_symbol",
    values = genes,
    mart = mart
  )
  setDT(bm)
  gene_lengths <- bm[hgnc_symbol != "", .(length_bp = max(transcript_length)), by = .(Gene = hgnc_symbol)]
  fwrite(gene_lengths, length_cache_path)
  cat("Cached gene lengths to", length_cache_path, "\n")
}

missing_genes <- setdiff(genes, gene_lengths$Gene)
if (length(missing_genes) > 0) {
  cat(length(missing_genes), "genes had no biomaRt length and will be dropped:\n")
  cat(paste(missing_genes, collapse = ", "), "\n")
}

len_map <- setNames(gene_lengths$length_bp, gene_lengths$Gene)
keep <- rownames(counts) %in% names(len_map)
counts <- counts[keep, ]
gene_len_kb <- len_map[rownames(counts)] / 1000
cat("Retained", nrow(counts), "genes with a known length for TPM.\n")

# ── 3. Counts -> TPM ────────────────────────────────────────────────────────
rpk <- counts / gene_len_kb                      # reads per kilobase
tpm <- t(t(rpk) / colSums(rpk) * 1e6)             # normalize each sample to sum to 1e6
log2tpm <- log2(tpm + 1)

# ── 4. Export ────────────────────────────────────────────────────────────────
fwrite(data.table(Gene = rownames(tpm), as.data.table(tpm)),
       file.path(out_dir, "nanostring-18patients-tpm.csv"))
fwrite(data.table(Gene = rownames(log2tpm), as.data.table(log2tpm)),
       file.path(out_dir, "nanostring-18patients-log2tpm.csv"))

cat("Done. Wrote", out_dir, "/nanostring-18patients-tpm.csv and -log2tpm.csv\n")
