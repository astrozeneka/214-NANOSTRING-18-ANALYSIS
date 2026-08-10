# 005_supervised_umap.R
# Supervised UMAP colored by CMS, mirroring ../211-ONCOGENES/109_plot_umap_v3.R
# exactly (same PCA30 -> uwot supervised UMAP recipe, same ellipse-label placement,
# same visual style), applied to the ComBat-merged NanoString + TCGA log2TPM matrix.

library(data.table)
library(uwot)
library(ggplot2)
library(ggrepel)

cms_colors <- c(CMS1 = "#E89D33", CMS2 = "#0072AC", CMS3 = "#D079A4", CMS4 = "#009D76")

out_dir <- "results/005_umap_supervised"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

target_weight <- 0.17
tw_tag <- sprintf("tw%03d", target_weight * 1000)

# ── 1. Load data ──────────────────────────────────────────────────────────────
df <- fread("results/004_batch_merging/nanostring_tcga_combat_log2tpm.csv")
mat <- as.matrix(df[, -1, with = FALSE])
rownames(mat) <- df[[1]]
mat <- t(mat)  # samples in rows, genes in columns

nano_cms <- fread("results/003_cms/CMS_labels_nanostring.csv")
tcga_cms <- fread("data/CMS-TCGA.tsv")
cms_dict <- setNames(tcga_cms$deepcms_subtype, tcga_cms$patient_id)
cms_dict[nano_cms$sample] <- nano_cms$CMS  # our labels take priority over TCGA's

# Batch: TCGA (public reference) vs Internal (our own NanoString samples)
# Matrix columns use the bare patient barcode (e.g. "TCGA-AA-3561"), matching
# patient_id, not sample_id which carries a "-01A" sample-type suffix.
batch <- setNames(rep("TCGA", nrow(tcga_cms)), tcga_cms$patient_id)

samples <- rownames(mat)

cms_labels <- cms_dict[samples]
keep <- cms_labels %in% names(cms_colors)  # drop NA, "", and "NOLBL"
mat <- mat[keep, ]
cms_labels <- factor(cms_labels[keep])
samples <- rownames(mat)
batch_labels <- ifelse(samples %in% names(batch), "TCGA", "Internal")

# ── 2. PCA30 → soft supervised UMAP (uwot) ───────────────────────────────────
expression_matrix <- scale(mat)
pca30 <- prcomp(expression_matrix, rank. = 30)$x

set.seed(48)
coords <- umap(
  X = pca30,             # Samples in rows, features/PCs in columns
  y = cms_labels,         # Factor of CMS labels
  n_components = 3,      # Output 3 UMAP dimensions
  target_weight = target_weight,    # Gentle label guidance (adjust between 0.05 - 0.2)
  n_neighbors = 10,       # Higher values preserve global continuum/proximity
  min_dist = 0.3,         # Prevents points from hyper-clumping
  metric = "cosine",
  n_threads = 1          # multi-threading breaks reproducibility even with a seed
)
colnames(coords) <- c("UMAP1", "UMAP2", "UMAP3")

umap_df <- data.table(SampleId = samples, coords, CMS = as.character(cms_labels), Batch = batch_labels)
fwrite(umap_df, file.path(out_dir, paste0("umap_supervised_cms_", tw_tag, ".csv")))

# ── 3. Plot ───────────────────────────────────────────────────────────────────
ellipse_level <- 0.80

# For each CMS group, place the label on the point of its ellipse border that sits
# farthest from every sample in the plot (i.e. the emptiest spot around the ring).
# This auto-adapts whenever the UMAP layout changes, instead of relying on
# hand-tuned angles that break as soon as the embedding parameters change.
#
# The border itself must match ggplot2's own stat_ellipse(type = "norm") geometry
# exactly, or the anchor point drifts off the drawn ring. Per ggplot2's source
# (StatEllipse$compute_group), the radius for type="norm" is
#   radius <- sqrt(dfn * qf(level, dfn, dfd)),  dfn = 2, dfd = n_group - 1
# an F-distribution correction for finite sample size — NOT the chi-square
# radius sqrt(qchisq(level, 2)), which is only its infinite-dfd limit and differs
# per group since each CMS group has a different n. Using the wrong radius (or a
# single shared one) is exactly what produced the gap between the label's arrow
# and the ellipse.
ellipse_border_point <- function(cx, cy, angle, width, height, label_angle) {
  t <- label_angle * pi / 180
  a <- width / 2
  b <- height / 2
  c(UMAP1 = cx + a * cos(t) * cos(angle) - b * sin(t) * sin(angle),
    UMAP2 = cy + a * cos(t) * sin(angle) + b * sin(t) * cos(angle))
}

best_label_point <- function(cx, cy, angle, width, height, all_x, all_y) {
  candidates <- seq(0, 350, by = 10)
  pts <- t(sapply(candidates, function(a) ellipse_border_point(cx, cy, angle, width, height, a)))
  min_dist <- sapply(seq_len(nrow(pts)), function(i) min(sqrt((all_x - pts[i, 1])^2 + (all_y - pts[i, 2])^2)))
  as.list(pts[which.max(min_dist), ])
}

label_nudge_dist <- 1.2  # extra push (data units) so the label box clears the ellipse entirely

cms_label_pos <- rbindlist(lapply(levels(cms_labels), function(cms) {
  d <- umap_df[CMS == cms]
  cx <- mean(d$UMAP1); cy <- mean(d$UMAP2)
  cov_mat <- cov(cbind(d$UMAP1, d$UMAP2))
  eig <- eigen(cov_mat)
  ord <- order(eig$values, decreasing = TRUE)
  eigenvalues <- eig$values[ord]
  eigenvectors <- eig$vectors[, ord]
  angle <- atan2(eigenvectors[2, 1], eigenvectors[1, 1])
  dfd <- nrow(d) - 1
  radius <- sqrt(2 * qf(ellipse_level, 2, dfd))  # exact stat_ellipse(type="norm") radius for this group's n
  width <- 2 * radius * sqrt(eigenvalues[1])
  height <- 2 * radius * sqrt(eigenvalues[2])
  pt <- best_label_point(cx, cy, angle, width, height, umap_df$UMAP1, umap_df$UMAP2)
  dx <- pt$UMAP1 - cx; dy <- pt$UMAP2 - cy
  d_norm <- sqrt(dx^2 + dy^2)
  data.table(CMS = cms, UMAP1 = pt$UMAP1, UMAP2 = pt$UMAP2,
             nudge_x = dx / d_norm * label_nudge_dist, nudge_y = dy / d_norm * label_nudge_dist)
}))

p <- ggplot(umap_df, aes(UMAP1, UMAP2, color = CMS, shape = Batch, alpha = Batch)) +
  geom_point(size = 2) +
  stat_ellipse(aes(shape = NULL, alpha = NULL), type = "norm", level = ellipse_level,
               linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
  geom_text_repel(data = umap_df[Batch == "Internal"], aes(label = SampleId),
                   size = 3, fontface = "bold", show.legend = FALSE, max.overlaps = Inf,
                   box.padding = 1.2, min.segment.length = 0,
                   segment.color = "grey40", segment.size = 0.4) +
  geom_label_repel(data = cms_label_pos, aes(label = CMS, shape = NULL, alpha = NULL),
                    nudge_x = cms_label_pos$nudge_x, nudge_y = cms_label_pos$nudge_y,
                    fontface = "bold", size = 4, box.padding = 0.5,
                    show.legend = FALSE, max.overlaps = Inf, min.segment.length = 0,
                    segment.color = "grey20", segment.size = 0.6, seed = 42) +
  scale_color_manual(values = cms_colors) +
  scale_shape_manual(values = c(TCGA = 16, Internal = 17)) +
  scale_alpha_manual(values = c(TCGA = 0.4, Internal = 0.9)) +
  guides(color = "none") +
  labs(title = paste0("Supervised UMAP (uwot, PCA30, target_weight = ", target_weight, ")")) +
  theme_minimal()

ggsave(file.path(out_dir, paste0("umap_supervised_cms_", tw_tag, ".png")), p, width = 8, height = 8, dpi = 150)
