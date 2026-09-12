# 101_supervised_umap_17_samples.R
# Same recipe as 005_supervised_umap.R (PCA30 -> uwot supervised UMAP, same
# ellipse-label placement, same visual style), but restricted to the 17
# NanoString samples that are common to both this project's UMAP and the
# 308-CMS-AND-ONCOPLOT heatmap (301_heatmap_nanostring.py):
#
#   Si02, Si07, Si16, Si20, Si25, Si26, Si27, Si30, Si32,
#   Si40, Si44, Si50, Si54, Si61, Si65, Si70, Si72
#
# Si27 is forced to CMS2 here: DeepCC on this project's TPM-normalized
# NanoString expression called it NA (unclassified at cutoff 0.5), but both
# the 308 project's geometric-mean-normalized NanoString call AND its
# independent RNA-seq-based call agree on CMS2, so we take that as the label.
#
# TCGA reference samples are left untouched (still plotted as background).

library(data.table)
library(uwot)
library(ggplot2)
library(ggrepel)

cms_colors <- c(CMS1 = "#E89D33", CMS2 = "#0072AC", CMS3 = "#D079A4", CMS4 = "#009D76")

out_dir <- "results/101_umap_supervised_17_samples"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Args: Rscript 101_supervised_umap_17_samples.R -tw <target_weight> -seed <seed> -el <ellipse_level>
args <- commandArgs(trailingOnly = TRUE)
parse_named_args <- function(args, defaults) {
  opts <- defaults
  i <- 1
  while (i <= length(args)) {
    flag <- sub("^-+", "", args[[i]])
    if (flag %in% names(defaults)) {
      opts[[flag]] <- args[[i + 1]]
      i <- i + 2
    } else {
      stop("Unknown argument: ", args[[i]])
    }
  }
  opts
}
opts <- parse_named_args(args, list(tw = "0.05", seed = "48", el = "0.80"))

target_weight <- as.numeric(opts$tw)
umap_seed <- as.integer(opts$seed)
ellipse_level <- as.numeric(opts$el)

tw_tag <- sprintf("tw%03d", target_weight * 1000)
seed_tag <- sprintf("seed%d", umap_seed)
el_tag <- sprintf("el%03d", ellipse_level * 100)

message(sprintf("[101_supervised_umap_17_samples] target_weight = %s, seed = %s, ellipse_level = %s",
                 target_weight, umap_seed, ellipse_level))

# ── 1. Load data ──────────────────────────────────────────────────────────────
df <- fread("results/004_batch_merging/nanostring_tcga_combat_log2tpm.csv")
mat <- as.matrix(df[, -1, with = FALSE])
rownames(mat) <- df[[1]]
mat <- t(mat)  # samples in rows, genes in columns

nano_cms <- fread("results/003_cms/CMS_labels_nanostring.csv")
tcga_cms <- fread("data/CMS-TCGA.tsv")
cms_dict <- setNames(tcga_cms$deepcms_subtype, tcga_cms$patient_id)
cms_dict[nano_cms$sample] <- nano_cms$CMS  # our labels take priority over TCGA's

# Manual override: DeepCC called Si27 NA under this project's TPM-normalized
# preprocessing, but the geometric-mean-normalized NanoString call and the
# independent RNA-seq call (308-CMS-AND-ONCOPLOT project) both agree on CMS2.
cms_dict["Si27.RCC"] <- "CMS2"

# Batch: TCGA (public reference) vs Internal (our own NanoString samples)
# Matrix columns use the bare patient barcode (e.g. "TCGA-AA-3561"), matching
# patient_id, not sample_id which carries a "-01A" sample-type suffix.
batch <- setNames(rep("TCGA", nrow(tcga_cms)), tcga_cms$patient_id)

# Restrict the Internal (NanoString) cohort to the 17 samples common to both
# this project's UMAP and the 308 project's heatmap. TCGA reference samples
# are unaffected.
allowed_internal <- c(
  "Si02.RCC", "Si07.RCC", "Si16.RCC", "Si20.RCC", "Si25.RCC", "Si26.RCC",
  "Si27.RCC", "Si30.RCC", "Si32.RCC", "Si40.RCC", "Si44.RCC", "Si50.RCC",
  "Si54.RCC", "Si61.RCC", "Si65.RCC", "Si70.RCC", "Si72.RCC"
)

samples <- rownames(mat)
is_tcga <- samples %in% names(batch)
keep_cohort <- is_tcga | samples %in% allowed_internal

cms_labels <- cms_dict[samples]
keep <- keep_cohort & cms_labels %in% names(cms_colors)  # drop NA, "", "NOLBL", and excluded Internal samples
mat <- mat[keep, ]
cms_labels <- factor(cms_labels[keep])
samples <- rownames(mat)
batch_labels <- ifelse(samples %in% names(batch), "TCGA", "Siriraj Cohort")

stopifnot(sum(batch_labels == "Siriraj Cohort") == length(allowed_internal))

# ── 2. PCA30 → soft supervised UMAP (uwot) ───────────────────────────────────
expression_matrix <- scale(mat)
pca30 <- prcomp(expression_matrix, rank. = 30)$x

set.seed(umap_seed)
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

# Match the sample-label style used in 308-CMS-AND-ONCOPLOT/301_heatmap_nanostring.py:
# "Si" + zero-padded 3-digit number (e.g. Si02.RCC -> Si002, Si27.RCC -> Si027).
# Only applies to Siriraj Cohort (NanoString) samples; TCGA barcodes are left as-is.
format_heatmap_style <- function(id) {
  paste0("Si", formatC(as.integer(sub("^Si0*(\\d+)\\..*$", "\\1", id)), width = 3, flag = "0"))
}
umap_df[, SampleLabel := ifelse(Batch == "Siriraj Cohort", format_heatmap_style(SampleId), SampleId)]

fwrite(umap_df, file.path(out_dir, paste0("umap_supervised_cms_", tw_tag, "_", seed_tag, "_", el_tag, ".csv")))

# ── 3. Plot ───────────────────────────────────────────────────────────────────
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
  geom_text_repel(data = umap_df[Batch == "Siriraj Cohort"], aes(label = SampleLabel),
                   size = 3, fontface = "bold", show.legend = FALSE, max.overlaps = Inf,
                   box.padding = 1.2, min.segment.length = 0,
                   segment.color = "grey40", segment.size = 0.4) +
  geom_label_repel(data = cms_label_pos, aes(label = CMS, shape = NULL, alpha = NULL),
                    nudge_x = cms_label_pos$nudge_x, nudge_y = cms_label_pos$nudge_y,
                    fontface = "bold", size = 4, box.padding = 0.5,
                    show.legend = FALSE, max.overlaps = Inf, min.segment.length = 0,
                    segment.color = "grey20", segment.size = 0.6, seed = 42) +
  scale_color_manual(values = cms_colors) +
  scale_shape_manual(values = c(TCGA = 16, `Siriraj Cohort` = 17)) +
  scale_alpha_manual(values = c(TCGA = 0.4, `Siriraj Cohort` = 0.9)) +
  guides(color = "none") +
  labs(title = paste0("Supervised UMAP, 17-sample cohort (uwot, PCA30, target_weight = ", target_weight, ", seed = ", umap_seed, ", ellipse_level = ", ellipse_level, ")")) +
  theme_minimal()

ggsave(file.path(out_dir, paste0("umap_supervised_cms_", tw_tag, "_", seed_tag, "_", el_tag, ".png")), p, width = 8, height = 8, dpi = 150)
