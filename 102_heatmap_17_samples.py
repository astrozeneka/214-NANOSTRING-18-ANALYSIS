import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list

# ══════════════════════════════════════════════════════════════════════════════
# GENE EXPRESSION PROFILE HEATMAP — NanoString  (bi-clustered + clinical)
# Port of 308-CMS-AND-ONCOPLOT/301_heatmap_nanostring.py, restricted to the
# same 17-sample cohort used in 101_supervised_umap_17_samples.R:
#
#   Si02, Si07, Si16, Si20, Si25, Si26, Si27, Si30, Si32,
#   Si40, Si44, Si50, Si54, Si61, Si65, Si70, Si72
#
# CMS + clinical annotations are read from the sibling 308-CMS-AND-ONCOPLOT
# project (same sibling-project pattern as 004_batch_merging.R's TCGA
# reference). CMS_compilation.csv there already calls Si27 = CMS2 (its
# geometric-mean-normalized NanoString call and its independent RNA-seq call
# agree), so no manual override is needed here — unlike
# 101_supervised_umap_17_samples.R, which relies on this project's own
# differently-normalized (TPM-based) CMS pipeline where Si27 came back NA.
#
# Raw NanoString counts are read from this project's own data/ directory
# (identical file content to the 308 project's copy).
# ══════════════════════════════════════════════════════════════════════════════

SIRIRAJ_DIR = "../308-CMS-AND-ONCOPLOT"

ALLOWED_SAMPLES = [
    "Si02", "Si07", "Si16", "Si20", "Si25", "Si26", "Si27", "Si30", "Si32",
    "Si40", "Si44", "Si50", "Si54", "Si61", "Si65", "Si70", "Si72",
]

# ── 0. CMS + clinical data (from the sibling 308 project) ────────────────────
cms_df   = pd.read_csv(os.path.join(SIRIRAJ_DIR, "cms/CMS_compilation.csv")).set_index("sample")
clin_raw = pd.read_csv(os.path.join(SIRIRAJ_DIR, "data/Clinical_data_compilation.csv"), header=0)
clin_raw = clin_raw.rename(columns={clin_raw.columns[0]: "sample_id"}).set_index("sample_id")

# ── 1. Load raw NanoString counts, restrict to the 17-sample cohort ──────────
raw = pd.read_csv("data/Nanostring rawdata for 18 patients.csv", index_col=0)
raw.columns = [c.replace(".RCC", "") for c in raw.columns]   # Si02.RCC → Si02
raw = raw[ALLOWED_SAMPLES]

# Si02 → Si-2 mapping
def _norm_sid(s):
    return "Si-" + str(int(s[2:]))

# ── 2. Geometric-mean normalisation (mirrors 101_cms_nanostring.R) ────────────
def geomean(col):
    pos = col[col > 0]
    return np.exp(np.mean(np.log(pos))) if len(pos) > 0 else np.nan

scale_factors = raw.apply(geomean, axis=0)          # one value per sample
ref           = np.exp(np.mean(np.log(scale_factors)))
norm_mat      = raw.div(scale_factors, axis=1) * ref

# ── 3. log2(normalised + 1) ───────────────────────────────────────────────────
log2norm = np.log2(norm_mat + 1)

# ── 4. Row-wise z-score ───────────────────────────────────────────────────────
# NanoString is already a curated panel (~784 genes) — keep all probes
def zscore_row(row):
    s = row.std()
    return (row - row.mean()) / s if s > 0 else row - row.mean()

mat_z = log2norm.apply(zscore_row, axis=1)

# ── 5. Clustering: rows (Ward) + columns sorted by CMS group ─────────────────
row_linkage = linkage(mat_z, method="ward", metric="euclidean")
row_order   = leaves_list(row_linkage)

_cms_rank = {"CMS1": 0, "CMS2": 1, "CMS3": 2, "CMS4": 3}

def _col_cms_key(c):
    v = cms_df.loc[c, "CMS"] if c in cms_df.index else "Unknown"
    v = v if not pd.isna(v) else "Unknown"
    return _cms_rank.get(v, 99)

col_order     = sorted(range(mat_z.shape[1]), key=lambda i: _col_cms_key(mat_z.columns[i]))
mat_clustered = mat_z.iloc[row_order, col_order]

# ── 6. Real clinical annotations ──────────────────────────────────────────────
samples = mat_clustered.columns.tolist()

_loc_map = {
    "rectum":          "Rectum",
    "rectosigmoid":    "Rectosigmoid",
    "sigmoid":         "Sigmoid",
    "descending":      "Descending",
    "splenic flexure": "Splenic flexure",
    "transverse":      "Transverse",
    "hepatic flexure": "Hepatic flexure",
    "ascending":       "Ascending",
}
_stage_map = {1: "I", 2: "II", 3: "III", 4: "IV"}
_grade_map  = {1: "G1", 2: "G2", 3: "G3", 4: "G4"}

def _safe_int(val):
    try:
        return int(float(val))
    except (TypeError, ValueError):
        return None

records = {}
for s in samples:
    sid_clin = _norm_sid(s)
    cms_val  = cms_df.loc[s, "CMS"] if s in cms_df.index else None
    cms_val  = None if pd.isna(cms_val) else cms_val

    if sid_clin in clin_raw.index:
        row      = clin_raw.loc[sid_clin]
        loc_raw  = str(row.get("Tumor location", "")).strip().lower()
        location = next((v for k, v in _loc_map.items() if k in loc_raw), "Unknown")
        stage    = _stage_map.get(_safe_int(row.get("Stage")), "Unknown")
        grade    = _grade_map.get(_safe_int(row.get("Grade[Original]")), "Unknown")
        org_raw  = str(row.get("Cancer organoid form", "")).strip()
        organoid = "Unknown" if org_raw in ("", "nan", "0") else org_raw.title()
    else:
        location = stage = grade = organoid = "Unknown"

    records[s] = {
        "CMS":           cms_val or "Unknown",
        "Location":      location,
        "Stage":         stage,
        "Grade":         grade,
        "Organoid form": organoid,
    }

clinical = pd.DataFrame(records).T.loc[samples]

# ── 7. Colour palettes ────────────────────────────────────────────────────────
cms_colors = {
    "CMS1": "#dca137",
    "CMS2": "#2f6fab", "CMS3": "#bf7ca4", "CMS4": "#459976",
}
stage_colors = {
    "I": "#e8f5e9", "II": "#a5d6a7", "III": "#388e3c", "IV": "#1b5e20",
}
location_colors = {
    "Rectum":           "#FFB000",
    "Rectosigmoid":     "#FE6100",
    "Sigmoid":          "#DC267F",
    "Descending":       "#9c20e8",
    "Splenic flexure":  "#785EF0",
    "Transverse":       "#648FFF",
    "Hepatic flexure":  "#26A69A",
    "Ascending":        "#2E7D32",
}
grade_colors = {
    "G1": "#fff9c4", "G2": "#f9a825", "G3": "#e65100", "G4": "#b71c1c",
}
organoid_colors = {
    "Compact":             "#4fc3f7",
    "Non Differentiation": "#ef6c00",
    "Grape-Liked": "#8e24aa",
    "Cystic": "#43a047"
}

ANNO_ROWS = [
    ("CMS",           cms_colors),
    ("Location",      location_colors),
    ("Stage",         stage_colors),
    # ("Grade",         grade_colors),
    ("Organoid form", organoid_colors),
]

if __name__ == "__main__":

    # ── 8. Layout ─────────────────────────────────────────────────────────────
    n_samples = mat_clustered.shape[1]
    n_anno    = len(ANNO_ROWS)

    fig = plt.figure(figsize=(max(10, n_samples * 0.55), 10))
    gs = gridspec.GridSpec(
        3, 2,
        height_ratios=[0.01, 0.106 * n_anno, 7],
        width_ratios=[0.8, n_samples * 0.55],
        hspace=0.02, wspace=0.02,
        left=0.06, right=0.78, top=0.95, bottom=0.08,
    )

    ax_corner  = fig.add_subplot(gs[0, 0])
    ax_col_den = fig.add_subplot(gs[0, 1])
    ax_anno    = fig.add_subplot(gs[1, 1])
    ax_row_den = fig.add_subplot(gs[2, 0])
    ax_heatmap = fig.add_subplot(gs[2, 1])

    for ax in [ax_corner, ax_anno]:
        ax.set_axis_off()
    ax_col_den.set_axis_off()

    # ── 9. Row dendrogram (left) ──────────────────────────────────────────────
    dendrogram(row_linkage, ax=ax_row_den, orientation="left",
               color_threshold=0, above_threshold_color="black",
               no_labels=True, link_color_func=lambda _: "black")
    ax_row_den.set_axis_off()

    # ── 10. Heatmap ───────────────────────────────────────────────────────────
    vmax = np.nanpercentile(np.abs(mat_clustered.values), 95)
    sns.heatmap(
        mat_clustered,
        ax=ax_heatmap,
        cmap="viridis",
        vmin=-vmax, vmax=vmax,
        xticklabels=True, yticklabels=False,
        cbar=False,
    )
    ax_heatmap.set_xlabel("")
    ax_heatmap.set_ylabel("")
    ax_heatmap.tick_params(axis="x", labelsize=8, rotation=45)
    ax_heatmap.set_xticklabels([
        "Si" + re.search(r"\d+", t.get_text()).group().zfill(3)
        for t in ax_heatmap.get_xticklabels()
    ])
    plt.setp(ax_heatmap.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    ax_cbar  = fig.add_axes([0.825, 0.25, 0.018, 0.10])
    mappable = ax_heatmap.collections[0]
    cbar     = fig.colorbar(mappable, cax=ax_cbar)
    cbar.set_label("Z-score", fontsize=8)
    cbar.ax.tick_params(labelsize=8)

    # ── 11. Annotation strips ─────────────────────────────────────────────────
    bar_h   = 1.0 / n_anno
    for strip_idx, (field, palette) in enumerate(ANNO_ROWS):
        y_bottom = 1.0 - (strip_idx + 1) * bar_h
        for col_idx, sample in enumerate(samples):
            val   = clinical.loc[sample, field]
            color = palette.get(val, "lightgray")
            rect  = Rectangle(
                (col_idx / n_samples, y_bottom),
                1.0 / n_samples, bar_h,
                transform=ax_anno.transAxes,
                facecolor=color, edgecolor="white", linewidth=0.3,
                clip_on=False,
            )
            ax_anno.add_patch(rect)
        ax_anno.text(
            -0.01, y_bottom + bar_h / 2, field,
            transform=ax_anno.transAxes,
            ha="right", va="center", fontsize=8, fontweight="bold",
        )

    ax_anno.set_xlim(0, 1)
    ax_anno.set_ylim(0, 1)
    ax_anno.set_axis_off()

    # ── 12. Legend ────────────────────────────────────────────────────────────
    legend_elements = []
    for i, (field, palette) in enumerate(ANNO_ROWS):
        if i > 0:
            legend_elements.append(Line2D([0], [0], color="none", label=""))
        legend_elements.append(Line2D([0], [0], color="none", label=f"{field}:"))
        for val, color in palette.items():
            legend_elements.append(Patch(facecolor=color, label=f"  {val}"))

    legend = fig.legend(
        handles=legend_elements,
        loc="upper left",
        bbox_to_anchor=(0.80, 0.95),
        frameon=False,
        fontsize=8,
        handlelength=1.2,
    )
    for text in legend.get_texts():
        if text.get_text().endswith(":"):
            text.set_fontweight("bold")

    # ── 13. Save ──────────────────────────────────────────────────────────────
    os.makedirs("results/102_heatmap_17_samples", exist_ok=True)
    out_path = "results/102_heatmap_17_samples/102_heatmap_nanostring_17_samples.png"
    plt.savefig(out_path, dpi=300, bbox_inches="tight", transparent=False)
    print(f"Saved: {out_path}")
