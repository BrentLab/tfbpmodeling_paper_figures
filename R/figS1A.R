# figS1A.R
#
# Venn diagram of TFs significant (DTO p <= 0.01) under TCC vs ChIP-chip
# binding data. Rendered via reticulate / matplotlib_venn; saves directly to
# file (no ggplot object returned).
#
# Usage (import into a notebook):
#   source(here("R/figS1A.R"))
#   ints <- .compute_s1_intersects(predictors_meta, response_meta, harb_meta, rr_meta)
#   save_figS1A("path/figS1A.svg", ints$harbison, ints$cc)
#
# `predictors_meta` — brent_nf_cc metadata (must have composite_binding col)
# `response_meta`   — mcisaac_oe metadata
# `harb_meta`       — Harbison ChIP-chip PSS metadata (data/{date}/harb_meta.csv)
# `rr_meta`         — rank-response metadata (data/{date}/rr_meta.csv)

library(here)
library(tidyverse)
library(reticulate)

# ---------------------------------------------------------------------------
# .compute_s1_intersects()
#
# Builds harbison_rr_intersect and cc_rr_intersect from raw metadata and
# rr_meta.
#
# Returns a named list: list(harbison = <tibble>, cc = <tibble>)
# Each tibble has columns: regulator_symbol, rank_25, dto_empirical_pvalue,
# binding_source.
# ---------------------------------------------------------------------------
.compute_s1_intersects <- function(predictors_meta, response_meta, harb_meta, rr_meta) {
    valid_harb_ids <- harb_meta %>%
        filter(regulator_symbol %in% unique(predictors_meta$regulator_symbol)) %>%
        filter(regulator_symbol %in% unique(response_meta$regulator_symbol)) %>%
        pull(id)

    harbison_rr_intersect <- rr_meta %>%
        filter(expression %in% response_meta$id) %>%
        filter(promotersetsig %in% valid_harb_ids) %>%
        select(regulator_symbol, rank_25, dto_empirical_pvalue) %>%
        mutate(binding_source = "ChIP-chip")

    cc_is_composite <- predictors_meta %>%
        mutate(is_composite = !is.na(composite_binding)) %>%
        select(id, is_composite)

    cc_rr_intersect <- rr_meta %>%
        filter(expression %in% response_meta$id) %>%
        filter(binding_source == "brent_nf_cc",
               regulator_symbol %in% harbison_rr_intersect$regulator_symbol) %>%
        left_join(cc_is_composite, by = c("promotersetsig" = "id")) %>%
        group_by(regulator_symbol) %>%
        reframe(
            rank_25 = if (any(is_composite, na.rm = TRUE)) {
                rank_25[which(is_composite)[1]]
            } else {
                max(rank_25, na.rm = TRUE)
            },
            dto_empirical_pvalue = if (any(is_composite, na.rm = TRUE)) {
                dto_empirical_pvalue[which(is_composite)[1]]
            } else {
                min(dto_empirical_pvalue, na.rm = TRUE)
            }
        ) %>%
        mutate(binding_source = "TCC")

    list(harbison = harbison_rr_intersect, cc = cc_rr_intersect)
}

# ---------------------------------------------------------------------------
# make_figS1A()
#
# Renders the Venn diagram via matplotlib_venn and saves to `path`. Returns
# the path invisibly.
#
# Parameters
# ----------
# harbison_rr_intersect  tibble — $harbison from .compute_s1_intersects()
# cc_rr_intersect        tibble — $cc from .compute_s1_intersects()
# path                   output SVG path
# pval_threshold         significance cutoff (default 0.01)
# ---------------------------------------------------------------------------
make_figS1A <- function(harbison_rr_intersect, cc_rr_intersect, path,
                         pval_threshold = 0.01) {
    py$tcc_syms <- cc_rr_intersect %>%
        filter(dto_empirical_pvalue <= pval_threshold) %>%
        pull(regulator_symbol)
    py$chipchip_syms <- harbison_rr_intersect %>%
        filter(dto_empirical_pvalue <= pval_threshold) %>%
        pull(regulator_symbol)
    py$out_path <- path

    py_run_string("
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

tcc      = set(tcc_syms)
chipchip = set(chipchip_syms)

fig, ax = plt.subplots(figsize=(6.2, 3.6), dpi=300)
v = venn2([tcc, chipchip], set_labels=('TCC', 'ChIP-chip'), ax=ax)

colors = {'10': '#3B4CC0', '01': '#FDE725', '11': '#1FA187'}
for k in ('10', '01', '11'):
    p = v.get_patch_by_id(k)
    if p:
        p.set_facecolor(colors[k])
        p.set_alpha(0.85)
        p.set_edgecolor('0.3')
        p.set_linewidth(1.2)

region_vals = {
    '10': len(tcc - chipchip),
    '01': len(chipchip - tcc),
    '11': len(tcc & chipchip),
}
for k, val in region_vals.items():
    lbl = v.get_label_by_id(k)
    if lbl:
        lbl.set_text(str(val))
        lbl.set_fontsize(11)
        lbl.set_color('#222')

for lbl in v.set_labels:
    if lbl:
        lbl.set_fontsize(12)
        lbl.set_color('#333')

ax.set_axis_off()
plt.tight_layout()
plt.savefig(out_path, bbox_inches='tight', facecolor='white')
plt.close(fig)
")

    invisible(path)
}

# ---------------------------------------------------------------------------
# save_figS1A()
# ---------------------------------------------------------------------------
save_figS1A <- function(path, harbison_rr_intersect, cc_rr_intersect,
                         pval_threshold = 0.01) {
    make_figS1A(harbison_rr_intersect, cc_rr_intersect, path, pval_threshold)
    invisible(path)
}

# ---------------------------------------------------------------------------
# Standalone entry point
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
    data_pull_date       <- "20250805"
    tfbpmodeling_version <- "1.1.0"

    use_virtualenv(here(".venv"), required = TRUE)

    input_dir <- here("data", data_pull_date)
    rr_meta         <- read_csv(file.path(input_dir, "rr_meta.csv"),          show_col_types = FALSE)
    predictors_meta <- read_csv(file.path(input_dir, "brent_nf_cc_meta.csv"), show_col_types = FALSE)
    response_meta   <- read_csv(file.path(input_dir, "mcisaac_oe_meta.csv"),  show_col_types = FALSE)
    harb_meta       <- read_csv(file.path(input_dir, "harb_meta.csv"),        show_col_types = FALSE)

    ints <- .compute_s1_intersects(predictors_meta, response_meta, harb_meta, rr_meta)

    out_path <- here("plots", data_pull_date, tfbpmodeling_version, "figS1A.svg")
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    save_figS1A(out_path, ints$harbison, ints$cc)
    message("Saved: ", out_path)
}
