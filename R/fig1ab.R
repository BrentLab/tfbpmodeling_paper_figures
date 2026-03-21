# fig1ab.R
#
# Faceted bar charts showing the distribution of TF binding-set sizes under two
# definitions: BH-FDR <= 0.05 (calling-cards p-value) and DTO rank-response
# binding set size.
#
# Usage (import into notebook):
#   source(here("R/fig1ab.R"))
#   plt <- make_fig1ab(tf_target_counts)
#
# `tf_target_counts` is a data frame with columns:
#   regulator_symbol  <chr>
#   n                 <int>   number of bound genes
#   source            <fct>   "FDR <= 0.05" | "DTO"
# It is built in the notebook as follows:
#   tf_target_counts <- input_data$predictors_list$raw_pvalue %>%
#     pivot_longer(-target_symbol, ...) %>%
#     group_by(regulator_symbol) %>%
#     mutate(bh_fdr = p.adjust(10**(-mlog10pvalue), method = "BH")) %>%
#     filter(bh_fdr <= 0.05) %>%
#     count(regulator_symbol, name = "n") %>%
#     mutate(source = "bh") %>%
#     bind_rows(
#       mcisaac_kem_predictor_rr %>%
#         filter(expression_source == "mcisaac_oe") %>%
#         select(regulator_symbol, binding_set_size) %>%
#         rename(n = binding_set_size) %>%
#         mutate(source = "dto")
#     ) %>%
#     group_by(regulator_symbol) %>%
#     filter(n() == 2) %>%               # keep only TFs with both sources
#     ungroup() %>%
#     mutate(source = factor(source,
#                            levels = c("bh", "dto"),
#                            labels = c("FDR <= 0.05", "DTO")))

library(here)
library(tidyverse)
source(here("R/theme_ptf.R"))

# Log2-scale bin boundaries and labels used across the notebook.
.log2_bins   <- c(0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, Inf)
.log2_labels <- c(
    "0", "1", "2\u20133", "4\u20137", "8\u201315", "16\u201331",
    "32\u201363", "64\u2013127", "128\u2013255", "256\u2013511",
    "512\u20131023", "1024+"
)

# ---------------------------------------------------------------------------
# make_fig1ab()
#
# Parameters
# ----------
# tf_target_counts   Data frame described in the header (before binning).
# bins               Numeric vector of cut-points (default: log2 scale).
# bin_labels         Character vector of bin labels (same length as gaps in
#                    bins). Defaults to .log2_labels.
# base_size          Font size passed to theme_ptf().
#
# Returns a ggplot object.
# ---------------------------------------------------------------------------
make_fig1ab <- function(
        tf_target_counts,
        bins       = .log2_bins,
        bin_labels = .log2_labels,
        base_size  = 16
) {
    binned <- tf_target_counts %>%
        mutate(
            bin = cut(
                n,
                breaks        = bins,
                labels        = bin_labels,
                include.lowest = TRUE,
                right          = FALSE
            )
        ) %>%
        group_by(source, bin) %>%
        tally()

    ggplot(binned, aes(bin, n)) +
        geom_bar(stat = "identity", fill = "black", color = "white") +
        xlab("Number of bound genes") +
        ylab("Number of TFs") +
        theme_ptf(base_size = base_size) +
        facet_wrap(vars(source), ncol = 1) +
        theme(
            panel.grid.major.x = element_blank(),
            axis.text.x        = element_text(angle = 45, hjust = 1),
            panel.grid.minor.y = element_blank(),
            strip.background   = element_rect(fill = "black"),
            strip.text         = element_text(color = "white"),
            aspect.ratio       = 1
        )
}

# ---------------------------------------------------------------------------
# save_fig1ab()
# ---------------------------------------------------------------------------
save_fig1ab <- function(path, tf_target_counts, width = 9, height = 9,
                         bg = "white", ...) {
    plt <- make_fig1ab(tf_target_counts, ...)
    ggsave(path, plt, bg = bg, width = width, height = height)
    invisible(plt)
}

# ---------------------------------------------------------------------------
# Standalone entry point
# Requires: data/rr_meta_20250805.csv, data/brent_nf_cc_meta_20250805.csv,
#           data/mcisaac_oe_meta_20250805.csv, data/predictors_list_20250805.rds
# (produced by prepare_data.R)
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
    data_pull_date       <- "20250805"
    tfbpmodeling_version <- "1.0.0"

    # run this separately, in general
    # source(here("R/create_predictors_response_lists.R"))

    rr_meta         <- read_csv(here("data/rr_meta_20250805.csv"))
    predictors_meta <- read_csv(here("data/brent_nf_cc_meta_20250805.csv"))
    response_meta   <- read_csv(here("data/mcisaac_oe_meta_20250805.csv"))
    input_data      <- pull_predictor_response_lists(pull_data = FALSE, date = "20250805")

    mcisaac_kem_predictor_rr <- rr_meta %>%
        filter(
            promotersetsig %in% predictors_meta$id,
            expression_source == "kemmeren_tfko" |
                (expression_source == "mcisaac_oe" & expression %in% response_meta$id)
        )

    tf_target_counts <- input_data$predictors_list$raw_pvalue %>%
        pivot_longer(-target_symbol,
                     names_to  = "regulator_symbol",
                     values_to = "mlog10pvalue") %>%
        group_by(regulator_symbol) %>%
        mutate(bh_fdr = p.adjust(10^(-mlog10pvalue), method = "BH")) %>%
        filter(bh_fdr <= 0.05) %>%
        ungroup() %>%
        count(regulator_symbol, name = "n") %>%
        mutate(source = "bh") %>%
        bind_rows(
            mcisaac_kem_predictor_rr %>%
                filter(expression_source == "mcisaac_oe") %>%
                select(regulator_symbol, binding_set_size) %>%
                dplyr::rename(n = binding_set_size) %>%
                mutate(source = "dto")
        ) %>%
        ungroup() %>%
        group_by(regulator_symbol) %>%
        filter(n() == 2) %>%
        ungroup() %>%
        mutate(source = factor(source,
                               levels = c("bh", "dto"),
                               labels = c("FDR <= 0.05", "DTO")))

    out_path <- here("plots", data_pull_date, tfbpmodeling_version, "fig1ab.svg")
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    save_fig1ab(out_path, tf_target_counts)
    message("Saved: ", out_path)
}
