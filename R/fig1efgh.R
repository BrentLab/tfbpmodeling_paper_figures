# fig1efgh.R
#
# Four-panel scatter plot for a single pTF showing the raw data with both a
# linear (blue) and cubic (red) fit, across all four predictor/response data
# transforms. R² annotations are placed in the top-right corner of each panel.
#
# Usage (import into notebook):
#   source(here("R/paper_figures/fig1efgh.R"))
#   plt <- make_fig1efgh(uni_modeling_results, regulator = "STB5")
#
# `uni_modeling_results` is the nested list produced by running fit_ols_model()
# in the notebook (lines ~580-605). Structure:
#   uni_modeling_results$degree1$<data_variant>$<regulator>  <- fit_ols_model result
#   uni_modeling_results$degree3$<data_variant>$<regulator>
# Each leaf contains: $fit (lm), $plot (ggplot), $leverage, etc.

library(here)
library(tidyverse)
library(patchwork)
source(here("R/theme_ptf.R"))

# ---------------------------------------------------------------------------
# Internal: build one scatter panel for a single (regulator, data_transform)
# ---------------------------------------------------------------------------
.make_efgh_panel <- function(
        uni_modeling_results,
        regulator,
        data_transform,
        plt_title  = "",
        plt_xlab   = data_transform,
        plt_ylab   = "Response",
        base_size  = 18
) {
    linear_obj <- uni_modeling_results$degree1[[data_transform]][[regulator]]
    cubic_obj  <- uni_modeling_results$degree3[[data_transform]][[regulator]]

    r2_linear <- summary(linear_obj$fit)$r.squared
    r2_cubic  <- summary(cubic_obj$fit)$r.squared

    # The stored plot has three layers: points, smooth, title label.
    # Remove layer 3 (the linear smooth label) and replace with cubic.
    plt <- linear_obj$plot
    plt$layers <- plt$layers[-3]

    plt +
        geom_smooth(
            method  = "lm",
            formula = y ~ poly(x, 3),
            se      = FALSE,
            color   = "#E41A1C"
        ) +
        annotate("text", x = Inf, y = Inf,
                 label    = sprintf("Linear R\u00b2: %.3f", r2_linear),
                 hjust    = 1.5, vjust = 2,
                 size     = 4.2, fontface = "italic", color = "#3366FF") +
        annotate("text", x = Inf, y = Inf,
                 label    = sprintf("Cubic R\u00b2: %.3f", r2_cubic),
                 hjust    = 1.5, vjust = 3.5,
                 size     = 4.2, fontface = "italic", color = "#E41A1C") +
        ggtitle(plt_title) +
        xlab(plt_xlab) +
        ylab(plt_ylab)
}

# ---------------------------------------------------------------------------
# make_fig1efgh()
#
# Parameters
# ----------
# uni_modeling_results   Nested list described in the header.
# regulator              pTF symbol to plot (default "STB5").
# panel_tags             Panel tag letters (default E–H to match paper).
# base_size              Font size passed to theme_ptf().
#
# Returns a patchwork object (4 panels, nrow = 1).
# ---------------------------------------------------------------------------
make_fig1efgh <- function(
        uni_modeling_results,
        regulator  = "STB5",
        panel_tags = c("E", "F", "G", "H"),
        base_size  = 18
) {
    panel_specs <- list(
        list(
            data_transform = "abs_effect_cc_effect",
            plt_ylab       = "|PR|",
            plt_xlab       = "Enrichment"
        ),
        list(
            data_transform = "abs_effect_cc_pval",
            plt_ylab       = "|PR|",
            plt_xlab       = expression(-~log(p[binding]))
        ),
        list(
            data_transform = "rank_rank",
            plt_ylab       = "Rank(|PR|)",
            plt_xlab       = expression(Rank(-~log(p[binding])))
        ),
        list(
            data_transform = "logrank_logrank",
            plt_ylab       = "LRR",
            plt_xlab       = "LRB"
        )
    )

    panels <- map(panel_specs, function(spec) {
        .make_efgh_panel(
            uni_modeling_results,
            regulator      = regulator,
            data_transform = spec$data_transform,
            plt_xlab       = spec$plt_xlab,
            plt_ylab       = spec$plt_ylab,
            base_size      = base_size
        ) + coord_cartesian(clip = "off")
    })

    wrap_plots(panels, nrow = 1) +
        plot_annotation(tag_levels = list(panel_tags))
}

# ---------------------------------------------------------------------------
# save_fig1efgh()
# ---------------------------------------------------------------------------
save_fig1efgh <- function(path, uni_modeling_results, width = 20, height = 4,
                           bg = "white", ...) {
    plt <- make_fig1efgh(uni_modeling_results, ...)
    ggsave(path, plt, bg = bg, width = width, height = height)
    invisible(plt)
}

# ---------------------------------------------------------------------------
# Standalone entry point
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
    # in general, run this once before running any of the fig scripts
    # source(here("R/create_predictors_response_lists.R"))
    source(here("R/fit_ols_model.R"))

    input_data <- pull_predictor_response_lists(pull_data = FALSE, date = "20250805")

    model_regulators <- intersect(
        colnames(input_data$predictors_list$log_rank_normalized),
        colnames(input_data$responses_list$log_rank_normalized)
    )
    model_regulators <- model_regulators[model_regulators != "target_symbol"]

    model_configs <- list(
        abs_effect_cc_effect = list(
            response  = input_data$responses_list$raw,
            predictor = input_data$predictors_list$raw_enrichment
        ),
        abs_effect_cc_pval = list(
            response  = input_data$responses_list$raw,
            predictor = input_data$predictors_list$raw_pvalue
        ),
        rank_rank = list(
            response  = input_data$responses_list$rank,
            predictor = input_data$predictors_list$rank
        ),
        logrank_logrank = list(
            response  = input_data$responses_list$log_rank_normalized,
            predictor = input_data$predictors_list$log_rank_normalized
        )
    )

    uni_modeling_results <- map(list(degree1 = FALSE, degree3 = TRUE), function(cubic_flag) {
        map(model_configs, function(cfg) {
            map(
                .x = set_names(model_regulators),
                .f = fit_ols_model,
                responses   = cfg$response,
                predictors  = cfg$predictor,
                cubic_model = cubic_flag
            )
        })
    })

    out_path <- here("plots/fig1efgh.svg")
    save_fig1efgh(out_path, uni_modeling_results)
    message("Saved: ", out_path)
}
