# fig1cd.R
#
# Faceted boxplot + jitter of variance-explained (R²) for linear vs cubic OLS
# models across four predictor/response data transforms. Legend is inset via
# cowplot.
#
# Usage (import into notebook):
#   source(here("R/fig1cd.R"))
#   plt <- make_fig1cd(r2_summary_df)
#
# `r2_summary_df` is a tidy data frame with columns:
#   model         <chr>   "degree1" | "degree3"
#   data_variant  <chr>   one of the four transform names (see notebook)
#   regulator     <chr>   pTF symbol
#   r2            <dbl>   R²
#   fstat_pval    <dbl>   F-statistic p-value
#
# It is built in the notebook (lines ~580-605) by running fit_ols_model()
# over all combinations of regulator × data_variant × degree.

library(here)
library(tidyverse)
library(cowplot)
source(here("R/theme_ptf.R"))

# ---------------------------------------------------------------------------
# make_fig1cd()
#
# Parameters
# ----------
# r2_summary_df   Tidy R² data frame described in the header.
# ylim            Numeric length-2 vector for y-axis clip (default 0–25%).
# base_size       Font size passed to theme_ptf().
# legend_x, legend_y  Fractional position for the inset legend (cowplot).
#
# Returns a ggdraw object.
# ---------------------------------------------------------------------------
make_fig1cd <- function(
        r2_summary_df,
        ylim      = c(0, 25),
        base_size = 16,
        legend_x  = 0.15,
        legend_y  = 0.70
) {
    plot_df <- r2_summary_df %>%
        mutate(
            data_variant = factor(
                data_variant,
                levels = c("abs_effect_cc_effect", "abs_effect_cc_pval",
                           "rank_rank", "logrank_logrank")
            ),
            model = factor(
                model,
                levels = c("degree1", "degree3"),
                labels = c("Linear", "Cubic")
            ),
            `Variance explained` = r2 * 100
        )

    p <- ggplot(
        plot_df,
        aes(x = data_variant, y = `Variance explained`, fill = data_variant)
    ) +
        geom_boxplot(
            position     = position_dodge(width = 0.75),
            width        = 0.9,
            outlier.shape = NA
        ) +
        geom_jitter(
            position  = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
            alpha     = 0.4,
            size      = 1.2,
            color     = "black"
        ) +
        scale_fill_discrete(
            labels = list(
                expression("|PR|" %~% Enrichment),
                expression("|PR|" %~% -~log(p[binding])),
                expression(Rank("|PR|") %~% Rank(-~log(p[binding]))),
                expression(LRR %~% LRB)
            )
        ) +
        theme_ptf(base_size = base_size) +
        facet_wrap(~model, scales = "free_x", nrow = 1) +
        theme(
            legend.text          = element_text(size = 12),
            legend.title         = element_blank(),
            panel.grid.major.x   = element_blank(),
            axis.ticks.x         = element_blank(),
            axis.text.x          = element_blank(),
            axis.title.x         = element_blank(),
            panel.grid.minor.y   = element_blank(),
            strip.background     = element_rect(fill = "black"),
            strip.text           = element_text(color = "white")
        ) +
        coord_cartesian(ylim = ylim)

    legend <- cowplot::get_legend(
        p + theme(
            legend.position   = "right",
            legend.background = element_rect(fill = alpha("white", 1.0))
        )
    )

    ggdraw() +
        draw_plot(p + theme(legend.position = "none")) +
        draw_grob(legend, x = legend_x, y = legend_y, width = 0.25, height = 0.25)
}

# ---------------------------------------------------------------------------
# save_fig1cd()
# ---------------------------------------------------------------------------
save_fig1cd <- function(path, r2_summary_df, width = 9, height = 7,
                         bg = "white", ...) {
    plt <- make_fig1cd(r2_summary_df, ...)
    ggsave(path, plt, bg = bg, width = width, height = height)
    invisible(plt)
}

# ---------------------------------------------------------------------------
# Standalone entry point
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
    # in general, just run this once before running any fig<> file
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

    model_degrees <- list(degree1 = FALSE, degree3 = TRUE)

    uni_modeling_results <- map(model_degrees, function(cubic_flag) {
        map(model_configs, function(cfg) {
            map(
                .x = set_names(model_regulators),
                .f = fit_ols_model,
                responses    = cfg$response,
                predictors   = cfg$predictor,
                cubic_model  = cubic_flag
            )
        })
    })

    r2_summary_df <- imap_dfr(uni_modeling_results, function(degree_list, model_name) {
        imap_dfr(degree_list, function(regulator_list, data_variant_name) {
            imap_dfr(regulator_list, function(result, regulator_name) {
                fstat <- summary(result$fit)$fstatistic
                tibble(
                    model        = model_name,
                    data_variant = data_variant_name,
                    regulator    = regulator_name,
                    r2           = summary(result$fit)$r.squared,
                    fstat_pval   = pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
                )
            })
        })
    })

    out_path <- here("plots/fig1cd.svg")
    save_fig1cd(out_path, r2_summary_df)
    message("Saved: ", out_path)
}
