# fig2bc.R
#
# Bootstrap coefficient distribution for a single pTF, comparing Stage 1
# (all-data, 98% CI) and Stage 2 (top-10%, 90% CI) lasso models.
#
# Usage (import into notebook):
#   source(here("R/fig2bc.R"))
#   all_data_coefs <- readRDS(here("data/bootstrap_coefs_all_data_20250805.rds"))
#   topn_coefs     <- readRDS(here("data/bootstrap_coefs_topn_20250805.rds"))
#   ci_all         <- readRDS(here("data/ci_df_all_data_20250805.rds"))
#   ci_topn        <- readRDS(here("data/ci_df_topn_20250805.rds"))
#   plt <- make_fig2bc(all_data_coefs, topn_coefs, ci_all, ci_topn, regulator = "SKO1")
#
# Usage (standalone):
#   Rscript R/paper_figures/fig2bc.R
#   # writes R/paper_figures/fig2bc.svg

library(here)
library(tidyverse)
library(cowplot)
source(here("R/theme_ptf.R"))

# ---------------------------------------------------------------------------
# Internal helper: extract coefficient draws from pre-computed R objects for
# one regulator.
#
# Parameters
# ----------
# bootstrap_coefs_df   Tibble of bootstrap draws (n_bootstraps x n_features).
#                      From bootstrap_coefs_all_data_{date}.rds[[regulator]].
# ci_df_for_regulator  Filtered ci_df rows for this regulator (columns:
#                      interactor, ci_lo, ci_hi, sign).
# interactor_levels    Optional character vector to fix factor order and ensure
#                      specific interactors appear even if not significant here.
#
# Returns a tidy data frame with columns:
#   interactors  <fct>  full "pTF:mTF" label
#   coef         <dbl>  bootstrap draw
#   significant  <fct>  "Significant" | "Not Significant"
# ---------------------------------------------------------------------------
.extract_coefs_df <- function(bootstrap_coefs_df, ci_df_for_regulator,
                              interactor_levels = NULL) {
    sig_coefs <- ci_df_for_regulator$interactor

    bootstrap_coefs_df %>%
        pivot_longer(everything(), names_to = "interactors", values_to = "coef") %>%
        mutate(
            interactors = factor(
                interactors,
                levels = interactor_levels %||% unique(interactors)
            )
        ) %>%
        filter(interactors %in% c(sig_coefs, interactor_levels)) %>%
        mutate(
            significant = factor(
                interactors %in% sig_coefs,
                levels = c(TRUE, FALSE),
                labels = c("Significant", "Not Significant")
            )
        )
}

# ---------------------------------------------------------------------------
# make_fig2bc()
#
# Parameters
# ----------
# all_data_coefs_list  Named list: regulator → tibble (n_bootstraps × n_features).
#                      From readRDS(here("data/bootstrap_coefs_all_data_{date}.rds")).
# topn_coefs_list      Same structure for the top-N stage.
# ci_df_all_data       Tidy tibble: regulator, interactor, ci_lo, ci_hi, sign.
#                      From readRDS(here("data/ci_df_all_data_{date}.rds")).
# ci_df_topn           Same for top-N stage (90% CI).
# regulator            pTF symbol, e.g. "SKO1".
# interactor_levels    Optional character vector of "pTF:mTF" names to display.
#                      Defaults to significant terms in all_data at 98% CI.
# xlim                 Numeric length-2 vector for x-axis limits.
# base_size            ggplot2 base font size passed to theme_ptf().
# legend_x, legend_y  Fractional position for the inset legend (cowplot).
#
# Returns
# -------
# A ggdraw object (cowplot). Can be passed directly to ggsave().
# ---------------------------------------------------------------------------
make_fig2bc <- function(
        all_data_coefs_list,
        topn_coefs_list,
        ci_df_all_data,
        ci_df_topn,
        regulator,
        interactor_levels = NULL,
        xlim      = c(-0.15, 0.20),
        base_size = 18,
        legend_x  = 0.62,
        legend_y  = 0.69
) {
    ci_all_reg  <- ci_df_all_data %>% filter(regulator == !!regulator)
    ci_topn_reg <- ci_df_topn     %>% filter(regulator == !!regulator)

    # If caller did not specify interactor_levels, use significant all-data terms.
    if (is.null(interactor_levels)) {
        interactor_levels <- ci_all_reg$interactor
    }

    coefs_all  <- .extract_coefs_df(
        all_data_coefs_list[[regulator]],
        ci_all_reg,
        interactor_levels = interactor_levels
    ) %>% mutate(stage = "all_data")

    coefs_topn <- .extract_coefs_df(
        topn_coefs_list[[regulator]],
        ci_topn_reg,
        interactor_levels = interactor_levels
    ) %>% mutate(stage = "topn")

    plot_df <- bind_rows(coefs_all, coefs_topn) %>%
        mutate(
            stage = factor(
                stage,
                levels = c("all_data", "topn"),
                labels = c("Stage 1",  "Stage 2")
            )
        ) %>%
        separate(interactors, into = c("tmp", "mTF"), sep = ":", remove = FALSE) %>%
        select(-tmp)

    p <- ggplot(
        plot_df,
        aes(coef, mTF, fill = stage, alpha = significant, linetype = significant)
    ) +
        geom_boxplot(
            aes(stage = forcats::fct_rev(stage)),
            outlier.alpha = TRUE,
            outlier.color = "black"
        ) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        scale_y_discrete(drop = FALSE) +
        scale_alpha_manual(
            values = c("Significant" = 1, "Not Significant" = 1),
            guide  = "none"
        ) +
        scale_fill_manual(
            values = c("Stage 1" = "#F8766D", "Stage 2" = "#00BA38")
        ) +
        guides(
            fill     = guide_legend(order = 1),
            linetype = guide_legend(order = 2)
        ) +
        labs(x = "Coefficient") +
        ggtitle(paste("pTF:", regulator)) +
        theme_ptf(base_size = base_size) +
        coord_cartesian(xlim = xlim) +
        theme(
            legend.box    = "horizontal",
            legend.title  = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank()
        )

    legend <- cowplot::get_legend(
        p + theme(
            legend.position   = "right",
            legend.box        = "horizontal",
            legend.background = element_rect(fill = alpha("white", 1.0))
        )
    )

    p_nolegend <- p + theme(legend.position = "none")

    ggdraw() +
        draw_plot(p_nolegend) +
        draw_grob(legend, x = legend_x, y = legend_y, width = 0.25, height = 0.25)
}

# ---------------------------------------------------------------------------
# save_fig2bc()
#
# Convenience wrapper around make_fig2bc() + ggsave().
# ---------------------------------------------------------------------------
save_fig2bc <- function(
        path,
        all_data_coefs_list,
        topn_coefs_list,
        ci_df_all_data,
        ci_df_topn,
        regulator,
        width  = 8,
        height = 5,
        bg     = "white",
        ...
) {
    plt <- make_fig2bc(all_data_coefs_list, topn_coefs_list,
                       ci_df_all_data, ci_df_topn, regulator, ...)
    ggsave(path, plt, bg = bg, width = width, height = height)
    invisible(plt)
}

# ---------------------------------------------------------------------------
# Standalone entry point
# Run via: Rscript R/fig2bc.R
# Requires data/bootstrap_coefs_all_data_20250805.rds and
#          data/bootstrap_coefs_topn_20250805.rds (produced by prepare_data.R).
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
    data_pull_date       <- "20250805"
    tfbpmodeling_version <- "1.0.0"

    all_data_coefs <- readRDS(here("data/bootstrap_coefs_all_data_20250805.rds"))
    topn_coefs     <- readRDS(here("data/bootstrap_coefs_topn_20250805.rds"))
    ci_all         <- readRDS(here("data/ci_df_all_data_20250805.rds"))
    ci_topn        <- readRDS(here("data/ci_df_topn_20250805.rds"))

    out_path <- here("plots", data_pull_date, tfbpmodeling_version, "fig2bc.svg")
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    save_fig2bc(
        path                = out_path,
        all_data_coefs_list = all_data_coefs,
        topn_coefs_list     = topn_coefs,
        ci_df_all_data      = ci_all,
        ci_df_topn          = ci_topn,
        regulator           = "SKO1",
        interactor_levels   = c("SKO1:GIS1", "SKO1:OPI1")
    )
    message("Saved: ", out_path)
}
