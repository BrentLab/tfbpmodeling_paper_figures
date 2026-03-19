# fig2d.R
#
# Distribution of the number of interactors per TF model across pipeline stages
# (Stage 1 = all-data 98% CI, Stage 2 = top-N 90% CI, Stage 3 = final lasso).
# Three visual variants via the `variant` argument:
#   "lineplot"  (default / canonical)
#   "barplot"   grouped bars
#   "boxplot"   boxplot + jitter
#
# Usage (import into notebook):
#   source(here("R/fig2d.R"))
#   plt <- make_fig2d(all_vs_topn_comp_df, stage3_consistent_df)
#
# `all_vs_topn_comp_df` columns: regulator_symbol, topn_n_interactors,
#   all_data_n_interactors, diff
#
# `stage3_consistent_df` is stage4_20250805_residuals filtered to
#   consistent == TRUE and complete.cases(), columns: model, interactor, ...

library(here)
library(tidyverse)
source(here("R/theme_ptf.R"))

# ---------------------------------------------------------------------------
# Internal: shared long-format data prep
# ---------------------------------------------------------------------------
.build_fig2d_long <- function(all_vs_topn_comp_df, stage3_consistent_df) {
    all_vs_topn_comp_df %>%
        select(regulator_symbol, topn_n_interactors, all_data_n_interactors) %>%
        rename(model = regulator_symbol) %>%
        left_join(
            stage3_consistent_df %>%
                group_by(model) %>%
                count(name = "stage4_n_interactors"),
            by = "model"
        ) %>%
        replace_na(list(stage4_n_interactors = 0)) %>%
        pivot_longer(-model, names_to = "stage", values_to = "n_interactors") %>%
        mutate(
            stage = factor(
                stage,
                levels = c("all_data_n_interactors",
                           "topn_n_interactors",
                           "stage4_n_interactors"),
                labels = c("stage1", "stage2", "stage3")
            )
        )
}

# ---------------------------------------------------------------------------
# make_fig2d()
# ---------------------------------------------------------------------------
make_fig2d <- function(
        all_vs_topn_comp_df,
        stage3_consistent_df,
        variant   = c("lineplot", "barplot", "boxplot"),
        base_size = 18
) {
    variant <- match.arg(variant)
    df_long <- .build_fig2d_long(all_vs_topn_comp_df, stage3_consistent_df)

    stage_colors <- c(
        "stage1" = "#F8766D",
        "stage2" = "#00BA38",
        "stage3" = "#619CFF"
    )

    if (variant == "lineplot") {
        df_plot <- df_long %>%
            count(stage, n_interactors) %>%
            complete(stage, n_interactors = 0:11, fill = list(n = 0)) %>%
            group_by(stage) %>%
            mutate(proportion = n / sum(n))

        ggplot(df_plot, aes(x = n_interactors, y = n, color = stage)) +
            geom_line(linewidth = 1.2, alpha = 0.8) +
            geom_point(size = 2.5, alpha = 0.9) +
            scale_color_manual(
                values = stage_colors,
                labels = c("Stage 1", "Stage 2", "Stage 3")
            ) +
            scale_x_continuous(breaks = 0:11, limits = c(0, 11)) +
            scale_y_continuous(breaks = seq(0, 60, 10)) +
            coord_cartesian(ylim = c(0, 60)) +
            theme_ptf(base_size = base_size) +
            theme(
                legend.position        = c(0.95, 0.95),
                legend.justification   = c("right", "top"),
                legend.title           = element_blank(),
                legend.background      = element_rect(fill = "white",
                                                      color = "black",
                                                      linewidth = 0.5),
                legend.margin          = margin(6, 6, 6, 6),
                panel.grid.minor.x     = element_blank(),
                panel.grid.minor.y     = element_blank()
            ) +
            labs(x = "Number of interactors", y = "Number of models")

    } else if (variant == "barplot") {
        df_plot <- df_long %>%
            count(stage, n_interactors) %>%
            complete(stage, n_interactors, fill = list(n = 0))

        ggplot(df_plot, aes(x = n_interactors, y = n, fill = stage)) +
            geom_col(
                position = position_dodge(preserve = "single"),
                width    = 0.7,
                color    = "black"
            ) +
            scale_fill_manual(
                values = stage_colors,
                breaks = c("stage1", "stage2", "stage3")
            ) +
            scale_x_continuous(breaks = seq(0, 11)) +
            theme_ptf(base_size = base_size) +
            theme(
                legend.position.inside = c(0.7, 0.7),
                panel.grid.minor       = element_blank(),
                legend.title           = element_blank()
            ) +
            labs(x = "Number of Interactors", y = "Count of Models")

    } else {
        ggplot(df_long, aes(x = stage, y = n_interactors, fill = stage)) +
            geom_boxplot(alpha = 0.7, outlier.shape = NA) +
            geom_point(
                position = position_jitter(width = 0.2, height = 0),
                alpha = 0.6, size = 1
            ) +
            scale_fill_manual(values = stage_colors) +
            theme_ptf(base_size = base_size) +
            theme(legend.position = "none", panel.grid.minor = element_blank()) +
            labs(x = "Stage", y = "Number of Interactors")
    }
}

# ---------------------------------------------------------------------------
# save_fig2d()
# ---------------------------------------------------------------------------
save_fig2d <- function(path, all_vs_topn_comp_df, stage3_consistent_df,
                        width = 8, height = 7, bg = "white", ...) {
    plt <- make_fig2d(all_vs_topn_comp_df, stage3_consistent_df, ...)
    ggsave(path, plt, bg = bg, width = width, height = height)
    invisible(plt)
}

# ---------------------------------------------------------------------------
# Standalone entry point
# Requires: data/ci_df_all_data_20250805.rds, data/ci_df_topn_20250805.rds,
#           data/stage4_results_residuals_20250805.rds
# (produced by prepare_data.R)
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
    ci_all  <- readRDS(here("data/ci_df_all_data_20250805.rds"))
    ci_topn <- readRDS(here("data/ci_df_topn_20250805.rds"))
    stage4  <- readRDS(here("data/stage4_results_residuals_20250805.rds"))

    # ---- Interactor counts per stage -------------------------------------
    significant_all_data_98_df <- ci_all %>%
        group_by(regulator_symbol = regulator) %>%
        summarise(all_data_n_interactors = n(), .groups = "drop")

    significant_topn_90_df <- ci_topn %>%
        group_by(regulator_symbol = regulator) %>%
        summarise(topn_n_interactors = n(), .groups = "drop")

    # ---- Flag consistent interactors (sign matches between stages) -------
    topn_signs <- ci_topn %>%
        select(model = regulator, interactor, topn_sign = sign)

    stage4_residuals <- stage4 %>%
        filter(!is.na(coef_interactor)) %>%
        left_join(topn_signs, by = c("model", "interactor")) %>%
        mutate(consistent = sign(coef_interactor) == topn_sign)

    # ---- Build function inputs -------------------------------------------
    all_vs_topn_comp_df <- significant_topn_90_df %>%
        full_join(significant_all_data_98_df, by = "regulator_symbol") %>%
        replace_na(list(all_data_n_interactors = 0L, topn_n_interactors = 0L)) %>%
        mutate(diff = topn_n_interactors - all_data_n_interactors)

    stage3_consistent_df <- stage4_residuals %>%
        filter(consistent, complete.cases(.))

    # ---- Save figure ------------------------------------------------------
    out_path <- here("plots/fig2d.svg")
    save_fig2d(out_path, all_vs_topn_comp_df, stage3_consistent_df)
    message("Saved: ", out_path)
}
