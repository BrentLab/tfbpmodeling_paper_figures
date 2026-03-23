# fig2d.R
#
# Distribution of the number of interactors per TF model across pipeline stages
# (Stage 1 = all-data 98% CI, Stage 2 = top-N 90% CI, Stage 3 = chosen column).
# Three visual variants via the `plot_variant` argument:
#   "lineplot"  (default / canonical)
#   "barplot"   grouped bars
#   "boxplot"   boxplot + jitter
#
# Usage (import into notebook):
#   source(here("R/fig2d.R"))
#   plt <- make_fig2d(stage_comp_df, stage3_col = "stage3_lassocv_n_interactors")
#
# `stage_comp_df` columns: regulator_symbol, all_data_n_interactors,
#   topn_n_interactors, stage3_lassocv_n_interactors,
#   stage3_bootstrap_n_interactors (may be NA if bootstrap CI not available)
#
# `stage3_col` selects which stage3 count column to use as the Stage 3 series.

library(here)
library(tidyverse)
library(glue)
source(here("R/theme_ptf.R"))

#' create the long format dataframe with columns model, stage, n_interactors
#'
#' @param stage_comp_df wide format with columns regulator_symbol,
#'    all_data_n_interactors, topn_n_interactors, stage3_lassocv_n_interactors,
#'    stage3_bootstrap_n_interactors
#' @param stage3_col which of the stage3 columns to use
#' @param include_0_interactors Default TRUE. Whether or not to include
#'     the 0 interactor value on in the dataframe
#'
#' @return a dataframe in long format with columns model
#'     (renamed regulator_symbol), stage, n_interactors
.build_fig2d_long <- function(stage_comp_df,
                              stage3_col,
                              include_0_interactors = TRUE) {
    message(glue::glue("Using column {stage3_col}"))
    out = stage_comp_df %>%
        select(model = regulator_symbol,
               all_data_n_interactors,
               topn_n_interactors,
               stage3_n_interactors = all_of(stage3_col)) %>%
        pivot_longer(-model, names_to = "stage", values_to = "n_interactors") %>%
        mutate(
            stage = factor(
                stage,
                levels = c("all_data_n_interactors",
                           "topn_n_interactors",
                           "stage3_n_interactors"),
                labels = c("stage1", "stage2", "stage3")
            )
        )

    if(include_0_interactors){
        out
    } else{
        filter(out, n_interactors != 0)
    }

}

# ---------------------------------------------------------------------------
# make_fig2d()
# ---------------------------------------------------------------------------
make_fig2d <- function(
        stage_comp_df,
        stage3_col    = c("stage3_lassocv_n_interactors", "stage3_bootstrap_n_interactors"),
        include_0_interactors = TRUE,
        plot_variant  = c("barplot", "lineplot", "boxplot"),
        base_size     = 18
) {
    plot_variant <- match.arg(plot_variant)
    stage3_col = match.arg(stage3_col)
    df_long <- .build_fig2d_long(stage_comp_df, stage3_col, include_0_interactors)

    stage_colors <- c(
        "stage1" = "#F8766D",
        "stage2" = "#00BA38",
        "stage3" = "#619CFF"
    )

    if (plot_variant == "lineplot") {
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

    } else if (plot_variant == "barplot") {
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
save_fig2d <- function(path, stage_comp_df,
                        width = 8, height = 7, bg = "white", ...) {
    plt <- make_fig2d(stage_comp_df, ...)
    ggsave(path, plt, bg = bg, width = width, height = height)
    invisible(plt)
}

# ---------------------------------------------------------------------------
# Standalone entry point — sources prepare_data.R to get stage_comp_df
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
    source(here("R/prepare_data.R"))

    out_path <- here("plots", data_pull_date, tfbpmodeling_version, "fig2d.svg")
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    save_fig2d(out_path, stage_comp_df, plot_variant="barplot", stage3_col = "stage3_bootstrap_n_interactors")
    message("Saved: ", out_path)
}
