# figS2A.R
#
# Two-panel supplementary figure illustrating the cubic-fit / interaction-term
# modeling approach for a single (pTF, mTF) pair:
#   Left panel  — LRR ~ cubic(LRB)  (log-rank response vs log-rank binding)
#   Right panel — Cubic residuals ~ pTF * mTF (LRB^2 interaction term)
#
# Usage (import into notebook):
#   source(here("R/figS2A.R"))
#   plt <- make_figS2A(input_data, pTF = "PDR3", mTF = "GAL4")
#
# `input_data` is the list returned by pull_predictor_response_lists() and
# must contain:
#   input_data$responses_list$log_rank_normalized
#   input_data$responses_list$residuals_normalized
#   input_data$predictors_list$log_rank_normalized

library(here)
library(tidyverse)
library(patchwork)
source(here("R/theme_ptf.R"))
source(here("R/fit_ols_model.R"))

# ---------------------------------------------------------------------------
# Internal: assemble the two-column data frame for one (pTF, mTF, response).
# ---------------------------------------------------------------------------
.build_figS2A_data <- function(input_data, pTF, mTF,
                                response_df_name = "log_rank_normalized") {
    create_input_df(
        pTF,
        input_data$responses_list[[response_df_name]],
        input_data$predictors_list$log_rank_normalized
    ) %>%
        as_tibble(rownames = "target_symbol") %>%
        filter(target_symbol != pTF) %>%
        mutate(interaction = !!rlang::sym(pTF) * !!rlang::sym(mTF)) %>%
        select(response, !!rlang::sym(pTF), interaction)
}

# ---------------------------------------------------------------------------
# make_figS2A()
#
# Parameters
# ----------
# input_data   List from pull_predictor_response_lists() (see header).
# pTF          Primary TF symbol (default "PDR3").
# mTF          Mediator TF symbol (default "GAL4").
# ylim_left    y-axis limits for the left (cubic) panel.
# ylim_right   y-axis limits for the right (interaction) panel.
# base_size    Font size passed to theme_ptf().
#
# Returns a patchwork object (2 panels).
# ---------------------------------------------------------------------------
make_figS2A <- function(
        input_data,
        pTF        = "PDR3",
        mTF        = "GAL4",
        ylim_left  = c(0, 1),
        ylim_right = c(-0.2, 0.75),
        base_size  = 18
) {
    df_lrr       <- .build_figS2A_data(input_data, pTF, mTF,
                                        response_df_name = "log_rank_normalized")
    df_residuals <- .build_figS2A_data(input_data, pTF, mTF,
                                        response_df_name = "residuals_normalized")

    left_panel <- ggplot(df_lrr, aes(!!rlang::sym(pTF), response)) +
        geom_point(size = 0.7, alpha = 0.7) +
        geom_smooth(
            method  = "lm",
            formula = y ~ poly(x, 3, raw = TRUE),
            se      = FALSE
        ) +
        coord_cartesian(ylim = ylim_left) +
        ylab(paste0(pTF, " LRR")) +
        xlab(paste0(pTF, " LRB")) +
        theme_ptf(base_size = base_size)

    right_panel <- ggplot(df_residuals, aes(x = interaction, y = response)) +
        geom_point(size = 0.7, alpha = 0.7) +
        coord_cartesian(ylim = ylim_right) +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
        ylab("Cubic Residuals") +
        xlab(bquote(.(pTF) * "*" * .(mTF) ~ (LRB^2))) +
        theme_ptf(base_size = base_size)

    left_panel + right_panel
}

# ---------------------------------------------------------------------------
# save_figS2A()
# ---------------------------------------------------------------------------
save_figS2A <- function(path, input_data, width = 9, height = 5,
                         bg = "white", ...) {
    plt <- make_figS2A(input_data, ...)
    ggsave(path, plt, bg = bg, width = width, height = height)
    invisible(plt)
}

# ---------------------------------------------------------------------------
# Standalone entry point
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
    # in general, just run this once before running any fig scripts
    # source(here("R/create_predictors_response_lists.R"))

    input_data <- pull_predictor_response_lists(pull_data = FALSE, date = "20250805")

    out_path <- here("plots/figS2A.svg")
    save_figS2A(out_path, input_data)
    message("Saved: ", out_path)
}
