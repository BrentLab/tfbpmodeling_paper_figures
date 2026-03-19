# fig2e.R
#
# Boxplot of cross-validated variance explained (%) for each TF model across
# four pipeline stages: Cubic Model (CV baseline), Stage 1 (all-data lasso),
# Stage 2 (top-N lasso), Stage 3 (final interaction lasso).
#
# Usage (import into notebook):
#   source(here("R/paper_figures/fig2e.R"))
#   plt <- make_fig2e(stages_comp_r2_df)
#
# `stages_comp_r2_df` columns: model <chr>, stage <fct>, var_expl <dbl>
# Built from joining CV R² across stages (see notebook lines ~2378-2400).

library(here)
library(tidyverse)
library(glue)
source(here("R/theme_ptf.R"))

# ---------------------------------------------------------------------------
# make_fig2e()
# ---------------------------------------------------------------------------
make_fig2e <- function(stages_comp_r2_df, base_size = 18) {
    ggplot(stages_comp_r2_df, aes(stage, var_expl)) +
        geom_boxplot(width = 0.9) +
        ylab("Variance explained (%)") +
        theme_ptf(base_size = base_size) +
        theme(
            axis.title.x = element_blank(),
            axis.text.x  = element_text(angle = 45, hjust = 1)
        )
}

# ---------------------------------------------------------------------------
# save_fig2e()
# ---------------------------------------------------------------------------
save_fig2e <- function(path, stages_comp_r2_df, width = 10, height = 7,
                        bg = "white", ...) {
    plt <- make_fig2e(stages_comp_r2_df, ...)
    ggsave(path, plt, bg = bg, width = width, height = height)
    invisible(plt)
}

# ---------------------------------------------------------------------------
# Helper functions (mirror notebook definitions; use input_data from the
# calling environment — loaded in the standalone block below).
# ---------------------------------------------------------------------------

stratification_classification <- function(binding_vector,
                                           bins = c(0, 8, 64, 512, Inf)) {
    if (length(bins) < 2) {
        warning("bins < 2; returning all ones.")
        return(rep(1L, length(binding_vector)))
    }
    labels       <- seq_len(length(bins) - 1)
    binding_rank <- rank(-binding_vector, ties.method = "min")
    binding_bin  <- cut(binding_rank, breaks = bins, labels = labels,
                        right = TRUE)
    as.integer(binding_bin)
}

# CV R² for the cubic (+ optional pre-perturbation) baseline model.
# References `input_data` and `preperturbation_expr` from the parent env.
cubic_model_cv_r2 <- function(pTF, pval_thres = 0.1, kfolds = 4, seed = 42) {
    input_df <- create_input_df(
        pTF,
        input_data$responses_list$log_rank_normalized,
        input_data$predictors_list$log_rank_normalized
    )

    prepert_vector <- preperturbation_expr %>%
        dplyr::select(target_symbol, !!rlang::sym(pTF)) %>%
        filter(target_symbol %in% rownames(input_df)) %>%
        dplyr::rename(red_median = !!rlang::sym(pTF)) %>%
        deframe()

    input_df$red_median <- prepert_vector[rownames(input_df)]

    formula_full <- as.formula(
        paste0("response ~ poly(", pTF, ", degree = 3) + red_median")
    )
    lmod_full    <- lm(formula_full, input_df)
    rm_pval      <- summary(lmod_full)$coefficients[, "Pr(>|t|)"][["red_median"]]

    final_formula <- if (rm_pval <= pval_thres) {
        message(glue("  {pTF}: using model with red_median"))
        formula_full
    } else {
        message(glue("  {pTF}: using model without red_median"))
        as.formula(paste0("response ~ poly(", pTF, ", degree = 3)"))
    }

    set.seed(seed)
    fold_classes <- stratification_classification(input_df[[pTF]])
    foldid       <- caret::createFolds(as.factor(fold_classes),
                                       k = kfolds, list = TRUE,
                                       returnTrain = FALSE)

    fold_rss <- numeric(length(foldid))
    i        <- 0L
    for (fold_name in names(foldid)) {
        i         <- i + 1L
        test_idx  <- foldid[[fold_name]]
        train_idx <- setdiff(seq_len(nrow(input_df)), test_idx)
        train_fit <- lm(final_formula, data = input_df[train_idx, , drop = FALSE])
        y_pred    <- predict(train_fit,
                             newdata = input_df[test_idx, , drop = FALSE])
        fold_rss[i] <- sum((input_df$response[test_idx] - y_pred)^2, na.rm = TRUE)
    }

    tss_full <- sum((input_df$response - mean(input_df$response))^2)
    1 - sum(fold_rss) / tss_full
}

# CV R² for a lasso-selected interaction formula.
# References `input_data` and `stratification_classification` from parent env.
r2_lasso_all_data_models <- function(pTF, interactors, cv = FALSE, kfolds = 4) {
    if (length(interactors) == 0) {
        return(tibble(
            model        = pTF,  train_rss    = NA_real_,
            train_r2     = NA_real_,  train_adj_r2 = NA_real_,
            test_rss     = NA_real_,  test_tss     = NA_real_,
            test_r2      = NA_real_,  n_terms      = 0L,
            fold         = NA_character_
        ))
    }

    formula <- as.formula(
        paste0("response ~ ", paste(interactors, collapse = " + "))
    )
    input_df <- create_input_df(
        pTF,
        input_data$responses_list$residuals_normalized,
        input_data$predictors_list$log_rank_normalized
    )
    n <- nrow(input_df)

    foldid <- if (cv) {
        fold_classes <- stratification_classification(input_df[[pTF]])
        caret::createFolds(as.factor(fold_classes), k = kfolds,
                           list = TRUE, returnTrain = FALSE)
    } else {
        list(Fold1 = integer(0))
    }

    message(glue("  Fitting: {paste(deparse(formula), collapse = '')}"))

    map_dfr(names(foldid), function(fold_name) {
        test_idx  <- foldid[[fold_name]]
        train_idx <- if (length(test_idx) == 0L) seq_len(n) else
            setdiff(seq_len(n), test_idx)

        fit        <- lm(formula, data = input_df[train_idx, , drop = FALSE])
        train_sum  <- summary(fit)
        train_rss  <- sum(residuals(fit)^2)

        if (length(test_idx) > 0L) {
            y_true   <- input_df$response[test_idx]
            y_pred   <- predict(fit, newdata = input_df[test_idx, , drop = FALSE])
            test_rss <- sum((y_true - y_pred)^2)
            test_tss <- sum((y_true - mean(y_true))^2)
            test_r2  <- if (test_tss > 0) 1 - test_rss / test_tss else NA_real_
        } else {
            test_rss <- NA_real_; test_tss <- NA_real_; test_r2 <- NA_real_
        }

        tibble(
            model        = pTF,
            train_rss    = train_rss,
            train_r2     = unname(train_sum$r.squared),
            train_adj_r2 = unname(train_sum$adj.r.squared),
            test_rss     = test_rss,
            test_tss     = test_tss,
            test_r2      = test_r2,
            n_terms      = length(interactors),
            fold         = fold_name
        )
    })
}

# ---------------------------------------------------------------------------
# Standalone entry point
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
    # in general, run this once before any fig scripts
    # source(here("R/create_predictors_response_lists.R"))
    source(here("R/fit_ols_model.R"))   # create_input_df

    ci_all  <- readRDS(here("data/ci_df_all_data_20250805.rds"))
    ci_topn <- readRDS(here("data/ci_df_topn_20250805.rds"))
    stage4  <- readRDS(here("data/stage4_results_residuals_20250805.rds"))

    # ---- Build sig coef lists from CI data (no Python needed) -----------
    significant_all_data_98_df <- ci_all %>%
        group_by(regulator_symbol = regulator) %>%
        summarise(interactors = list(interactor), .groups = "drop") %>%
        mutate(n_interactors = map_int(interactors, length))

    significant_topn_90_df <- ci_topn %>%
        group_by(regulator_symbol = regulator) %>%
        summarise(interactors = list(interactor), .groups = "drop") %>%
        mutate(n_interactors = map_int(interactors, length))

    # ---- Flag consistent interactors ------------------------------------
    topn_signs <- ci_topn %>%
        dplyr::select(model = regulator, interactor, topn_sign = sign)

    stage4_residuals <- stage4 %>%
        filter(!is.na(coef_interactor)) %>%
        left_join(topn_signs, by = c("model", "interactor")) %>%
        mutate(consistent = sign(coef_interactor) == topn_sign)

    # ---- Load R data ----------------------------------------------------
    input_data <- pull_predictor_response_lists(pull_data = FALSE,
                                                date = "20250805")

    preperturbation_expr <- read_csv(here("data/red_median_wide_20250805.csv"))

    # ---- Cubic CV R² ----------------------------------------------------
    message("Running cubic CV models...")
    cubic_model_cv_r2_df <- tibble(
        model    = significant_all_data_98_df$regulator_symbol,
        var_expl = unlist(map(model, cubic_model_cv_r2)) * 100,
        stage    = "cubic_model"
    )

    # ---- TSS from original response -------------------------------------
    original_response_tss <- input_data$responses_list$log_rank_normalized %>%
        pivot_longer(-target_symbol, names_to = "model", values_to = "lrr_norm") %>%
        filter(target_symbol != model) %>%
        group_by(model) %>%
        reframe(tss = var(lrr_norm) * (n() - 1))

    # ---- Stage 1 CV R² --------------------------------------------------
    message("Running Stage 1 (all-data) CV models...")
    rss_tss_res_all_data <- map2(
        significant_all_data_98_df$regulator_symbol,
        significant_all_data_98_df$interactors,
        r2_lasso_all_data_models, cv = TRUE) %>%
        bind_rows() %>%
        group_by(model) %>%
        reframe(sum_test_rss = sum(test_rss)) %>%
        left_join(original_response_tss, by = "model") %>%
        mutate(var_expl = (1 - sum_test_rss / tss) * 100, stage = "stage1")

    # ---- Stage 2 CV R² --------------------------------------------------
    message("Running Stage 2 (top-N) CV models...")
    rss_tss_res_topn <- map2(
        significant_topn_90_df$regulator_symbol,
        significant_topn_90_df$interactors,
        r2_lasso_all_data_models, cv = TRUE
    ) %>%
        bind_rows() %>%
        group_by(model) %>%
        reframe(sum_test_rss = sum(test_rss)) %>%
        left_join(original_response_tss, by = "model") %>%
        mutate(var_expl = (1 - sum_test_rss / tss) * 100, stage = "stage2")

    # ---- Stage 3 CV R² --------------------------------------------------
    message("Running Stage 3 (final lasso) CV models...")
    stage3_input_df <- stage4_residuals %>%
        filter(consistent, complete.cases(.)) %>%
        group_by(model) %>%
        reframe(interactors = list(interactor))

    rss_tss_res_stage3 <- map2(
        stage3_input_df$model,
        stage3_input_df$interactors,
        r2_lasso_all_data_models, cv = TRUE) %>%
        bind_rows() %>%
        group_by(model) %>%
        reframe(sum_test_rss = sum(test_rss)) %>%
        left_join(original_response_tss, by = "model") %>%
        mutate(var_expl = (1 - sum_test_rss / tss) * 100, stage = "stage3")

    # ---- Assemble stages_comp_r2_df -------------------------------------
    stages_comp_r2_df <- cubic_model_cv_r2_df %>%
        dplyr::select(model, var_expl) %>%
        rename(cubic_model = var_expl) %>%
        left_join(rss_tss_res_all_data %>% dplyr::select(model, var_expl) %>%
                      rename(stage1 = var_expl), by = "model") %>%
        left_join(rss_tss_res_topn %>% dplyr::select(model, var_expl) %>%
                      rename(stage2 = var_expl), by = "model") %>%
        left_join(rss_tss_res_stage3 %>% dplyr::select(model, var_expl) %>%
                      rename(stage3 = var_expl), by = "model") %>%
        mutate(
            stage2 = ifelse(is.na(stage2), stage1, stage2),
            stage3 = ifelse(is.na(stage3), stage2, stage3)
        ) %>%
        pivot_longer(-model, names_to = "stage", values_to = "var_expl") %>%
        mutate(stage = factor(
            stage,
            levels = c("cubic_model", "stage1", "stage2", "stage3"),
            labels = c("Cubic Model", "Stage 1", "Stage 2", "Stage 3")
        ))

    # ---- Save figure ----------------------------------------------------
    out_path <- here("plots/fig2e.svg")
    save_fig2e(out_path, stages_comp_r2_df)
    message("Saved: ", out_path)
}
