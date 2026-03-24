# prepare_data.R
#
# Loads all serialized data from data/<data_pull_date>/ and
# data/<data_pull_date>/<tfbpmodeling_version>/ into memory, applies
# predictor/response transforms, and builds all derived objects needed as
# inputs to the fig scripts. No htcf_ref mount or Python required.
#
# Prerequisites: the target data directories must already exist and be
# populated. If data/<data_pull_date>/ is missing, run:
#   Rscript R/prepare_input_data.R
# If data/<data_pull_date>/<tfbpmodeling_version>/ is missing, run:
#   Rscript R/prepare_results_data.R   # for tfbpmodeling >= 1.0.0
#   Rscript R/prepare_results_data.R      # for legacy results
#
# `data_pull_date` is the date string identifying which versioned data pull
# from the yeast database to use (e.g. "20250805").
#
# `tfbpmodeling_version` is the version of the tfbpmodeling software that was
# used to produce the results being loaded (e.g. "1.0.0").
#
# Usage (source at the top of a notebook or analysis script):
#   data_pull_date       <- "20250805"   # optional — defaults set below
#   variant              <- "residuals"
#   tfbpmodeling_version <- "1.0.0"
#   source(here("R/prepare_data.R"))
#
# Objects created in the calling environment:
#   input_data               list(predictors_list, responses_list) — transformed
#   rr_meta                  tibble — rank-response metadata
#   predictors_meta          tibble — brent_nf_cc sample metadata
#   response_meta            tibble — mcisaac_oe sample metadata
#   cc_usable_meta           tibble — calling-cards replicate metadata
#   cc_usable_data           tibble — calling-cards replicate data (long)
#   preperturbation_expr     tibble — red_median_wide expression matrix
#   mcisaac_kem_predictor_rr tibble — rr_meta filtered to usable pTF/response pairs
#   ci_df_all_data           tibble — Stage 1 (all-data) CI intervals
#   ci_df_topn               tibble — Stage 2 (top-N) CI intervals
#   bootstrap_coefs_all_data named list — Stage 1 bootstrap coef draws per TF
#   bootstrap_coefs_topn     named list — Stage 2 bootstrap coef draws per TF
#   stage3_results_df        tibble — final-stage interactor significance results
#   stage_comp_df            tibble — interactor counts per TF across all stages
#                              (cols: regulator_symbol, all_data_n_interactors,
#                               topn_n_interactors, stage3_lassocv_n_interactors,
#                               stage3_bootstrap_n_interactors)

library(here)
library(tidyverse)
library(glue)

# =============================================================================
# Parameters
# =============================================================================
# if (!exists("data_pull_date"))       data_pull_date       <- "20250805"
# if (!exists("variant"))              variant              <- "residuals"
# if (!exists("tfbpmodeling_version")) tfbpmodeling_version <- "1.0.0"

data_pull_date       <- "20250805"
variant              <- "residuals"
tfbpmodeling_version <- "1.0.0"

input_dir   <- here("data", data_pull_date)
results_dir <- here("data", data_pull_date, tfbpmodeling_version, variant)

if (!dir.exists(input_dir)) {
    stop(glue(
        "Input data directory not found: {input_dir}\n",
        "Run `Rscript R/prepare_input_data.R` first."
    ))
}
if (!dir.exists(results_dir)) {
    stop(glue(
        "Results directory not found: {results_dir}\n",
        "Run `Rscript R/prepare_results_data.R` (or prepare_results_data.R) first.\n",
        "Expected path: data/{data_pull_date}/{tfbpmodeling_version}/{variant}/"
    ))
}

# =============================================================================
# Input data: predictors and responses (load raw, apply transforms)
# =============================================================================
message("Loading predictors and responses...")

theoretical_max <- log10(6066)

predictors_list <- readRDS(file.path(input_dir, "predictors_list.rds"))
responses_list  <- readRDS(file.path(input_dir, "responses_list.rds"))

# ---- predictor transforms --------------------------------------------------
predictors_list$rank <- predictors_list$rank %>%
    pivot_longer(-target_symbol, names_to = "regulator", values_to = "rank") %>%
    group_by(regulator) %>%
    mutate(rank = -rank + max(rank)) %>%
    pivot_wider(names_from = regulator, values_from = rank)

predictors_list$raw_pvalue <- predictors_list$raw_pvalue %>%
    pivot_longer(-target_symbol, names_to = "regulator", values_to = "raw_pvalue") %>%
    group_by(regulator) %>%
    mutate(raw_pvalue = -log10(raw_pvalue + 1e-300)) %>%
    pivot_wider(names_from = regulator, values_from = raw_pvalue)

predictors_list$log_rank_standardized <- predictors_list$log_rank %>%
    mutate(across(where(is.numeric), ~ .x / sd(.x, na.rm = TRUE)))

predictors_list$log_rank_normalized <- predictors_list$log_rank %>%
    mutate(across(where(is.numeric), ~ .x / theoretical_max))

# ---- response transforms ---------------------------------------------------
responses_list$raw <- responses_list$raw %>%
    mutate(across(where(is.numeric), abs))

responses_list$rank <- responses_list$rank %>%
    pivot_longer(-target_symbol, names_to = "regulator", values_to = "rank") %>%
    group_by(regulator) %>%
    mutate(rank = -rank + max(rank)) %>%
    pivot_wider(names_from = regulator, values_from = rank)

responses_list$log_rank_normalized <- responses_list$log_rank %>%
    mutate(across(
        where(is.numeric) & !matches("target_symbol"),
        ~ .x / theoretical_max
    ))

responses_list$log_rank_standardized <- responses_list$log_rank %>%
    mutate(across(where(is.numeric), ~ .x / sd(.x, na.rm = TRUE)))

input_data <- list(
    predictors_list = predictors_list,
    responses_list  = responses_list
)

message("  Done.")

# =============================================================================
# Metadata
# =============================================================================
message("Loading metadata...")

rr_meta              <- read_csv(file.path(input_dir, "rr_meta.csv"),          show_col_types = FALSE)
predictors_meta      <- read_csv(file.path(input_dir, "brent_nf_cc_meta.csv"), show_col_types = FALSE)
response_meta        <- read_csv(file.path(input_dir, "mcisaac_oe_meta.csv"),  show_col_types = FALSE)
cc_usable_meta       <- read_csv(file.path(input_dir, "cc_usable_meta.csv"),   show_col_types = FALSE)
cc_usable_data       <- readRDS(file.path(input_dir, "cc_usable_data.rds"))
preperturbation_expr <- read_csv(file.path(input_dir, "red_median_wide.csv"),  show_col_types = FALSE)

# Rank-response pairs for usable pTF / expression combinations
mcisaac_kem_predictor_rr <- rr_meta %>%
    filter(
        promotersetsig %in% predictors_meta$id,
        expression_source == "kemmeren_tfko" |
            (expression_source == "mcisaac_oe" & expression %in% response_meta$id)
    )

message("  Done.")

# =============================================================================
# Results: CI intervals and bootstrap coefs
# =============================================================================
message(glue("Loading results from data/{data_pull_date}/{tfbpmodeling_version}/{variant}/..."))

ci_df_all_data           <- readRDS(file.path(results_dir, "ci_df_all_data.rds"))
ci_df_topn               <- readRDS(file.path(results_dir, "ci_df_topn.rds"))
bootstrap_coefs_all_data <- readRDS(file.path(results_dir, "bootstrap_coefs_all_data.rds"))
bootstrap_coefs_topn     <- readRDS(file.path(results_dir, "bootstrap_coefs_topn.rds"))

# Stage results: v1 produces both stage3_lassocv_results.rds and
# ci_df_stage3_bootstrap_lassocv.rds; legacy produces only stage3_lassocv_results.rds.
ci_stage3_bootstrap_lassocv_path <- file.path(results_dir, "ci_df_stage3_bootstrap_lassocv.rds")
stage_results_lassocv_path       <- file.path(results_dir, "stage3_lassocv_results.rds")

if (file.exists(ci_stage3_bootstrap_lassocv_path)) {
    ci_df_stage3_bootstrap <- readRDS(ci_stage3_bootstrap_lassocv_path)
    message("  Loaded ci_df_stage3_bootstrap_lassocv.rds")
    message("  Loaded stage3_lassocv_results.rds")
} else {
    message("No stage3 results discovered for {variant}")
}

stage3_results_df <- readRDS(stage_results_lassocv_path)
message("  Loaded stage3_lassocv_results.rds (non-bootstrapped)")

message("  Done.")

# =============================================================================
# Derived objects (cheap — no model fitting)
# =============================================================================
message("Building derived objects...")

# Interactor counts per TF for Stage 1 and Stage 2
significant_all_data_df <- ci_df_all_data %>%
    count(regulator, name = "all_data_n_interactors") %>%
    dplyr::rename(regulator_symbol = regulator)

significant_topn_df <- ci_df_topn %>%
    count(regulator, name = "topn_n_interactors") %>%
    dplyr::rename(regulator_symbol = regulator)

# Flag consistent interactors: coef sign in stage_results matches topn CI sign
topn_signs <- ci_df_topn %>%
    select(model = regulator, interactor, topn_sign = sign)

stage3_consistent_df <- stage3_results_df %>%
    filter(!is.na(coef_interactor)) %>%
    left_join(topn_signs, by = c("model", "interactor")) %>%
    mutate(consistent = sign(coef_interactor) == topn_sign) %>%
    filter(consistent, complete.cases(.))

surviving_stage3_lassocv_df = stage3_consistent_df %>%
    count(model, name = "stage3_lassocv_n_interactors") %>%
    dplyr::rename(regulator_symbol = model)

surviving_stage3_bootstrap_df = tryCatch(
    ci_df_stage3_bootstrap %>%
        mutate(interaction_term = str_detect(interactor, ":")) %>%
        filter(interaction_term) %>%
        count(regulator, name = "stage3_bootstrap_n_interactors") %>%
        dplyr::rename(regulator_symbol = regulator),
    error = function(e) tibble(regulator_symbol = character(), stage3_bootstrap_n_interactors = integer())
)

all_model_names <- tibble(
    regulator_symbol = intersect(colnames(input_data$responses_list$log_rank_standardized),
                                 colnames((input_data$predictors_list$raw_pvalue)))[-1]
)

stage_comp_df <- all_model_names %>%
    left_join(significant_topn_df) %>%
    left_join(significant_all_data_df) %>%
    left_join(surviving_stage3_bootstrap_df) %>%
    left_join(surviving_stage3_lassocv_df) %>%
    replace_na(list(all_data_n_interactors = 0L,
                    topn_n_interactors = 0L,
                    stage3_lassocv_n_interactors = 0L,
                    stage3_bootstrap_n_interactors = 0L)) %>%
    dplyr::relocate(regulator_symbol,
                    all_data_n_interactors,
                    topn_n_interactors,
                    stage3_lassocv_n_interactors,
                    stage3_bootstrap_n_interactors) %>%
    arrange(desc(all_data_n_interactors))

message("  Done.")

# =============================================================================
message("=== prepare_data.R complete ===")
message(glue("Data loaded from data/{data_pull_date}/ and data/{data_pull_date}/{tfbpmodeling_version}/{variant}/"))
