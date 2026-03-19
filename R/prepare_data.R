# prepare_data.R
#
# One-time data extraction script. Reads from the htcf_ref mount and external
# Python APIs; serializes all data into data/ at the repo root so that figure
# scripts can run without the remote mount or a Python environment.
#
# Run via:
#   Rscript R/prepare_data.R
#
# Or source with custom params:
#   date    <- "20250805"
#   variant <- "residuals"
#   source(here("R/prepare_data.R"))
#
# Produces (in data/):
#   Section 1 (pure R):       predictors_list_{date}.rds, responses_list_{date}.rds
#   Section 2 (pure R):       ci_df_all_data_{date}.rds, ci_df_topn_{date}.rds
#   Section 3 (pure R):       stage4_results_{variant}_{date}.rds
#   Section 4 (reticulate):   bootstrap_coefs_all_data_{date}.rds,
#                             bootstrap_coefs_topn_{date}.rds
#   Section 5 (reticulate):   rr_meta_{date}.csv, brent_nf_cc_meta_{date}.csv,
#                             mcisaac_oe_meta_{date}.csv
#   Section 6 (reticulate):   cc_usable_meta_{date}.csv, cc_usable_data_{date}.rds
#   Section 7 (pure R):       red_median_wide_{date}.csv

library(here)
library(tidyverse)
library(glue)
library(jsonlite)

results_basepath <- "/home/chase/htcf_ref/data/yeast_database_modelling/results"
results_variant  <- glue("results_{date}_{variant}")
results_dir      <- file.path(results_basepath, date, results_variant)
htcf_data_dir    <- glue("~/htcf_ref/data/yeast_database_modelling/pull_data_{date}/data")

dir.create(here("data"), showWarnings = FALSE)

# =============================================================================
# Run directly (source these lines first, then source the script)
# =============================================================================
date    <- "20250805"
variant <- "residuals"
pull_data <- FALSE # or TRUE

# =============================================================================
# create_predictors_response_lists
# =============================================================================

library(tidyverse)
library(here)
library(glue)

# eg pull_data = FALSE, date = "20250801"
pull_predictor_response_lists <- function(pull_data, date) {
    theoretical_max <- log10(6066)
    input_basedir <- glue("~/htcf_ref/data/yeast_database_modelling/pull_data_{date}")

    if (pull_data) {
        predictors_list <- list(
            raw_pvalue = read_csv(file.path(input_basedir, glue("data/predictors_brent_nf_cc_mcisaac_oe_raw_pvalue_{date}.csv"))),
            raw_enrichment = read_csv(file.path(input_basedir, glue("data/predictors_brent_nf_cc_mcisaac_oe_raw_enrichment_{date}.csv"))),
            rank = read_csv(file.path(input_basedir, glue("data/predictors_brent_nf_cc_mcisaac_oe_rank_{date}.csv"))),
            log_rank = read_csv(file.path(input_basedir, glue("data/predictors_brent_nf_cc_mcisaac_oe_{date}.csv")))
        )
        saveRDS(predictors_list, here(glue("data/predictors_list_{date}.rds")))

        responses_list <- list(
            raw = read_csv(file.path(input_basedir, glue("data/response_brent_nf_cc_mcisaac_oe_raw_{date}.csv"))),
            rank = read_csv(file.path(input_basedir, glue("data/response_brent_nf_cc_mcisaac_oe_abs_{date}.csv"))),
            log_rank = read_csv(file.path(input_basedir, glue("data/response_brent_nf_cc_mcisaac_oe_{date}.csv")))
        )

        norm_resid_paths <- list.files(
            glue("~/htcf_ref/data/yeast_database_modelling/pull_data_{date}/response_frames_residuals_normalized"),
            full.names = TRUE
        )
        names(norm_resid_paths) <- str_remove(basename(norm_resid_paths), ".csv")
        norm_resid_dfs <- map(norm_resid_paths, read_csv)
        responses_list$residuals_normalized <- reduce(norm_resid_dfs, left_join, by = "target_symbol")

        saveRDS(responses_list, here(glue("data/responses_list_{date}.rds")))
    } else {
        predictors_list <- readRDS(here(glue("data/predictors_list_{date}.rds")))
        responses_list <- readRDS(here(glue("data/responses_list_{date}.rds")))
    }

    # Transform predictors
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

    # Transform responses
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

    list(
        predictors_list = predictors_list,
        responses_list = responses_list
    )
}

# =============================================================================
# Section 1: Predictor / Response Data (pure R)
# =============================================================================
message("=== Section 1: predictor/response data ===")
pull_predictor_response_lists(pull_data = pull_data, date = date)
message("  Saved: data/predictors_list_{date}.rds, data/responses_list_{date}.rds")

# =============================================================================
# Section 2: CI Data from result_obj.json (pure R)
#
# result_obj.json format: {"98.0": {"ACA1:OPI1": [-0.066, -0.027], ...}}
# Only significant interactions are stored — each entry is already filtered.
# =============================================================================
message("=== Section 2: CI data from result_obj.json ===")

parse_ci_json <- function(result_dir_name, ci_level) {
    json_paths <- list.files(
        results_dir,
        pattern     = "result_obj.json",
        recursive   = TRUE,
        full.names  = TRUE
    )
    json_paths <- json_paths[grepl(result_dir_name, json_paths, fixed = TRUE)]

    map_dfr(json_paths, function(p) {
        regulator <- basename(dirname(dirname(p)))
        ci_data   <- read_json(p)[[ci_level]]
        if (is.null(ci_data) || length(ci_data) == 0) return(NULL)
        tibble(
            regulator  = regulator,
            interactor = names(ci_data),
            ci_lo      = map_dbl(ci_data, 1),
            ci_hi      = map_dbl(ci_data, 2)
        )
    }) %>%
        mutate(sign = if_else(ci_lo > 0, 1L, -1L))
}

ci_df_all_data <- parse_ci_json("all_data_result_object", "98.0")
ci_df_topn     <- parse_ci_json("topn_result_object",     "90.0")

saveRDS(ci_df_all_data, here(glue("data/ci_df_all_data_{date}.rds")))
saveRDS(ci_df_topn,     here(glue("data/ci_df_topn_{date}.rds")))
message(glue("  Saved: data/ci_df_all_data_{date}.rds ({nrow(ci_df_all_data)} rows)"))
message(glue("  Saved: data/ci_df_topn_{date}.rds ({nrow(ci_df_topn)} rows)"))

# =============================================================================
# Section 3: Stage 4 Results from interactor_vs_main_result.json (pure R)
# =============================================================================
message("=== Section 3: stage4 results ===")

models <- list.files(results_dir, include.dirs = TRUE, full.names = FALSE)

stage4_list <- map(models, function(model) {
    json_path <- file.path(results_dir, model, "interactor_vs_main_result.json")
    if (file.exists(json_path)) {
        df        <- read_json(json_path, simplifyVector = TRUE)
        df$model  <- model
    } else {
        df <- data.frame(
            interactor       = NA_character_,
            variant          = NA_character_,
            r2_lasso_model   = NA_real_,
            coef_interactor  = NA_real_,
            coef_main_effect = NA_real_,
            model            = model
        )
    }
    df
})

stage4_df <- bind_rows(stage4_list) %>% as_tibble()

saveRDS(stage4_df, here(glue("data/stage4_results_{variant}_{date}.rds")))
message(glue("  Saved: data/stage4_results_{variant}_{date}.rds ({nrow(stage4_df)} rows)"))

# =============================================================================
# Section 4: Bootstrap Coefficient Draws (reticulate — requires .venv + pkl)
# =============================================================================
message("=== Section 4: bootstrap coefficient draws (requires Python) ===")

library(reticulate)
use_virtualenv(here(".venv"), required = TRUE)

py_run_string(glue("
import pickle
from pathlib import Path
import pandas as pd

results_dir = '{results_dir}'
models      = [p.name for p in Path(results_dir).iterdir() if p.is_dir()]

def load_coefs(result_dir_name):
    out = {{}}
    for reg in models:
        pkl_path = Path(results_dir, reg, result_dir_name, 'result_obj.pkl')
        if pkl_path.exists():
            with open(pkl_path, 'rb') as f:
                bootstrap_coefs_df, alpha_list = pickle.load(f)
            out[reg] = bootstrap_coefs_df
    return out

all_data_coefs_py = load_coefs('all_data_result_object')
topn_coefs_py     = load_coefs('topn_result_object')
"))

all_data_coefs <- map(py$all_data_coefs_py, as_tibble)
topn_coefs     <- map(py$topn_coefs_py, as_tibble)

saveRDS(all_data_coefs, here(glue("data/bootstrap_coefs_all_data_{date}.rds")))
saveRDS(topn_coefs,     here(glue("data/bootstrap_coefs_topn_{date}.rds")))
message(glue("  Saved: data/bootstrap_coefs_all_data_{date}.rds ({length(all_data_coefs)} regulators)"))
message(glue("  Saved: data/bootstrap_coefs_topn_{date}.rds ({length(topn_coefs)} regulators)"))

# =============================================================================
# Section 5: Rank-Response Metadata (reticulate — tfbpapi)
# =============================================================================
message("=== Section 5: rank-response metadata (requires tfbpapi) ===")

py_run_string("
import asyncio
from tfbpapi import RankResponseAPI
loop = asyncio.get_event_loop()
rr_api = RankResponseAPI()
rr_res  = loop.run_until_complete(rr_api.read())
")

rr_meta <- py$rr_res$metadata %>% as_tibble()
write_csv(rr_meta, here(glue("data/rr_meta_{date}.csv")))
message(glue("  Saved: data/rr_meta_{date}.csv ({nrow(rr_meta)} rows)"))

# Copy predictor and response metadata from htcf_ref
file.copy(
    file.path(htcf_data_dir, glue("brent_nf_cc_meta_{date}.csv")),
    here(glue("data/brent_nf_cc_meta_{date}.csv")),
    overwrite = TRUE
)
file.copy(
    file.path(htcf_data_dir, glue("mcisaac_oe_meta_{date}.csv")),
    here(glue("data/mcisaac_oe_meta_{date}.csv")),
    overwrite = TRUE
)
message(glue("  Copied: data/brent_nf_cc_meta_{date}.csv, data/mcisaac_oe_meta_{date}.csv"))

# =============================================================================
# Section 6: Calling-Cards Replicate Data (reticulate — tfbpapi)
# =============================================================================
message("=== Section 6: calling-cards replicate data (requires tfbpapi) ===")

py_run_string("
import asyncio
from tfbpapi import PromoterSetSigAPI
loop    = asyncio.get_event_loop()
pss_api = PromoterSetSigAPI()
pss_api.pop_params()
pss_api.push_params({'source_name': 'brent_nf_cc', 'data_usable': 'pass'})
cc_usable_res = loop.run_until_complete(pss_api.read(retrieve_files=True))
")

cc_usable_meta <- py$cc_usable_res$metadata %>%
    as_tibble() %>%
    filter(!is.na(single_binding))

cc_usable_data <- map(cc_usable_meta$id, ~ {
    py$cc_usable_res$data[[as.character(.)]] %>% as_tibble()
}) %>% bind_rows()

write_csv(cc_usable_meta, here(glue("data/cc_usable_meta_{date}.csv")))
saveRDS(cc_usable_data,   here(glue("data/cc_usable_data_{date}.rds")))
message(glue("  Saved: data/cc_usable_meta_{date}.csv ({nrow(cc_usable_meta)} replicates)"))
message(glue("  Saved: data/cc_usable_data_{date}.rds ({nrow(cc_usable_data)} rows)"))

# =============================================================================
# Section 7: Pre-Perturbation Expression (pure R — file copy)
# =============================================================================
message("=== Section 7: pre-perturbation expression ===")

file.copy(
    file.path(htcf_data_dir, "red_median_wide.csv"),
    here(glue("data/red_median_wide_{date}.csv")),
    overwrite = TRUE
)
message(glue("  Copied: data/red_median_wide_{date}.csv"))

# =============================================================================
message("=== prepare_data.R complete ===")
message("All data objects written to data/")
