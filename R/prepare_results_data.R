# prepare_results_data.R
#
# Extracts and serializes results data objects from tfbpmodeling v1 output,
# which introduces a stage3_lassocv_bootstrap stage as the final stage
# (replacing the old interactor_vs_main stage3). Reads from a local results
# directory and writes parsed RDS files into
# data/<data_pull_date>/<tfbpmodeling_version>/. Section 4 requires the .venv
# Python virtualenv (pickle).
#
# `data_pull_date` is the date string identifying which versioned data pull from
# the yeast database to use (e.g. "20250805"). It determines the output
# subdirectory where parsed results are written.
#
# `tfbpmodeling_version` is the version of the tfbpmodeling software that was
# used to produce the results (e.g. "1.0.0"). This determines the output
# subdirectory so that results from different software versions do not collide.
#
# `variant` is the data variant (e.g. "residuals"). It becomes a subdirectory
# under the version dir so filenames within each variant dir are standardized
# (no _{variant} suffix needed).
#
# Run via:
#   Rscript R/prepare_results_data.R
#
# Or source with custom params:
#   data_pull_date       <- "20250805"
#   variant              <- "residuals"
#   tfbpmodeling_version <- "1.0.0"
#   results_src_dir      <- here("data/v1_0_0_stage3_bootstrap_lassocv_results")
#   source(here("R/prepare_results_data.R"))
#
# Produces (in data/<data_pull_date>/<tfbpmodeling_version>/<variant>/):
#   Section 2:                ci_df_all_data.rds, ci_df_topn.rds,
#                             ci_df_stage3_bootstrap_lassocv.rds
#   Section 3:                stage3_lassocv_results.rds
#   Section 4 (reticulate):   bootstrap_coefs_all_data.rds,
#                             bootstrap_coefs_topn.rds,
#                             bootstrap_coefs_stage3_bootstrap_lassocv.rds

library(here)
library(tidyverse)
library(glue)
library(jsonlite)

# =============================================================================
# Parameters (set these before sourcing, or edit defaults here)
# =============================================================================
if (!exists("data_pull_date"))       data_pull_date       <- "20250805"
if (!exists("variant"))              variant              <- "residuals"
if (!exists("tfbpmodeling_version")) tfbpmodeling_version <- "1.1.0"
if (!exists("results_src_dir"))      results_src_dir      <- here("data/v1.1.0/results_residuals")

# Expand ~ so paths are absolute when passed to Python
results_src_dir <- path.expand(results_src_dir)

stopifnot(dir.exists(results_src_dir))

data_dir <- here("data", data_pull_date, tfbpmodeling_version, variant)

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Section 2: CI Data from pre-computed significant JSON files (pure R)
#
# File format: {"TF:interactor": [lo, hi], ...}
# One file per TF per stage, e.g. all_data_significant_98-0.json
# =============================================================================
message("=== Section 2: CI data from significant JSON files ===")

parse_ci_json <- function(results_dir, filename_pattern) {
    json_paths <- list.files(
        results_dir,
        pattern    = filename_pattern,
        recursive  = TRUE,
        full.names = TRUE
    )

    result <- map_dfr(json_paths, function(p) {
        regulator <- basename(dirname(p))
        ci_data   <- read_json(p)
        if (length(ci_data) == 0) return(NULL)
        tibble(
            regulator  = regulator,
            interactor = names(ci_data),
            ci_lo      = map_dbl(ci_data, 1),
            ci_hi      = map_dbl(ci_data, 2)
        )
    })
    if (nrow(result) == 0) return(tibble(regulator = character(), interactor = character(),
                                         ci_lo = double(), ci_hi = double(), sign = integer()))
    mutate(result, sign = if_else(ci_lo > 0, 1L, -1L))
}

ci_df_all_data       <- parse_ci_json(results_src_dir, "^all_data_significant_98-0\\.json$")
ci_df_topn           <- parse_ci_json(results_src_dir, "^topn_significant_90-0\\.json$")
ci_df_stage3_bootstrap_lassocv <- parse_ci_json(results_src_dir, "^stage3_lassocv_bootstrap_significant_98-0\\.json$")

saveRDS(ci_df_all_data,       file.path(data_dir, "ci_df_all_data.rds"))
saveRDS(ci_df_topn,           file.path(data_dir, "ci_df_topn.rds"))
saveRDS(ci_df_stage3_bootstrap_lassocv, file.path(data_dir, "ci_df_stage3_bootstrap_lassocv.rds"))
message(glue("  Saved: data/{data_pull_date}/{tfbpmodeling_version}/{variant}/ci_df_all_data.rds ({nrow(ci_df_all_data)} rows)"))
message(glue("  Saved: data/{data_pull_date}/{tfbpmodeling_version}/{variant}/ci_df_topn.rds ({nrow(ci_df_topn)} rows)"))
message(glue("  Saved: data/{data_pull_date}/{tfbpmodeling_version}/{variant}/ci_df_stage3_bootstrap_lassocv.rds ({nrow(ci_df_stage3_bootstrap_lassocv)} rows)"))

# =============================================================================
# Section 3: Stage 3 LassoCV Significance Results
#
# stage3_lassocv_significance_results.json: array of objects with fields:
#   interactor, variant, r2_lasso_model, coef_interactor, coef_main_effect
# This replaces the old interactor_vs_main_result.json (stage3).
# =============================================================================
message("=== Section 3: stage3 lassocv significance results ===")

tfs <- list.files(results_src_dir, include.dirs = TRUE, full.names = FALSE)

stage3_lassocv_list <- map(tfs, function(tf) {
    json_path <- file.path(results_src_dir, tf, "stage3_lassocv_significance_results.json")
    if (file.exists(json_path)) {
        df       <- read_json(json_path, simplifyVector = TRUE)
        df$model <- tf
    } else {
        df <- data.frame(
            interactor       = NA_character_,
            variant          = NA_character_,
            r2_lasso_model   = NA_real_,
            coef_interactor  = NA_real_,
            coef_main_effect = NA_real_,
            model            = tf
        )
    }
    df
})

stage3_lassocv_df <- bind_rows(stage3_lassocv_list) %>% as_tibble()

saveRDS(stage3_lassocv_df, file.path(data_dir, "stage3_lassocv_results.rds"))
message(glue("  Saved: data/{data_pull_date}/{tfbpmodeling_version}/{variant}/stage3_lassocv_results.rds ({nrow(stage3_lassocv_df)} rows)"))

# =============================================================================
# Section 4: Bootstrap Coefficient Draws (reticulate — requires .venv + pkl)
#
# Loads pkl from all three stage result object directories.
# Each pkl is a tuple: (bootstrap_coefs_df, alpha_list)
# =============================================================================
message("=== Section 4: bootstrap coefficient draws (requires Python) ===")

library(reticulate)
use_virtualenv(here(".venv"), required = TRUE)

py_run_string(glue("
import pickle
from pathlib import Path
import pandas as pd

results_src_dir = '{results_src_dir}'
tfs = [p.name for p in Path(results_src_dir).iterdir() if p.is_dir()]

def load_coefs(result_dir_name):
    out = {{}}
    for tf in tfs:
        pkl_path = Path(results_src_dir, tf, result_dir_name, 'result_obj.pkl')
        if pkl_path.exists():
            with open(pkl_path, 'rb') as f:
                bootstrap_coefs_df, alpha_list = pickle.load(f)
            out[tf] = bootstrap_coefs_df
    return out

all_data_coefs_py       = load_coefs('all_data_result_object')
topn_coefs_py           = load_coefs('topn_result_object')
stage3_bootstrap_lassocv_coefs_py = load_coefs('stage3_lassocv_bootstrap_result_object')
"))

all_data_coefs       <- map(py$all_data_coefs_py, as_tibble)
topn_coefs           <- map(py$topn_coefs_py, as_tibble)
stage3_bootstrap_lassocv_coefs <- map(py$stage3_bootstrap_lassocv_coefs_py, as_tibble)

saveRDS(all_data_coefs,       file.path(data_dir, "bootstrap_coefs_all_data.rds"))
saveRDS(topn_coefs,           file.path(data_dir, "bootstrap_coefs_topn.rds"))
saveRDS(stage3_bootstrap_lassocv_coefs, file.path(data_dir, "bootstrap_coefs_stage3_bootstrap_lassocv.rds"))
message(glue("  Saved: data/{data_pull_date}/{tfbpmodeling_version}/{variant}/bootstrap_coefs_all_data.rds ({length(all_data_coefs)} TFs)"))
message(glue("  Saved: data/{data_pull_date}/{tfbpmodeling_version}/{variant}/bootstrap_coefs_topn.rds ({length(topn_coefs)} TFs)"))
message(glue("  Saved: data/{data_pull_date}/{tfbpmodeling_version}/{variant}/bootstrap_coefs_stage3_bootstrap_lassocv.rds ({length(stage3_bootstrap_lassocv_coefs)} TFs)"))

# =============================================================================
message("=== prepare_results_data.R complete ===")
message(glue("All results data objects written to data/{data_pull_date}/{tfbpmodeling_version}/{variant}/"))
