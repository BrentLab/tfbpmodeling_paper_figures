# prepare_input_data.R
#
# One-time pull script. Reads raw CSVs from the ~/htcf_ref mount and external
# Python APIs and saves them verbatim into data/<data_pull_date>/. No transforms
# are applied here — all transforms live in prepare_data.R.
# Sections 5 and 6 require the .venv Python virtualenv (tfbpapi).
#
# `data_pull_date` is the date string identifying which versioned data pull from
# the yeast database to use (e.g. "20250805"). Files must exist at:
#   ~/htcf_ref/data/yeast_database_modelling/pull_data_<data_pull_date>/
#
# Run via:
#   Rscript R/prepare_input_data.R
#
# Or source with custom params:
#   data_pull_date <- "20250805"
#   source(here("R/prepare_input_data.R"))
#
# Produces (in data/<data_pull_date>/):
#   Section 1 (pure R):       predictors_list.rds, responses_list.rds
#   Section 5 (reticulate):   rr_meta.csv, brent_nf_cc_meta.csv,
#                             mcisaac_oe_meta.csv
#   Section 6 (reticulate):   cc_usable_meta.csv, cc_usable_data.rds
#   Section 7 (pure R):       red_median_wide.csv

library(here)
library(tidyverse)
library(glue)

# =============================================================================
# Parameters (set these before sourcing, or edit defaults here)
# =============================================================================
if (!exists("data_pull_date")) data_pull_date <- "20250805"

# data_pull_date <- "20250805"

htcf_data_dir <- glue("~/htcf_ref/data/yeast_database_modelling/pull_data_{data_pull_date}/data")
data_dir      <- here("data", data_pull_date)

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Section 1: Predictor / Response Data (pure R — pull from htcf_ref)
# =============================================================================
message("=== Section 1: predictor/response data ===")

input_basedir <- glue("~/htcf_ref/data/yeast_database_modelling/pull_data_{data_pull_date}")

predictors_list <- list(
    raw_pvalue = read_csv(file.path(input_basedir, glue("data/predictors_brent_nf_cc_mcisaac_oe_raw_pvalue_{data_pull_date}.csv"))),
    raw_enrichment = read_csv(file.path(input_basedir, glue("data/predictors_brent_nf_cc_mcisaac_oe_raw_enrichment_{data_pull_date}.csv"))),
    rank = read_csv(file.path(input_basedir, glue("data/predictors_brent_nf_cc_mcisaac_oe_rank_{data_pull_date}.csv"))),
    log_rank = read_csv(file.path(input_basedir, glue("data/predictors_brent_nf_cc_mcisaac_oe_{data_pull_date}.csv")))
)
saveRDS(predictors_list, file.path(data_dir, "predictors_list.rds"))

responses_list <- list(
    raw = read_csv(file.path(input_basedir, glue("data/response_brent_nf_cc_mcisaac_oe_raw_{data_pull_date}.csv"))),
    rank = read_csv(file.path(input_basedir, glue("data/response_brent_nf_cc_mcisaac_oe_abs_{data_pull_date}.csv"))),
    log_rank = read_csv(file.path(input_basedir, glue("data/response_brent_nf_cc_mcisaac_oe_{data_pull_date}.csv")))
)

norm_resid_paths <- list.files(
    glue("~/htcf_ref/data/yeast_database_modelling/pull_data_{data_pull_date}/response_frames_residuals_normalized"),
    full.names = TRUE
)
names(norm_resid_paths) <- str_remove(basename(norm_resid_paths), ".csv")
norm_resid_dfs <- map(norm_resid_paths, read_csv)
responses_list$residuals_normalized <- reduce(norm_resid_dfs, left_join, by = "target_symbol")

saveRDS(responses_list, file.path(data_dir, "responses_list.rds"))
message(glue("  Saved: data/{data_pull_date}/predictors_list.rds, data/{data_pull_date}/responses_list.rds"))

# =============================================================================
# Section 5: Rank-Response Metadata (reticulate — tfbpapi)
# =============================================================================
message("=== Section 5: rank-response metadata (requires tfbpapi) ===")

library(reticulate)
use_virtualenv(here(".venv"), required = TRUE)

py_run_string("
import asyncio
from tfbpapi import RankResponseAPI
loop = asyncio.get_event_loop()
rr_api = RankResponseAPI()
rr_res  = loop.run_until_complete(rr_api.read())
")

rr_meta <- py$rr_res$metadata %>% as_tibble()
write_csv(rr_meta, file.path(data_dir, "rr_meta.csv"))
message(glue("  Saved: data/{data_pull_date}/rr_meta.csv ({nrow(rr_meta)} rows)"))

file.copy(
    file.path(htcf_data_dir, glue("brent_nf_cc_meta_{data_pull_date}.csv")),
    file.path(data_dir, "brent_nf_cc_meta.csv"),
    overwrite = TRUE
)
file.copy(
    file.path(htcf_data_dir, glue("mcisaac_oe_meta_{data_pull_date}.csv")),
    file.path(data_dir, "mcisaac_oe_meta.csv"),
    overwrite = TRUE
)
message(glue("  Copied: data/{data_pull_date}/brent_nf_cc_meta.csv, data/{data_pull_date}/mcisaac_oe_meta.csv"))

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

write_csv(cc_usable_meta, file.path(data_dir, "cc_usable_meta.csv"))
saveRDS(cc_usable_data,   file.path(data_dir, "cc_usable_data.rds"))
message(glue("  Saved: data/{data_pull_date}/cc_usable_meta.csv ({nrow(cc_usable_meta)} replicates)"))
message(glue("  Saved: data/{data_pull_date}/cc_usable_data.rds ({nrow(cc_usable_data)} rows)"))

# =============================================================================
# Section 7: Pre-Perturbation Expression (pure R — file copy)
# =============================================================================
message("=== Section 7: pre-perturbation expression ===")

file.copy(
    file.path(htcf_data_dir, "red_median_wide.csv"),
    file.path(data_dir, "red_median_wide.csv"),
    overwrite = TRUE
)
message(glue("  Copied: data/{data_pull_date}/red_median_wide.csv"))

# =============================================================================
message("=== prepare_input_data.R complete ===")
message(glue("All raw input files written to data/{data_pull_date}/"))
