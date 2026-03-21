# compare_results.R
#
# Load and compare two versions of per-regulator stage results from the
# tfbpmodeling pipeline.
#
# ---------------------------------------------------------------------------
# Input structure for load_stage_results():
#
#   A named list keyed by regulator symbol. Each entry must contain exactly
#   two sublists (names are arbitrary — first = version A, second = version B).
#   Each sublist must contain exactly three file paths IN ORDER:
#     [[1]]  all_data_significant_98-0.json              (Stage 1)
#     [[2]]  topn_significant_90-0.json                  (Stage 2)
#     [[3]]  stage3_lassocv_significance_results.json    (Stage 3)
#
#   Example:
#     list(
#       AFT1 = list(
#         old = list("path/v0/AFT1/all_data...json",
#                    "path/v0/AFT1/topn...json",
#                    "path/v0/AFT1/stage3...json"),
#         new = list("path/v1/AFT1/all_data...json",
#                    "path/v1/AFT1/topn...json",
#                    "path/v1/AFT1/stage3...json")
#       )
#     )
#
# ---------------------------------------------------------------------------
# Output of load_stage_results():
#
#   Same top-level structure (keyed by regulator), but each version sublist
#   is now named list(stage1 = <parsed>, stage2 = <parsed>, stage3 = <parsed>).
#
# ---------------------------------------------------------------------------
# Comparison functions:
#
#   compare_stage1(a, b)  — Stage 1 (all-data 98% CI) interactor sets
#   compare_stage2(a, b)  — Stage 2 (top-N 90% CI) interactor sets
#   compare_stage3(a, b)  — Stage 3 (lasso-CV) interactor sets + coef sign
#
#   Each takes the stage-N parsed JSON for version A and version B for a
#   single regulator and returns a tibble summarising interactor-level
#   agreement.

library(tidyverse)
library(jsonlite)

# ---------------------------------------------------------------------------
# load_stage_results()
# ---------------------------------------------------------------------------
load_stage_results <- function(path_list) {
    stage_names <- c("stage1", "stage2", "stage3")

    imap(path_list, function(reg_entry, reg_name) {
        if (length(reg_entry) != 2)
            stop(glue::glue(
                "Regulator '{reg_name}': expected 2 version sublists, got {length(reg_entry)}"
            ))

        imap(reg_entry, function(version_paths, version_name) {
            if (length(version_paths) != 3)
                stop(glue::glue(
                    "Regulator '{reg_name}', version '{version_name}': ",
                    "expected 3 file paths, got {length(version_paths)}"
                ))

            setNames(
                map(version_paths, function(p) {
                    p_expanded <- path.expand(p)
                    if (!file.exists(p_expanded)) {
                        message(glue::glue("File not found (returning empty result): {p_expanded}"))
                        return(list())
                    }
                    jsonlite::read_json(p_expanded)
                }),
                stage_names
            )
        })
    })
}

# ---------------------------------------------------------------------------
# Internal: shared set-comparison logic for stage1 and stage2
# ---------------------------------------------------------------------------
.compare_interactor_sets <- function(a, b) {
    interactors_a   <- as.character(names(a))
    interactors_b   <- as.character(names(b))
    all_interactors <- union(interactors_a, interactors_b)

    tibble(
        interactor = all_interactors,
        in_a       = all_interactors %in% interactors_a,
        in_b       = all_interactors %in% interactors_b
    ) %>%
        mutate(
            agreement = case_when(
                in_a & in_b   ~ "both",
                in_a & !in_b  ~ "only_a",
                !in_a & in_b  ~ "only_b"
            )
        ) %>%
        arrange(agreement, interactor)
}

# ---------------------------------------------------------------------------
# compare_stage1()  — Stage 1 (all-data 98% CI)
#
#   a, b: parsed stage1 JSON for a single regulator (named list keyed by
#         interactor, values are CI objects)
# ---------------------------------------------------------------------------
compare_stage1 <- function(a, b) {
    .compare_interactor_sets(a, b)
}

# ---------------------------------------------------------------------------
# compare_stage2()  — Stage 2 (top-N 90% CI)
#
#   a, b: parsed stage2 JSON for a single regulator
# ---------------------------------------------------------------------------
compare_stage2 <- function(a, b) {
    .compare_interactor_sets(a, b)
}

# ---------------------------------------------------------------------------
# compare_stage3()  — Stage 3 (lasso-CV significance results)
#
#   a, b: parsed stage3 JSON for a single regulator (named list keyed by
#         interactor; each entry expected to have $coef_interactor)
#
#   In addition to set agreement, compares the sign of coef_interactor for
#   interactors present in both versions.
# ---------------------------------------------------------------------------
compare_stage3 <- function(a, b) {
    base_df <- .compare_interactor_sets(a, b)

    both <- base_df %>% filter(agreement == "both") %>% pull(interactor)

    coef_df <- tibble(
        interactor   = both,
        coef_a       = map_dbl(both, ~ {
            val <- a[[.x]]$coef_interactor
            if (is.null(val)) NA_real_ else as.numeric(val)
        }),
        coef_b       = map_dbl(both, ~ {
            val <- b[[.x]]$coef_interactor
            if (is.null(val)) NA_real_ else as.numeric(val)
        })
    ) %>%
        mutate(sign_agree = case_when(
            is.na(coef_a) | is.na(coef_b) ~ NA,
            sign(coef_a) == sign(coef_b)   ~ TRUE,
            TRUE                            ~ FALSE
        ))

    base_df %>%
        left_join(coef_df, by = "interactor")
}
