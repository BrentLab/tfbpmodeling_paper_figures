# fig1b1.R
#
# Boxplot comparing Pearson correlation *between replicates* (same TF, different
# CC experiments) vs *between TFs* (across-TF predictor correlation), both
# measured on the -log10(p-value) data transform.
#
# Usage (import into notebook):
#   source(here("R/fig1b1.R"))
#   plt <- make_fig1b1(long_corr_pearson_res, replicate_corr_comp_df)
#
# `long_corr_pearson_res` is built from input_data (see workflow notebook,
# lines ~163-174): upper-triangle Pearson correlations across TF predictors,
# long format with columns: var1, var2, corr, data_transform.
#
# `replicate_corr_comp_df` is built from calling-cards replicate data via
# the Python API (lines ~394-402): pairwise Pearson/Spearman correlations
# across CC replicates, with columns: regulator_symbol, record_id_1,
# record_id_2, method, corr, data_transform.

library(here)
library(tidyverse)
source(here("R/theme_ptf.R"))

# ---------------------------------------------------------------------------
# make_fig1b1()
#
# Parameters
# ----------
# long_corr_pearson_res   Long-format predictor Pearson correlation data frame.
# replicate_corr_comp_df  Long-format replicate Pearson correlation data frame.
# data_transform_label    Which data_transform level to use for both inputs
#                         (default: "-log10(pvalue)").
# base_size               Font size passed to theme_ptf().
#
# Returns a ggplot object.
# ---------------------------------------------------------------------------
make_fig1b1 <- function(
        long_corr_pearson_res,
        replicate_corr_comp_df,
        data_transform_label = "-log10(pvalue)",
        base_size = 18
) {
    predictor_corrs <- long_corr_pearson_res %>%
        mutate(
            data_transform = factor(
                data_transform,
                levels = c("enrichment", "pval", "lrb"),
                labels = c("Abs(enrichment)", "-log10(pvalue)", "LRB")
            )
        ) %>%
        filter(var1 != var2, data_transform == data_transform_label) %>%
        select(corr) %>%
        mutate(source = "predictors")

    replicate_corrs <- replicate_corr_comp_df %>%
        filter(method == "pearson", data_transform == data_transform_label) %>%
        select(corr) %>%
        mutate(source = "replicates")

    plot_df <- bind_rows(predictor_corrs, replicate_corrs) %>%
        mutate(
            source = factor(
                source,
                levels = c("replicates", "predictors"),
                labels = c("Between\nreplicates", "Between\nTFs")
            )
        )

    ggplot(plot_df, aes(source, corr)) +
        geom_boxplot(width = 0.9) +
        coord_cartesian(ylim = c(0, 1)) +
        ylab("Pearson Correlation") +
        theme_ptf(base_size = base_size) +
        theme(
            axis.title.x       = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_blank()
        )
}

# ---------------------------------------------------------------------------
# save_fig1b1()
# ---------------------------------------------------------------------------
save_fig1b1 <- function(path, long_corr_pearson_res, replicate_corr_comp_df,
                         width = 4, height = 4, bg = "white", ...) {
    plt <- make_fig1b1(long_corr_pearson_res, replicate_corr_comp_df, ...)
    ggsave(path, plt, bg = bg, width = width, height = height)
    invisible(plt)
}

# ---------------------------------------------------------------------------
# Standalone entry point
# Requires: data/cc_usable_meta_20250805.csv, data/cc_usable_data_20250805.rds,
#           data/predictors_list_20250805.rds
# (produced by prepare_data.R)
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
    library(corrr)
    # in general, run this first
    # source(here("R/create_predictors_response_lists.R"))

    pval_epsilon <- 1e-300
    input_data   <- pull_predictor_response_lists(pull_data = FALSE, date = "20250805")
    cc_usable_meta <- read_csv(here("data/cc_usable_meta_20250805.csv"))
    cc_usable_data <- readRDS(here("data/cc_usable_data_20250805.rds"))

    # ----- predictor correlations -------------------------------------------
    corr_mat_and_vis <- function(df, plot_title = "", row_order = NULL,
                                 method = "pearson") {
        data_matrix <- df %>% select(-target_symbol) %>% as.matrix()
        cor_matrix  <- cor(data_matrix, method = method)
        if (!is.null(row_order)) {
            cor_matrix  <- cor_matrix[row_order, row_order]
            cluster_rows <- FALSE; cluster_cols <- FALSE
        } else {
            cluster_rows <- TRUE;  cluster_cols <- TRUE
        }
        pheatmap_res <- pheatmap::pheatmap(
            cor_matrix, clustering_method = "ward.D2",
            cluster_rows = cluster_rows, cluster_cols = cluster_cols,
            silent = TRUE
        )
        list(mat = cor_matrix, plot = pheatmap_res)
    }

    get_long_corr <- function(corr_mat, data_desc) {
        corr_mat[lower.tri(corr_mat)] <- NA
        corr_mat %>%
            as_tibble(rownames = "var1") %>%
            pivot_longer(-var1, names_to = "var2", values_to = "corr") %>%
            filter(!is.na(corr)) %>%
            mutate(data_transform = data_desc)
    }

    init_lrb    <- corr_mat_and_vis(input_data$predictors_list$log_rank_normalized)
    lrb_order   <- rownames(init_lrb$mat)[init_lrb$plot$tree_row$order]
    cor_results <- list(
        enrichment = corr_mat_and_vis(input_data$predictors_list$raw_enrichment, row_order = lrb_order),
        pval       = corr_mat_and_vis(input_data$predictors_list$raw_pvalue,     row_order = lrb_order),
        lrb        = corr_mat_and_vis(input_data$predictors_list$log_rank_normalized, row_order = lrb_order)
    )
    long_corr_pearson_res <- map(names(cor_results), ~ get_long_corr(cor_results[[.]]$mat, .)) %>%
        bind_rows()

    # ----- replicate correlations -------------------------------------------
    spread_tied_ranks <- function(adjusted_rank) {
        N <- length(adjusted_rank)
        ranked_min <- rank(adjusted_rank, ties.method = "min")
        ranked_avg <- rank(adjusted_rank, ties.method = "average")
        spread_rank <- numeric(N)
        for (val in sort(unique(adjusted_rank))) {
            tie_idx <- which(adjusted_rank == val)
            n_ties  <- length(tie_idx)
            if (n_ties == 1) {
                spread_rank[tie_idx] <- ranked_min[tie_idx[1]]
            } else {
                spread_rank[tie_idx] <- ranked_min[tie_idx[1]] +
                    (ranked_avg[tie_idx[1]] / N) * n_ties
            }
        }
        spread_rank
    }

    cc_ranked <- cc_usable_data %>%
        filter(target_symbol %in% unique(input_data$predictors_list$log_rank$target_symbol)) %>%
        select(regulator_symbol, record_id, target_symbol, callingcards_enrichment, poisson_pval) %>%
        group_by(record_id) %>%
        mutate(primary_rank = rank(poisson_pval, ties.method = "min")) %>%
        group_by(record_id, primary_rank) %>%
        mutate(adjusted_rank = if (n() == 1) primary_rank else
            primary_rank + rank(-callingcards_enrichment, ties.method = "average") / 1e6) %>%
        ungroup() %>%
        group_by(record_id) %>%
        mutate(
            spread_rank = spread_tied_ranks(adjusted_rank),
            lrb         = -log10(spread_rank) + log10(max(spread_rank, na.rm = TRUE)),
            mlog10pval  = -log10(poisson_pval + pval_epsilon)
        ) %>%
        ungroup()

    compute_replicate_correlations <- function(value_col, df) {
        df %>%
            select(regulator_symbol, record_id, target_symbol, !!sym(value_col)) %>%
            group_by(regulator_symbol) %>%
            filter(n_distinct(record_id) > 1) %>%
            ungroup() %>%
            group_by(regulator_symbol) %>%
            nest() %>%
            mutate(
                wide  = map(data, ~ pivot_wider(.x, names_from = record_id,
                                               values_from = !!sym(value_col))),
                corrs = map(wide, ~ {
                    dat     <- select(.x, -target_symbol)
                    pearson  <- correlate(dat, method = "pearson",
                                         use = "pairwise.complete.obs") %>%
                        shave() %>% corrr::stretch(na.rm = TRUE) %>% rename(pearson = r)
                    spearman <- correlate(dat, method = "spearman",
                                         use = "pairwise.complete.obs") %>%
                        shave() %>% corrr::stretch(na.rm = TRUE) %>% rename(spearman = r)
                    full_join(pearson, spearman, by = c("x", "y"))
                })
            ) %>%
            select(regulator_symbol, corrs) %>%
            unnest(corrs) %>%
            rename(record_id_1 = x, record_id_2 = y) %>%
            pivot_longer(-c(regulator_symbol, record_id_1, record_id_2),
                         names_to = "method", values_to = "corr") %>%
            mutate(data_transform = value_col)
    }

    replicate_corr_comp_df <- map(
        c("callingcards_enrichment", "mlog10pval", "lrb"),
        compute_replicate_correlations, cc_ranked
    ) %>%
        bind_rows() %>%
        ungroup() %>%
        mutate(data_transform = factor(
            data_transform,
            levels = c("callingcards_enrichment", "mlog10pval", "lrb"),
            labels = c("Abs(enrichment)", "-log10(pvalue)", "LRB")
        ))

    out_path <- here("plots/fig1b1.svg")
    save_fig1b1(out_path, long_corr_pearson_res, replicate_corr_comp_df)
    message("Saved: ", out_path)
}
