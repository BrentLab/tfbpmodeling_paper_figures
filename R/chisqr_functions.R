library(tidyverse)
library(broom)
library(glue)
library(here)
library(patchwork)

## binning function

classify_mTF_bins <- function(df, mTF, quantile_cutoff, cap_n = NULL) {
    bottom_cutoff <- quantile(df[[mTF]], probs = quantile_cutoff, na.rm = TRUE)
    top_cutoff <- quantile(df[[mTF]], probs = 1 - quantile_cutoff, na.rm = TRUE)

    if (bottom_cutoff > 0.5 || top_cutoff > 0.5) {
        warning(glue::glue("Bottom/top cutoff > 0.5: bottom={bottom_cutoff}, top={top_cutoff}"))
    }

    df <- df %>%
        mutate(mTF_bin = case_when(
            .data[[mTF]] <= bottom_cutoff ~ "Bottom",
            .data[[mTF]] >= top_cutoff ~ "Top",
            TRUE ~ "Middle"
        ))

    if (!is.null(cap_n)) {
        bottom_df <- df %>%
            filter(mTF_bin == "Bottom") %>%
            arrange(.data[[mTF]]) %>%
            slice_head(n = cap_n)

        top_df <- df %>%
            filter(mTF_bin == "Top") %>%
            arrange(desc(.data[[mTF]])) %>%
            slice_head(n = cap_n)

        middle_df <- df %>%
            filter(mTF_bin == "Middle")

        df <- bind_rows(bottom_df, top_df, middle_df)
    }

    return(df)
}

## Classify pTF

classify_pTF_bin <- function(df, pTF, threshold) {
    df %>%
        mutate(pTF_bin = cut(.data[[pTF]],
                             breaks = c(-Inf, threshold, Inf),
                             labels = c("Low", "High")))
}

## linear model

compute_linear_slope <- function(
        predictors,
        response,
        subset_index,
        pTF
) {
    x <- predictors %>%
        filter(target_symbol %in% subset_index) %>%
        select(!!rlang::sym(pTF)) %>%
        pull()

    y <- response %>%
        filter(target_symbol %in% subset_index) %>%
        select(!!rlang::sym(pTF)) %>%
        rename(response = !!rlang::sym(pTF)) %>%
        pull()
    model <- lm(y ~ x)
    coef(model)[2]  # slope
}

## compare CI

is_strictly_positive <- function(ci) ci[1] > 0 && ci[2] > 0
is_strictly_negative <- function(ci) ci[1] < 0 && ci[2] < 0
crosses_zero <- function(ci) ci[1] * ci[2] < 0

compare_confidence_intervals <- function(main_ci, interaction_ci) {
    if (is.null(main_ci)) return("NA")
    if (is.null(interaction_ci)) return("NA")
    if (crosses_zero(main_ci) || crosses_zero(interaction_ci)) return("NA")

    if (is_strictly_positive(main_ci) && is_strictly_positive(interaction_ci)) {
        return(if (main_ci[1] > interaction_ci[1]) "main effect greater" else "main effect less")
    }
    if (is_strictly_negative(main_ci) && is_strictly_negative(interaction_ci)) {
        return(if (main_ci[2] < interaction_ci[2]) "main effect greater" else "main effect less")
    }

    main_mag <- if (is_strictly_negative(main_ci)) abs(main_ci[1]) else main_ci[1]
    inter_mag <- if (is_strictly_negative(interaction_ci)) abs(interaction_ci[1]) else interaction_ci[1]

    if (main_mag > inter_mag) "main effect greater" else "main effect less"
}

## main function

get_interaction_term_chi_squared <- function(
        significant_predictors_list,
        predictors_df,
        response_df,
        threshold = 0.6,
        quantile_cutoff = 0.10) {
    stopifnot(quantile_cutoff <= 0.5)

    # note that the final output shape changed
    interaction_terms <- map(significant_predictors_list, ~.$interactor) # names(significant_predictors_list)[str_detect(names(significant_predictors_list), ":")]

    if(length(interaction_terms) == 0){
        return(NA)
    }

    pTF = str_split(interaction_terms[[1]], ":")[[1]][[1]]

    predictors_df_copy <- predictors_df %>%
        filter(target_symbol != pTF)

    response_df_copy <- response_df %>%
        filter(target_symbol != pTF)

    results <- list()
    for (term in interaction_terms) {
        parts <- str_split(term, ":", simplify = TRUE)
        pTF <- parts[1]
        mTF <- parts[2]

        df <- predictors_df_copy %>%
            select(target_symbol, all_of(c(pTF, mTF))) %>%
            drop_na()

        df <- classify_mTF_bins(df, mTF, quantile_cutoff)
        df <- classify_pTF_bin(df, pTF, threshold)

        low_index <- df %>%
            filter(tolower(mTF_bin) == "bottom") %>%
            pull(target_symbol)

        high_index <- df %>%
            filter(tolower(mTF_bin) == "top") %>%
            pull(target_symbol)

        low_slope <- compute_linear_slope(predictors_df_copy, response_df_copy, low_index, pTF)
        high_slope <- compute_linear_slope(predictors_df_copy, response_df_copy, high_index, pTF)

        df_filtered <- df %>%
            filter(tolower(mTF_bin) %in% c("bottom", "top"))

        contingency <- table(df_filtered$pTF_bin, df_filtered$mTF_bin)
        chi2 <- suppressWarnings(chisq.test(contingency, correct = FALSE))

        expected <- as.data.frame.matrix(chi2$expected)
        # observed_gt_expected <- tryCatch({
        #     contingency["High", "Bottom"] > expected["High", "Bottom"]
        # }, error = function(e) message(glue::glue("error calculating observed_gt_expected: {e}")))

        # first, figure out if "Bottom" mtf or "Top" mtf has the greater slope
        selected_side_mtf = if (low_slope > high_slope) "Bottom" else "Top"
        # Next, check whether High is greater than expected
        unbalanced_high_ptf = tryCatch({
            contingency["High", selected_side_mtf] > expected["High", selected_side_mtf]
        }, error = function(e) message(glue::glue("error calculating __: {e}")))

        # NOTE: when using residuals as response, and if the pTF is not in the
        # original model matrix, then it will always be NA
        compare_result <- compare_confidence_intervals(
            significant_predictors_list[[pTF]],
            significant_predictors_list[[term]])

        # Low_mTF_High_pTF_observed_gt_expected = observed_gt_expected,

        results[[term]] <- list(
            chi2_stat = chi2$statistic,
            pvalue = chi2$p.value,
            dof = chi2$parameter,
            expected_freq = expected,
            raw_table = as.data.frame.matrix(contingency),
            unbalanced_high_ptf = unbalanced_high_ptf,
            main_effect_compare = compare_result,
            ptf_by_low_mtf_slope = low_slope,
            ptf_by_high_mtf_slope = high_slope,
            mtf_low_index = low_index,
            mtf_high_index = high_index
        )
    }

    results
}

get_interaction_term_chi_squared_stage3_bootstrap <- function(
        ci_df,
        predictors_df,
        response_df,
        threshold = 0.6,
        quantile_cutoff = 0.10) {
    stopifnot(quantile_cutoff <= 0.5)

    # ci_df is the subset of ci_df_stage3_bootstrap for one regulator
    # interaction terms are rows where interactor contains ":"
    interaction_terms <- ci_df %>%
        filter(str_detect(interactor, ":")) %>%
        pull(interactor)

    if (length(interaction_terms) == 0) {
        return(NA)
    }

    pTF <- str_split(interaction_terms[[1]], ":")[[1]][[1]]

    predictors_df_copy <- predictors_df %>%
        filter(target_symbol != pTF)

    response_df_copy <- response_df %>%
        filter(target_symbol != pTF)

    # Look up [ci_lo, ci_hi] for a given interactor name; NULL if absent
    # (bootstrap may have dropped the main effect entirely)
    get_ci <- function(interactor_name) {
        row <- ci_df %>% filter(interactor == interactor_name)
        if (nrow(row) == 0) return(NULL)
        c(row$ci_lo, row$ci_hi)
    }

    results <- list()
    for (term in interaction_terms) {
        parts <- str_split(term, ":", simplify = TRUE)
        pTF <- parts[1]
        mTF <- parts[2]

        df <- predictors_df_copy %>%
            select(target_symbol, all_of(c(pTF, mTF))) %>%
            drop_na()

        df <- classify_mTF_bins(df, mTF, quantile_cutoff)
        df <- classify_pTF_bin(df, pTF, threshold)

        low_index <- df %>%
            filter(tolower(mTF_bin) == "bottom") %>%
            pull(target_symbol)

        high_index <- df %>%
            filter(tolower(mTF_bin) == "top") %>%
            pull(target_symbol)

        low_slope  <- compute_linear_slope(predictors_df_copy, response_df_copy, low_index,  pTF)
        high_slope <- compute_linear_slope(predictors_df_copy, response_df_copy, high_index, pTF)

        df_filtered <- df %>%
            filter(tolower(mTF_bin) %in% c("bottom", "top"))

        contingency <- table(df_filtered$pTF_bin, df_filtered$mTF_bin)
        chi2 <- suppressWarnings(chisq.test(contingency, correct = FALSE))
        expected <- as.data.frame.matrix(chi2$expected)

        selected_side_mtf <- if (low_slope > high_slope) "Bottom" else "Top"
        unbalanced_high_ptf <- tryCatch({
            contingency["High", selected_side_mtf] > expected["High", selected_side_mtf]
        }, error = function(e) message(glue::glue("error calculating unbalanced_high_ptf: {e}")))

        # pTF main effect CI may be NULL if bootstrap dropped it — handled by compare_confidence_intervals
        compare_result <- compare_confidence_intervals(
            get_ci(pTF),
            get_ci(term)
        )

        results[[term]] <- list(
            chi2_stat             = chi2$statistic,
            pvalue                = chi2$p.value,
            dof                   = chi2$parameter,
            expected_freq         = expected,
            raw_table             = as.data.frame.matrix(contingency),
            unbalanced_high_ptf   = unbalanced_high_ptf,
            main_effect_compare   = compare_result,
            ptf_by_low_mtf_slope  = low_slope,
            ptf_by_high_mtf_slope = high_slope,
            mtf_low_index         = low_index,
            mtf_high_index        = high_index
        )
    }

    results
}

parse_interactor_pair <- function(interactor_pair) {
    tryCatch({
        split_interactor_pair <- str_split(interactor_pair, ":", simplify = TRUE)

        if (ncol(split_interactor_pair) < 2) {
            stop(glue::glue("Interactor pair '{interactor_pair}' does not contain a valid ':' separator."))
        }

        pTF <- split_interactor_pair[1]
        mTF <- split_interactor_pair[2]
        list(pTF = pTF, mTF = mTF)
    }, error = function(e) {
        stop(glue::glue("Failed to parse interactor pair '{interactor_pair}': {e$message}"))
    })
}

extract_indices <- function(results, pTF, interactor_pair) {
    tryCatch({
        res_entry <- results[[pTF]][[interactor_pair]]
        if (is.null(res_entry)) {
            stop(glue::glue("No entry found for pTF = '{pTF}' and interactor = '{interactor_pair}'"))
        }

        if (is.null(res_entry$mtf_low_index) || is.null(res_entry$mtf_high_index)) {
            stop(glue::glue("Missing 'mtf_low_index' or 'mtf_high_index' in results[['{pTF}']][['{interactor_pair}']]"))
        }

        list(
            low = res_entry$mtf_low_index,
            high = res_entry$mtf_high_index
        )
    }, error = function(e) {
        stop(glue::glue("Failed to extract indices for pTF = '{pTF}', interactor = '{interactor_pair}': {e$message}"))
    })
}

plot_linear_fits <- function(
        predictors,
        response,
        results,
        interactor_str) {

    interactor_pair = parse_interactor_pair(interactor_str)

    indicies = extract_indices(results, interactor_pair$pTF, interactor_str)


    # Styling
    red_color <- "#D21E32"
    blue_color <- "#3396C8"
    text_size <- 22

    mtf_styled <- str_to_sentence(interactor_pair$mTF)
    ptf_styled <- str_to_sentence(interactor_pair$pTF)

    # Function to compute regression fit and CI
    compute_fit_df <- function(predictors, response, subset_index) {

        # create a dataframe with two columns (response, pTF) which is
        # subsetted by the subset_index, the list of target_symbols chosen by
        # mTF
        df = predictors %>%
            filter(target_symbol %in% subset_index) %>%
            select(target_symbol, !!rlang::sym(interactor_pair$pTF)) %>%
            dplyr::rename(x = !!rlang::sym(interactor_pair$pTF)) %>%
            left_join(
                response %>%
                    filter(target_symbol %in% subset_index) %>%
                    select(target_symbol, !!rlang::sym(interactor_pair$pTF)) %>%
                    dplyr::rename(y = !!rlang::sym(interactor_pair$pTF))) %>%
            select(-target_symbol)

        lm_fit <- lm(y ~ x, data = df)
        cubic_poly_fit = lm(y~poly(x, 3, raw=TRUE), data=df)

        new_data <- tibble(x = seq(0, 1, length.out = 100))

        pred <- predict(lm_fit, newdata = new_data)
        pred_cubic = predict(cubic_poly_fit, newdata=new_data)

        fit_df <- bind_cols(new_data,tibble(linear_fit = pred), tibble(cubic_fit=pred_cubic))

        mean_mTF = predictors %>%
            filter(target_symbol %in% subset_index) %>%
            pull(!!rlang::sym(interactor_pair$mTF)) %>%
            mean()

        list(
            fit_df = fit_df,
            scatter_df = df,
            mean_mTF = mean_mTF)
    }

    # Generate data for both subsets
    low_data <- compute_fit_df(predictors, response, indicies$low)
    high_data <- compute_fit_df(predictors, response, indicies$high)

    # Generate data for both subsets
    low_label <- paste0("Weakest ", mtf_styled, " binding\n(mean = ", round(low_data$mean_mTF, 2), ")")
    high_label <- paste0("Strongest ", mtf_styled, " binding\n(mean = ", round(high_data$mean_mTF, 2), ")")

    low_data$scatter_df$group <- low_label
    high_data$scatter_df$group <- high_label
    low_data$fit_df$group <- low_label
    high_data$fit_df$group <- high_label


    scatter_df <- bind_rows(low_data$scatter_df, high_data$scatter_df)


    fit_df <- bind_rows(low_data$fit_df, high_data$fit_df)

    # Ensure consistent factor levels in both frames
    group_levels <- c(low_label, high_label)
    scatter_df$group <- factor(scatter_df$group, levels = group_levels)
    fit_df$group <- factor(fit_df$group, levels = group_levels)

    # Plot
    p <- ggplot() +
        geom_point(data = scatter_df, aes(x = x, y = y, color = group), size = 1.5, alpha = 1.0) +
        geom_line(data = fit_df, aes(x = x, y = linear_fit), color = "black", linewidth = 1.1) +
        # geom_line(data = fit_df, aes(x = x, y = cubic_fit), color = "purple", linetype="dashed", linewidth = 1.1) +
        facet_wrap(~group) +
        scale_color_manual(values = c(blue_color, red_color)) +
        labs(
            x = paste0(ptf_styled, " binding (LRB)"),
            y = paste0(ptf_styled, " response (LRR)")
        ) +
        theme_minimal(base_size = text_size) +
        theme(
            strip.text = element_text(size = text_size),
            axis.text = element_text(color = "black", size = text_size * 0.8),
            axis.title = element_text(color = "black", size = text_size),
            legend.position = "none"
        ) +
        coord_cartesian(xlim = c(-0.02, 1.02), ylim = c(-0.02, 1.02)) +
        scale_x_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0))

    return(p)
}
