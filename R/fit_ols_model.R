# Construct input data for a specific regulator
create_input_df <- function(regulator,
                            response_df,
                            predictors,
                            feature_col = "target_symbol") {
    expected_cols <- c(feature_col, regulator)

    response_df = response_df %>%
        select(all_of(c(feature_col, regulator)))

    if (!setequal(colnames(response_df), expected_cols)) {
        stop(sprintf("`response_df` must have exactly two columns: '%s' and '%s'", feature_col, regulator), call. = FALSE)
    }

    stopifnot(all(c(feature_col, regulator) %in% colnames(predictors)))

    out_df <- response_df %>%
        dplyr::rename(response = !!rlang::sym(regulator)) %>%
        dplyr::left_join(predictors, by = feature_col) %>%
        dplyr::filter(.data[[feature_col]] != regulator) %>%
        as.data.frame()

    row_ids <- out_df[[feature_col]]

    out_df = out_df %>%
        dplyr::select(-dplyr::all_of(feature_col))

    rownames(out_df) <- row_ids

    # Count number of rows before and after filtering
    n_before <- nrow(out_df)

    # Remove rows with any NA/NaN/Inf
    out_df <- out_df %>%
        dplyr::filter(if_all(everything(), ~ is.finite(.)))

    n_after <- nrow(out_df)

    if (n_before > n_after) {
        message(glue::glue(
            "Dropped {n_before - n_after} rows with NA/NaN/Inf for regulator '{regulator}'"
        ))
    }

    out_df
}

# Create interaction model formula
create_interaction_formula <- function(ptf, mtf_list, response_var = "response", no_const_term = TRUE, p = 1) {
    # Main effect term
    linear_term <- if (p >= 1) ptf else NULL

    # Polynomial terms beyond degree 1, wrapped in I()
    poly_terms <- if (p > 1) paste0("I(", ptf, "^", 2:p, ")") else NULL

    # Interaction terms
    interactors <- paste(ptf, mtf_list, sep = ":")

    # Constant term
    base <- if (no_const_term) "-1" else "1"

    # Assemble all terms
    all_terms <- c(base, linear_term, poly_terms, interactors)

    # Create and return formula
    as.formula(paste(response_var, "~", paste(all_terms, collapse = " + ")))
}


# example
# out = fit_ols_model("STB5", responses_list$rank, predictors_list$raw_enrichment, cubic_model=TRUE)

fit_ols_model <- function(regulator, responses, predictors, cubic_model = FALSE, top600 = FALSE) {

    if (top600) {
        selected_target_symbols <- predictors_list$rank %>%
            dplyr::select(target_symbol, !!rlang::sym(regulator)) %>%
            filter(target_symbol != regulator) %>%
            arrange(desc(!!rlang::sym(regulator))) %>%
            slice_head(n = 600) %>%
            pull(target_symbol)
    } else {
        selected_target_symbols <- predictors$target_symbol
    }

    response_df <- responses %>%
        dplyr::select(target_symbol, !!rlang::sym(regulator)) %>%
        filter(target_symbol %in% selected_target_symbols)

    predictors_df <- predictors %>%
        filter(target_symbol %in% selected_target_symbols)

    x <- create_input_df(regulator, response_df, predictors_df)

    if (!all(c("response", regulator) %in% colnames(x))) {
        stop(glue::glue("Missing required columns for {regulator}"))
    }

    if (cubic_model) {
        formula <- as.formula(paste("response ~ poly(", regulator, ", degree = 3)", sep = ""))
    } else {
        formula <- as.formula(glue::glue("response ~ {regulator}"))
    }

    message(glue::glue("Fitting using formula: {paste(deparse(formula), collapse = '')}"))
    fit <- lm(formula, data = x)

    fit_df = broom::augment(fit) %>%
        mutate(target_symbol = rownames(x))

    leverage_vector = fit_df %>%
        select(target_symbol, .hat) %>%
        arrange(desc(.hat)) %>%
        deframe()

    plot <- ggplot(x, aes(x = !!rlang::sym(regulator), y = response)) +
        geom_point(alpha = 0.5, size = 1.3)

    if (cubic_model) {
        # Fit raw cubic to fitted values for derivative analysis
        x_vals <- x[[regulator]]
        fitted_vals <- predict(fit)
        raw_fit <- lm(fitted_vals ~ poly(x_vals, degree = 3, raw = TRUE))
        coef_raw <- coef(raw_fit)

        first_derivative <- function(x) {
            coef_raw[2] + 2 * coef_raw[3] * x + 3 * coef_raw[4] * x^2
        }

        second_derivative <- function(x) {
            2 * coef_raw[3] + 6 * coef_raw[4] * x
        }

        a <- 3 * coef_raw[4]
        b <- 2 * coef_raw[3]
        c <- coef_raw[2]
        crit_points <- polyroot(c(c, b, a))
        crit_points_real <- Re(crit_points[Im(crit_points) == 0])
        x_range <- range(x[[regulator]], na.rm = TRUE)
        crit_points_real <- crit_points_real[crit_points_real >= x_range[1] & crit_points_real <= x_range[2]]

        plot <- plot +
            geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE) +
            labs(title = "Cubic Polynomial Fit")

    } else {
        first_derivative <- NULL
        second_derivative <- NULL
        crit_points_real <- NULL

        plot <- plot +
            geom_smooth(method = "lm", se = FALSE) +
            labs(title = "Linear Fit")
    }

    plot <- plot +
        annotate(
            "text",
            x = Inf, y = Inf,
            label = sprintf("R^2 = %.3f", summary(fit)$r.squared),
            hjust = 1.1, vjust = 1.5,
            size = 4,
            fontface = "italic"
        ) +
        labs(
            x = paste0(regulator, " Binding"),
            y = paste0(regulator, " Perturbation Response")
        ) +
        theme_minimal(base_size = 14) +
        theme(
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12)
        )

    list(
        fit = fit,
        first_derivative = first_derivative,
        second_derivative = second_derivative,
        critical_points = crit_points_real,
        plot = plot,
        leverage = leverage_vector
    )
}
