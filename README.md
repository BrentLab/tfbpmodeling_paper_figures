# tfbpmodeling_paper_figures

Publication figures for the tfbpmodeling yeast TF regulatory network paper.
Refactored from `~/projects/glmnetstrap/tmp/workflow_result_analysis_20250805.Rmd`.

## Structure

```
R/              Figure scripts and shared helpers
data/           Serialized data objects (gitignored)
plots/          Output SVGs (gitignored)
```

## Data pipeline

Run once to populate `data/` from the htcf mount. Requires Python/reticulate
for sections 4–6.

```r
date      <- "20250805"
variant   <- "residuals"
pull_data <- TRUE
source(here("R/prepare_data.R"))
```

Set `pull_data <- FALSE` to reload from already-serialized RDS files without
hitting the htcf mount.

## Figure scripts

Each `R/fig*.R` exports `make_<fig>(...)` returning a ggplot/patchwork object
and `save_<fig>(path, ...)` which calls `ggsave()`. Source the script and pass
pre-built data objects:

```r
source(here("R/fig1ab.R"))
plt <- make_fig1ab(tf_target_counts)
```

Shared helpers sourced by the figure scripts:

| File | Provides |
|------|---------|
| `R/theme_ptf.R` | `theme_ptf()` ggplot theme |
| `R/fit_ols_model.R` | `create_input_df()`, `fit_ols_model()` |
| `R/create_predictors_response_lists.R` | `pull_predictor_response_lists()` |

## Figure index

| Script | Figure |
|--------|--------|
| `fig1ab.R` | TF binding-set size distributions (FDR vs DTO) |
| `fig1b1.R` | Replicate vs TF-TF Pearson correlation |
| `fig1cd.R` | R² linear vs cubic across four data transforms |
| `fig1efgh.R` | Scatter panels with linear and cubic fits for a pTF |
| `fig2bc.R` | Bootstrap coefficient boxplots, Stage 1 vs Stage 2 |
| `fig2d.R` | Interactor counts per stage |
| `fig2e.R` | CV variance explained across pipeline stages |
| `figS2A.R` | PDR3×GAL4 cubic fit and interaction scatter |
