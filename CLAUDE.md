# CLAUDE.md — tfbpmodeling_paper_figures

## What this repo is

Standalone R project containing all publication figure scripts for the
**tfbpmodeling** yeast TF regulatory network paper. This is a refactored
version of
`/home/chase/projects/glmnetstrap/tmp/workflow_result_analysis_20250805.Rmd`. 
It is important that as we refactor that we check for accuracy of this repos
results and process against that notebook.

**Reference repo:** `~/projects/glmnetstrap` — the glmnetstrap R package plus
the full analysis workspace in `tmp/`. Consult it for scientific context,
original notebook code (`tmp/workflow_result_analysis_20250805.Rmd`), and the
source of truth for any logic not yet fully migrated here.

---

## Repo structure

```
R/              Figure scripts + data pipeline (one file per figure panel)
data/           Serialized data objects (gitignored; produced by prepare_data.R)
```

---

## Figure script contract

Each `R/fig*.R` follows this pattern:

- **`make_<fig>(...)`** — takes pre-loaded R data objects as arguments, returns
  a ggplot / ggdraw / patchwork object. Importable into a notebook.
- **`save_<fig>(path, ...)`** — calls `make_<fig>()` + `ggsave()`.
- **Standalone block** (`if (sys.nframe() == 0L)`) — loads from `data/` (pure R,
  no Python required) and writes an SVG to `R/` or a future `figures/` output dir.

Scripts that still require Python (via reticulate): **`prepare_data.R`** only.
All figure scripts themselves should be pure R after the data pipeline runs.

| Script | Figure | Key inputs |
|--------|--------|-----------|
| `fig1ab.R` | TF binding-set size bar charts (FDR vs DTO) | `tf_target_counts` |
| `fig1b1.R` | Replicate vs TF-TF Pearson correlation boxplot | `long_corr_pearson_res`, `replicate_corr_comp_df` |
| `fig1cd.R` | R² linear vs cubic, 4 data transforms | `r2_summary_df` |
| `fig1efgh.R` | Scatter panels linear+cubic fits for a pTF | `uni_modeling_results`, `regulator` |
| `fig2bc.R` | Bootstrap coef boxplots, Stage 1 vs Stage 2 | `all_data_coefs_list`, `topn_coefs_list`, `ci_df_all_data`, `ci_df_topn`, `regulator` |
| `fig2d.R` | Interactor counts per stage (lineplot/barplot/boxplot) | `stage_comp_df` |
| `fig2e.R` | CV variance explained across stages | `stages_comp_r2_df` |
| `figS2A.R` | PDR3×GAL4 cubic + interaction scatter | `input_data`, `pTF`, `mTF` |

---

## Data pipeline

`R/prepare_data.R` is the one-time extraction script. It requires the
`~/htcf_ref/` mount and (for some sections) the `.venv` Python virtualenv.
Run it once per data pull date to populate `data/`:

```r
# Default: date = "20250805", variant = "residuals"
Rscript R/prepare_data.R

# Or source with custom params:
date <- "20250901"
source(here("R/prepare_data.R"))
```

| File | Requires Python |
|------|----------------|
| `data/predictors_list_{date}.rds` | No |
| `data/responses_list_{date}.rds` | No |
| `data/ci_df_all_data_{date}.rds` | No |
| `data/ci_df_topn_{date}.rds` | No |
| `data/stage4_results_{variant}_{date}.rds` | No |
| `data/bootstrap_coefs_all_data_{date}.rds` | **Yes** (pkl) |
| `data/bootstrap_coefs_topn_{date}.rds` | **Yes** (pkl) |
| `data/rr_meta_{date}.csv` | **Yes** (tfbpapi) |
| `data/brent_nf_cc_meta_{date}.csv` | No (file copy) |
| `data/mcisaac_oe_meta_{date}.csv` | No (file copy) |
| `data/cc_usable_meta_{date}.csv` | **Yes** (tfbpapi) |
| `data/cc_usable_data_{date}.rds` | **Yes** (tfbpapi) |
| `data/red_median_wide_{date}.csv` | No (file copy) |

---

## Known issues to resolve (migration from glmnetstrap)

These path references still point to the old glmnetstrap workspace and must be
updated as part of the ongoing refactor:

1. **`source(here("tmp/theme_ptf.R"))`** — `theme_ptf.R` lives at
   `~/projects/glmnetstrap/tmp/theme_ptf.R`. Needs to be copied into `R/` and
   source paths updated across all figure scripts.

2. **`source(here("tmp/fit_ols_model.R"))`** — used by `fig2e.R` and `figS2A.R`.
   Lives at `~/projects/glmnetstrap/tmp/fit_ols_model.R`. Same treatment.

3. **`devtools::load_all(here())`** — in some standalone blocks, calls
   `devtools::load_all()` on the *glmnetstrap* package (for `create_input_df`
   etc.). This repo is not a package. These calls should be replaced with direct
   `source()` of the relevant helpers once they are copied in.

4. **`source(here("R/create_predictors_response_lists.R"))`** —
   the canonical version of this file is now
   `~/projects/glmnetstrap/R/create_predictors_response_lists.R`.
   It should be copied to `R/` and source paths updated.

5. **Output paths** — standalone blocks write SVGs to
   `here("R/fig*.svg")`. Should be updated to
   `here("figures/fig*.svg")` or similar once a `figures/` output dir is agreed on.
