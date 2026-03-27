# tfbpmodeling_paper_figures

Publication figures for the tfbpmodeling yeast TF regulatory network paper.  

## Structure

```
R/              Figure scripts and shared helpers
data/           Serialized data objects (gitignored)
plots/          Output SVGs (gitignored)
```

The `data/` directory stores the results from tfbpmodeling in `data/vx.x.x`, 
eg `data/v1.1.0`. The processed results are stored under the data freeze 
date and the processed results are stored under the datafreeze date and
the tfbpmodeling version, eg `data/20250805` and `data/20250805/1.0.0`. See
`R/prepare_input_data.R` (note: this will not run because the database will not 
be up. However, all of the data that this produces is already stored in `data`, 
so there is no need) and `R/prepare_results_data.R`.

## Notebooks and figure scripts

The notebook `paper_figures.Rmd` and `chisqr.Rmd` both produce figures that
are used in the paper. With the virtual environment installed via
`reticulate` (use `requirements.txt`), you can run each to produce all of the
figures in `plots/`. You can also run the figure creation scripts individually.
Scroll down to the bottom of each file in `R/` for code that will run in the
console to individually create a figure.
