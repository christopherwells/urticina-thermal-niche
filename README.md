# Behavioral metabolic suppression and the thermal niche of an undescribed cold-water *Urticina* sea anemone

Reproducible analysis and figures for a study of an undescribed cold-water intertidal *Urticina* sea anemone in the Northwest Atlantic (Massachusetts through Newfoundland). The paper asks whether a sedentary marine ectotherm's ability to close and retract under thermal stress biases the thermal limits inferred from standard respirometry. Three analyses feed that question: a negative binomial encounter-rate regression and a maximum-entropy (MaxEnt) species distribution model, both fit on an H3 hexagonal grid built from effort-corrected iNaturalist observations, and a set of Bayesian thermal performance curves fit to closed-chamber respirometry across seven temperature treatments (1-30 °C, 18 individuals). The distributional models locate the species' realized thermal envelope; the performance curves show that behavioral closure, rather than physiological failure, drives the apparent high-temperature decline in metabolic rate, so nominal CT~max~ estimates overshoot the observed lethal bracket.

This repository is the trimmed, archival version of the full analysis. It starts from analysis-ready data and contains only the three focal analyses (encounter-rate regression, SDM, and thermal performance curves) plus the code that draws the manuscript figures. The upstream steps that download and clean raw environmental rasters and the raw iNaturalist export are not included, because those inputs are large and not intended for redistribution.

## Contents

| File | What it is |
|------|------------|
| `analysis.qmd` | The whole analysis in one document: the negative binomial encounter-rate regression, the MaxEnt SDM and the concordance between them, and the Bayesian nonlinear mixed-effects thermal performance curves (four TPC shapes crossed with three behavioral structures, model comparison by LOO-IC, and the behavioral decomposition). Reads `data/clean/`, writes the fitted distributional models and predictions back to `data/clean/`. The `brms` sampling is kept as documentation but not re-run on render; the thermal-performance section loads the shipped fits instead. Render with `quarto render analysis.qmd`. |
| `make_figures.R` | Draws manuscript Figures 1-5 and writes them to `figures/`. Run with `Rscript make_figures.R`. |
| `R/utils.R` | Shared configuration, I/O helpers (`load_output()` / `save_output()`), map themes, and the hex-suitability extractor. Sourced by both `analysis.qmd` and `make_figures.R`. |
| `data/clean/` | Analysis-ready inputs (see [Data files](#data-files)). |

Render `analysis.qmd` first, then run `make_figures.R`; see [Reproducing the analysis](#reproducing-the-analysis).

## Data provenance

The analysis-ready data in `data/clean/` are derived products.

- **Species occurrences and sampling effort** come from [iNaturalist](https://www.inaturalist.org), whose observations are contributed by the public and released under Creative Commons licenses. The 357 research-grade *Urticina* presence points, the observer-effort tallies, and the catch-per-unit-effort (CPUE) surface on the H3 grid are all built from that public record.
- **Environmental predictors** are hex-level zonal summaries of published climate and ocean layers: sea surface temperature, salinity, sea ice, chlorophyll, pH, and bottom temperature from **Bio-ORACLE v3.0**; cloud cover and wind speed from **ERA5** reanalysis; air temperature, precipitation, and vapor pressure deficit from **CHELSA v2.1**; and solar radiation from **WorldClim 2.1**.
- **Respirometry** measurements (`anemone_oxygen.xlsx`, `anemone_weights.xlsx`) are the authors' own laboratory data.

The raw environmental rasters are **not redistributed here** — download them from the sources above if you need to rebuild the predictor tables. What this repository ships is the small set of extracted, collinearity-screened tables that the models actually consume.

Code is released under the MIT License; the data files under CC-BY-4.0. See `LICENSE`.

## Reproducing the analysis

Built and tested with **R 4.6.0**. Rendering does not re-sample the Bayesian TPCs, so a CmdStan toolchain is not needed to reproduce the analysis. Regenerating the fits from scratch (rather than using the shipped `data/clean/06-tpc-fits.rds`) does require a working [CmdStan](https://mc-stan.org/cmdstanr/) toolchain via `cmdstanr`.

Key package versions:

| Package | Version | | Package | Version |
|---------|---------|---|---------|---------|
| sf | 1.1-0 | | brms | 2.23.0 |
| terra | 1.9-11 | | cmdstanr | 0.9.0 |
| h3jsr | 1.3.1 | | loo | 2.9.0 |
| glmmTMB | 1.1.14 | | ggplot2 | 4.0.3 |
| maxnet | 0.1.4 | | cowplot | 1.2.0 |
| ENMeval | 2.0.5.2 | | patchwork | 1.3.2 |
| blockCV | 3.2-0 | | viridis | 0.6.5 |
| DHARMa | 0.4.7 | | kableExtra | 1.4.0 |
| spdep | 1.4-2 | | tidyterra | 1.1.0 |
| mgcv | 1.9-4 | | rnaturalearthhires | 1.0.0.9000 |

`rnaturalearthhires` is on GitHub rather than CRAN: `remotes::install_github("ropensci/rnaturalearthhires")`. The figure script also uses `ragg` and `png` for text measurement and the optional organism inset.

Run the pieces in this order:

```sh
quarto render analysis.qmd   # fits NB + MaxEnt, loads the shipped TPC fits, writes 04-*/05-* into data/clean/
Rscript make_figures.R        # writes Figures 1-5 to figures/
```

Notes:

- **`analysis.qmd` must run before `make_figures.R`.** Figures 1-3 and 5 depend on the fitted NB and MaxEnt objects, which `analysis.qmd` writes into `data/clean/` (`04-nb-model.rds`, `04-scale-params.rds`, `04-nb-selected-predictors.rds`, `05-maxent-model.rds`, `05-sdm-prediction.rds`). Figure 4 reads the fitted TPC objects in `data/clean/06-tpc-fits.rds`, which are shipped, so it does not require re-sampling.
- **The thermal performance curves are not re-sampled on render.** The fitted `brms` objects are shipped in `data/clean/06-tpc-fits.rds`, and the thermal-performance section of `analysis.qmd` loads them; every downstream quantity (model comparison, T~opt~, CT~max~, the behavioral multiplier) is derived from the loaded fits. The full sampling code (fifteen `brms` models across four TPC shapes and three behavioral structures) is kept in the `brms-fit` chunk as documentation but is marked `eval: false`, which is what keeps the render fast. To regenerate the fits from scratch, run that chunk; a per-model cache is written to `data/clean/06-brms-cache/` and reused thereafter.
- **Figure 1's organism inset** is a photograph that is not archived here. Without it, `make_figures.R` writes the study-area map without the inset. To include the inset, place the photo at `figures/fig1-organism-inset.png` before running the script.

## Data files

| File | Description |
|------|-------------|
| `01-urticina-obs.rds` | 357 research-grade iNaturalist *Urticina* presence points (`sf`, POINT). |
| `01-study-extent.rds` | Study-area extent polygon (`sf`). |
| `02-cpue-hex.rds` | H3 resolution-5 coastal hex grid with observer effort, *Urticina* counts, and CPUE (`sf`). |
| `02-h3-grid.rds` | H3 resolution-5 coastal hex grid: geometry and hex addresses (`sf`). |
| `03-predictors-uncorr.rds` | Hex grid with the collinearity-screened environmental predictors used by the models (`sf`). |
| `03-predictors-hex.rds` | Hex grid with the full environmental predictor set before collinearity screening (`sf`). |
| `03-sst-raster.rds` | Summer/winter sea surface temperature layers used for the TPC-versus-SST overlay (`SpatRaster`, packed). |
| `06-tpc-fits.rds` | Fitted `brms` TPC objects, the best-model selection, posterior-derived parameters, and the behavioral decomposition (list). |
| `anemone_oxygen.xlsx` | Closed-chamber respirometry: experimental and control (blank) chamber measurements. |
| `anemone_weights.xlsx` | Anemone dry-mass measurements used to compute mass-specific rates. |

Rendering `analysis.qmd` adds the derived model outputs (`04-nb-model.rds`, `04-nb-summary.rds`, `04-scale-params.rds`, `04-nb-selected-predictors.rds`, `05-maxent-model.rds`, `05-sdm-eval.rds`, `05-sdm-prediction.rds`) to `data/clean/`; these are regenerated from the shipped inputs and are not part of the archived data set.

## Archive

Archived at Zenodo: https://doi.org/10.5281/zenodo.XXXXXXX

## Cite

Edgar, Penfold, Martinez, and Wells, 2026. Behavioral metabolic suppression confounds thermal performance estimates and climate vulnerability assessments in a marine ectotherm. Data and code archived at Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX
