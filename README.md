# Mechanism-Based Flood Risk Analysis

R implementation of the mechanism-based flood-frequency framework used in the manuscript:

> **A Mechanism-Based Framework for Compound Flood Risk Assessment Using Flood-Generating Mechanism Attribution and Multivariate Frequency Analysis**

The workflow is deliberately ordered so that each stage produces the inputs required by the next:

1. Flood-generating mechanism event-date extraction
2. Circular seasonality statistics
3. Dominant-mechanism identification using unconstrained quadratic optimization
4. Marginal-distribution fitting and copula-based multivariate flood-frequency analysis

The framework considers the following mechanisms:

- `D1_AMX_Pr`: annual maximum precipitation
- `D2_AMX_MD_Pr`: multi-day precipitation (2-day precipitation)
- `D3_CSMP`: compound soil-moisture–precipitation mechanism
- `D4_CWSP`: compound wind-speed–precipitation mechanism

For the soil-moisture mechanism, **antecedent soil moisture is explicitly retained**: the candidate soil-moisture extreme is selected only when high precipitation occurs within the prescribed forward window. This preserves the antecedent catchment-condition interpretation used in the methodology.

## Repository structure

```text
mechanism_based_flood_analysis/
├── README.md
├── LICENSE
├── .gitignore
├── config/
│   └── analysis_config.R
├── R/
│   ├── 00_run_all.R
│   ├── 01_extract_flood_mechanisms.R
│   ├── 02_circular_statistics.R
│   ├── 03_dominant_mechanism_optimization.R
│   ├── 04_multivariate_flood_frequency.R
│   └── utils.R
├── data/
│   ├── input/
│   │   └── README.md
│   └── output/
│       └── .gitkeep
└── docs/
    └── methodology.md
```

## Required R packages

```r
install.packages(c(
  "data.table",
  "lubridate",
  "zoo",
  "hydroTSM",
  "hydromad",
  "lmom",
  "goftest",
  "VineCopula"
))
```

## Input data

The scripts assume that each hydrometeorological variable is supplied as a daily CSV with:

- rows = dates
- columns = catchments/stations
- identical column ordering across variables

Required inputs:

- daily precipitation
- daily soil moisture
- daily wind speed
- annual maximum streamflow/flood series with corresponding dates

See `data/input/README.md`.

## Run

Edit `config/analysis_config.R`, then run:

```r
source("R/00_run_all.R")
```

or execute the four numbered scripts sequentially.

## Important implementation choices

### Marginal distributions

Eleven L-moment candidate distributions are evaluated:

- Gumbel
- GEV
- Log-Pearson III
- Weibull
- Exponential
- Gamma
- Generalized Logistic
- Generalized Normal
- Kappa
- Log-normal III
- Normal

The distribution with the smallest Cramér–von Mises statistic is selected separately for each marginal.

### Copula

The implementation uses the Archimedean copula families required by the Kendall-distribution formulation:

- Clayton
- Gumbel
- Frank

The copula is selected by `VineCopula::BiCopSelect`.

### Return-period scenarios

The multivariate analysis calculates:

- AND
- OR
- Kendall

For D1/D2 the analysis is bivariate (streamflow + precipitation mechanism). For D3 it is trivariate (streamflow + precipitation + antecedent soil moisture), and for D4 it is trivariate (streamflow + precipitation + wind speed).

For a trivariate case, inclusion–exclusion is used for the AND/OR probabilities.

## Reproducibility

The analysis is deterministic apart from any package-level numerical optimization behaviour. A seed is set in `00_run_all.R`.

## Outputs

The main outputs are written to `data/output/`:

- mechanism event dates
- circular-statistics vectors and concentration
- optimization weights and dominant mechanisms
- selected marginal distributions
- selected copula families and parameters
- AND/OR/Kendall return periods
- discharge estimates corresponding to requested return periods
