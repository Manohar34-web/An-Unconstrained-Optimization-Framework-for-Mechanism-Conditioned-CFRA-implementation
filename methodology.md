# Computational methodology

## 1. Flood-generating mechanism identification

Four candidate flood-generating mechanisms are represented:

1. annual maximum precipitation (`D1_AMX_Pr`);
2. annual maximum multi-day precipitation (`D2_AMX_MD_Pr`);
3. compound soil-moisture–precipitation (`D3_CSMP`);
4. compound wind-speed–precipitation (`D4_CWSP`).

February 29 is removed so that annual blocks contain 365 days. The original
analysis period and initial-record trimming are retained through the
configuration file.

### D3: antecedent soil-moisture condition

The soil-moisture mechanism is not treated as a simple soil-moisture maximum.
The algorithm ranks high-soil-moisture days and searches forward for a
precipitation event exceeding the annual positive-precipitation percentile
threshold. The first qualifying precipitation concurrence identifies the
representative compound event.

This preserves the physical interpretation of soil moisture as an
**antecedent catchment condition** rather than an independent flood trigger.

### D4: wind–precipitation concurrence

Wind-speed extremes are ranked from highest to lowest. A wind extreme is
retained only when precipitation exceeds the specified threshold within the
defined concurrence window. This follows the requirement that high wind speed
alone is insufficient to constitute the flood-generating mechanism.

## 2. Circular statistics

For each mechanism and annual maximum flood series, the day of occurrence is
transformed into an angle:

`theta = 2*pi*D/m`

where `D` is the day of occurrence and `m = 365`.

The mean cosine and sine components are calculated as:

`xbar = mean(cos(theta))`

`ybar = mean(sin(theta))`

and the concentration parameter is:

`R = sqrt(xbar^2 + ybar^2)`.

The circular mean date is recovered from `atan2(ybar, xbar)`.

## 3. Dominant mechanism

The flood seasonality vector is represented by:

`v_f = (x_f, y_f)'`

and each mechanism by:

`v_i = (x_i, y_i)'`.

The mechanism vectors form:

`V = [v_1, ..., v_n]`.

The relative importance vector is obtained from the unconstrained problem:

`min_alpha ||v_f - V alpha||^2`.

The equivalent quadratic-programming form is:

`min 1/2 alpha' (2 V'V) alpha - (2 v_f'V) alpha`.

The largest optimized coefficient is identified as the dominant mechanism.

## 4. Marginal distributions

Each marginal is fitted using L-moments. The candidate set contains:

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

The Cramér–von Mises statistic is used to select the best-fitting marginal.

## 5. Copula-based flood-frequency analysis

The selected marginals are transformed to uniform pseudo-observations and a
bivariate Archimedean copula is selected from:

- Clayton
- Gumbel
- Frank

The analysis evaluates:

- AND dependence;
- OR dependence;
- Kendall dependence.

For the AND case:

`P(X>x, Y>y) = 1 - F_X(x) - F_Y(y) + C(F_X(x), F_Y(y))`.

For the OR case:

`P(X>x or Y>y) = 1 - C(F_X(x), F_Y(y))`.

For the Kendall case, the Kendall distribution of the copula is evaluated and
converted to the corresponding exceedance probability.

Return period is calculated as:

`T = 1 / p`

where `p` is the relevant annual exceedance probability.
