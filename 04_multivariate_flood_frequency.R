# ============================================================
# 04 - Multivariate flood-frequency analysis
#
# Marginal distributions:
# Gumbel, GEV, Log-Pearson III, Weibull, Exponential, Gamma,
# Generalized Logistic, Generalized Normal, Kappa, LogNormal III,
# Normal
#
# Copulas:
# Clayton, Gumbel, Frank
#
# Scenarios:
# AND, OR, Kendall
#
# D1/D2 -> bivariate: Streamflow + mechanism
# D3     -> trivariate: Streamflow + precipitation + antecedent soil moisture
# D4     -> trivariate: Streamflow + precipitation + wind speed
# ============================================================

source("R/utils.R")
source("config/analysis_config.R")

events <- readRDS(file.path(cfg$output_dir, "mechanism_events.rds"))
dominant <- readRDS(file.path(cfg$output_dir, "dominant_mechanisms.rds"))

streamflow <- safe_numeric_matrix(
  data.table::fread(cfg$streamflow_file, header = TRUE)
)

if (ncol(streamflow) != ncol(events$D1_value)) {
  stop("Streamflow columns do not match mechanism columns.")
}

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

fit_marginal_triplet <- function(x, y, z) {
  fx <- fit_lmom_candidates(x)
  fy <- fit_lmom_candidates(y)
  fz <- fit_lmom_candidates(z)

  list(x = fx, y = fy, z = fz)
}

fit_marginal_pair <- function(x, y) {
  list(x = fit_lmom_candidates(x), y = fit_lmom_candidates(y))
}

pseudo_pair <- function(fit, x, y) {
  u <- fit$x$cdf(x)
  v <- fit$y$cdf(y)

  keep <- is.finite(u) & is.finite(v)
  if (sum(keep) < 10L) stop("Too few valid observations for copula fitting.")

  uv <- VineCopula::pobs(cbind(
    pmin(pmax(u[keep], 1e-8), 1 - 1e-8),
    pmin(pmax(v[keep], 1e-8), 1 - 1e-8)
  ))

  list(u = uv[, 1], v = uv[, 2])
}

fit_arch_copula <- function(u, v) {
  fit_bicopula(u, v)
}

copula_cdf <- function(u, v, cop) {
  VineCopula::BiCopCDF(
    u, v,
    family = cop$family,
    par = cop$par,
    par2 = cop$par2,
    check.pars = TRUE
  )
}

# Bivariate AND/OR/Kendall return periods.
bivariate_scenarios <- function(fit_x, fit_y, cop, qx, qy) {
  ux <- pmin(pmax(fit_x$cdf(qx), 1e-8), 1 - 1e-8)
  uy <- pmin(pmax(fit_y$cdf(qy), 1e-8), 1 - 1e-8)

  Cxy <- copula_cdf(ux, uy, cop)

  p_and <- 1 - ux - uy + Cxy
  p_or  <- 1 - Cxy

  K <- kendall_distribution_archimedean(Cxy, cop$family, cop$par)
  p_kendall <- 1 - K

  c(
    AND = return_from_probability(p_and),
    OR = return_from_probability(p_or),
    Kendall = return_from_probability(p_kendall)
  )
}

# Trivariate construction follows the structure of the original analysis:
# 1. fit C(Y,Z)
# 2. construct W = C(Y,Z)
# 3. fit C(X,W)
# 4. evaluate C3(X,Y,Z) = C(X, C(Y,Z))
#
# Inclusion-exclusion is then used for the AND/OR scenarios.
trivariate_fit <- function(x, y, z, fx, fy, fz) {

  ux <- fx$cdf(x)
  uy <- fy$cdf(y)
  uz <- fz$cdf(z)

  keep <- is.finite(ux) & is.finite(uy) & is.finite(uz)

  ux <- pmin(pmax(ux[keep], 1e-8), 1 - 1e-8)
  uy <- pmin(pmax(uy[keep], 1e-8), 1 - 1e-8)
  uz <- pmin(pmax(uz[keep], 1e-8), 1 - 1e-8)

  u <- VineCopula::pobs(cbind(ux, uy, uz))
  ux <- u[, 1]; uy <- u[, 2]; uz <- u[, 3]

  cop_xy <- fit_arch_copula(ux, uy)
  cop_yz <- fit_arch_copula(uy, uz)
  cop_zx <- fit_arch_copula(uz, ux)

  C_yz <- copula_cdf(uy, uz, cop_yz)
  C_zx <- copula_cdf(uz, ux, cop_zx)
  C_xy <- copula_cdf(ux, uy, cop_xy)

  # Outer copula links X with C(Y,Z).
  cop_x_yz <- fit_arch_copula(ux, pmin(pmax(C_yz, 1e-8), 1 - 1e-8))

  list(
    cop_xy = cop_xy,
    cop_yz = cop_yz,
    cop_zx = cop_zx,
    cop_x_yz = cop_x_yz
  )
}

trivariate_scenarios <- function(fx, fy, fz, fit, qx, qy, qz) {

  ux <- pmin(pmax(fx$cdf(qx), 1e-8), 1 - 1e-8)
  uy <- pmin(pmax(fy$cdf(qy), 1e-8), 1 - 1e-8)
  uz <- pmin(pmax(fz$cdf(qz), 1e-8), 1 - 1e-8)

  Cxy <- copula_cdf(ux, uy, fit$cop_xy)
  Cyz <- copula_cdf(uy, uz, fit$cop_yz)
  Czx <- copula_cdf(uz, ux, fit$cop_zx)

  Cxyz <- copula_cdf(
    ux,
    pmin(pmax(Cyz, 1e-8), 1 - 1e-8),
    fit$cop_x_yz
  )

  # P(X>x,Y>y,Z>z)
  p_and <- 1 - ux - uy - uz +
    Cxy + Cyz + Czx - Cxyz

  # P(at least one exceeds threshold)
  p_or <- 1 - Cxyz

  # Kendall scenario uses the outer copula.
  K <- kendall_distribution_archimedean(
    Cxyz,
    fit$cop_x_yz$family,
    fit$cop_x_yz$par
  )
  p_kendall <- 1 - K

  c(
    AND = return_from_probability(p_and),
    OR = return_from_probability(p_or),
    Kendall = return_from_probability(p_kendall)
  )
}

site_results <- vector("list", ncol(streamflow))

for (j in seq_len(ncol(streamflow))) {

  cat("Flood-frequency analysis for catchment", j, "\n")

  dom <- dominant$dominant_mechanism[
    match(j, dominant$catchment)
  ]

  if (!is.character(dom) || is.na(dom)) next

  sf <- streamflow[, j]

  # ----------------------------------------------------------
  # D1: Streamflow + annual maximum precipitation
  # D2: Streamflow + 2-day precipitation
  # ----------------------------------------------------------
  if (dom == "D1_AMX_Pr") {

    mech <- events$D1_value[, j]
    keep <- is.finite(sf) & is.finite(mech)

    sf <- sf[keep]
    mech <- mech[keep]

    if (length(sf) < 10L) next

    fm <- fit_marginal_pair(sf, mech)
    uv <- pseudo_pair(fm, sf, mech)
    cop <- fit_arch_copula(uv$u, uv$v)

    rows <- lapply(cfg$return_periods, function(T) {
      p <- 1 - 1 / T
      qx <- fm$x$quantile(p)
      qy <- fm$y$quantile(p)

      rp <- bivariate_scenarios(fm$x, fm$y, cop, qx, qy)

      data.frame(
        catchment = j,
        dominant_mechanism = dom,
        analysis_type = "bivariate",
        base_return_period = T,
        streamflow_quantile = qx,
        mechanism_quantile = qy,
        third_variable_quantile = NA_real_,
        AND_JRP = rp["AND"],
        OR_JRP = rp["OR"],
        Kendall_JRP = rp["Kendall"],
        AND_probability = 1 / rp["AND"],
        OR_probability = 1 / rp["OR"],
        Kendall_probability = 1 / rp["Kendall"]
      )
    })

    dist_row <- data.frame(
      catchment = j,
      dominant_mechanism = dom,
      analysis_type = "bivariate",
      streamflow_distribution = fm$x$selected,
      second_distribution = fm$y$selected,
      third_distribution = NA_character_,
      copula_xy = cop$family,
      copula_yz = NA_integer_,
      copula_zx = NA_integer_,
      copula_outer = NA_integer_,
      copula_parameter_xy = cop$par,
      copula_parameter_yz = NA_real_,
      copula_parameter_zx = NA_real_,
      copula_parameter_outer = NA_real_
    )

  } else if (dom == "D2_AMX_MD_Pr") {

    mech <- events$D2_value[, j]
    keep <- is.finite(sf) & is.finite(mech)

    sf <- sf[keep]
    mech <- mech[keep]

    if (length(sf) < 10L) next

    fm <- fit_marginal_pair(sf, mech)
    uv <- pseudo_pair(fm, sf, mech)
    cop <- fit_arch_copula(uv$u, uv$v)

    rows <- lapply(cfg$return_periods, function(T) {
      p <- 1 - 1 / T
      qx <- fm$x$quantile(p)
      qy <- fm$y$quantile(p)

      rp <- bivariate_scenarios(fm$x, fm$y, cop, qx, qy)

      data.frame(
        catchment = j,
        dominant_mechanism = dom,
        analysis_type = "bivariate",
        base_return_period = T,
        streamflow_quantile = qx,
        mechanism_quantile = qy,
        third_variable_quantile = NA_real_,
        AND_JRP = rp["AND"],
        OR_JRP = rp["OR"],
        Kendall_JRP = rp["Kendall"],
        AND_probability = 1 / rp["AND"],
        OR_probability = 1 / rp["OR"],
        Kendall_probability = 1 / rp["Kendall"]
      )
    })

    dist_row <- data.frame(
      catchment = j,
      dominant_mechanism = dom,
      analysis_type = "bivariate",
      streamflow_distribution = fm$x$selected,
      second_distribution = fm$y$selected,
      third_distribution = NA_character_,
      copula_xy = cop$family,
      copula_yz = NA_integer_,
      copula_zx = NA_integer_,
      copula_outer = NA_integer_,
      copula_parameter_xy = cop$par,
      copula_parameter_yz = NA_real_,
      copula_parameter_zx = NA_real_,
      copula_parameter_outer = NA_real_
    )

  # ----------------------------------------------------------
  # D3: Streamflow + precipitation + antecedent soil moisture
  # ----------------------------------------------------------
  } else if (dom == "D3_CSMP") {

    x <- sf
    y <- events$D3_pr_value[, j]
    z <- events$D3_sm_value[, j]

    keep <- is.finite(x) & is.finite(y) & is.finite(z)
    x <- x[keep]; y <- y[keep]; z <- z[keep]

    if (length(x) < 10L) next

    fm <- fit_marginal_triplet(x, y, z)

    tf <- trivariate_fit(
      x, y, z,
      fm$x, fm$y, fm$z
    )

    rows <- lapply(cfg$return_periods, function(T) {
      p <- 1 - 1 / T
      qx <- fm$x$quantile(p)
      qy <- fm$y$quantile(p)
      qz <- fm$z$quantile(p)

      rp <- trivariate_scenarios(
        fm$x, fm$y, fm$z,
        tf, qx, qy, qz
      )

      data.frame(
        catchment = j,
        dominant_mechanism = dom,
        analysis_type = "trivariate_CSMP",
        base_return_period = T,
        streamflow_quantile = qx,
        mechanism_quantile = qy,
        third_variable_quantile = qz,
        AND_JRP = rp["AND"],
        OR_JRP = rp["OR"],
        Kendall_JRP = rp["Kendall"],
        AND_probability = 1 / rp["AND"],
        OR_probability = 1 / rp["OR"],
        Kendall_probability = 1 / rp["Kendall"]
      )
    })

    dist_row <- data.frame(
      catchment = j,
      dominant_mechanism = dom,
      analysis_type = "trivariate_CSMP",
      streamflow_distribution = fm$x$selected,
      second_distribution = fm$y$selected,
      third_distribution = fm$z$selected,
      copula_xy = tf$cop_xy$family,
      copula_yz = tf$cop_yz$family,
      copula_zx = tf$cop_zx$family,
      copula_outer = tf$cop_x_yz$family,
      copula_parameter_xy = tf$cop_xy$par,
      copula_parameter_yz = tf$cop_yz$par,
      copula_parameter_zx = tf$cop_zx$par,
      copula_parameter_outer = tf$cop_x_yz$par
    )

  # ----------------------------------------------------------
  # D4: Streamflow + precipitation + wind speed
  # ----------------------------------------------------------
  } else if (dom == "D4_CWSP") {

    x <- sf
    y <- events$D4_pr_value[, j]
    z <- events$D4_ws_value[, j]

    keep <- is.finite(x) & is.finite(y) & is.finite(z)
    x <- x[keep]; y <- y[keep]; z <- z[keep]

    if (length(x) < 10L) next

    fm <- fit_marginal_triplet(x, y, z)

    tf <- trivariate_fit(
      x, y, z,
      fm$x, fm$y, fm$z
    )

    rows <- lapply(cfg$return_periods, function(T) {
      p <- 1 - 1 / T
      qx <- fm$x$quantile(p)
      qy <- fm$y$quantile(p)
      qz <- fm$z$quantile(p)

      rp <- trivariate_scenarios(
        fm$x, fm$y, fm$z,
        tf, qx, qy, qz
      )

      data.frame(
        catchment = j,
        dominant_mechanism = dom,
        analysis_type = "trivariate_CWSP",
        base_return_period = T,
        streamflow_quantile = qx,
        mechanism_quantile = qy,
        third_variable_quantile = qz,
        AND_JRP = rp["AND"],
        OR_JRP = rp["OR"],
        Kendall_JRP = rp["Kendall"],
        AND_probability = 1 / rp["AND"],
        OR_probability = 1 / rp["OR"],
        Kendall_probability = 1 / rp["Kendall"]
      )
    })

    dist_row <- data.frame(
      catchment = j,
      dominant_mechanism = dom,
      analysis_type = "trivariate_CWSP",
      streamflow_distribution = fm$x$selected,
      second_distribution = fm$y$selected,
      third_distribution = fm$z$selected,
      copula_xy = tf$cop_xy$family,
      copula_yz = tf$cop_yz$family,
      copula_zx = tf$cop_zx$family,
      copula_outer = tf$cop_x_yz$family,
      copula_parameter_xy = tf$cop_xy$par,
      copula_parameter_yz = tf$cop_yz$par,
      copula_parameter_zx = tf$cop_zx$par,
      copula_parameter_outer = tf$cop_x_yz$par
    )

  } else {
    next
  }

  site_results[[j]] <- list(
    distributions = dist_row,
    return_periods = data.table::rbindlist(rows, fill = TRUE)
  )
}

site_results <- site_results[!vapply(site_results, is.null, logical(1))]

dist_table <- data.table::rbindlist(
  lapply(site_results, `[[`, "distributions"),
  fill = TRUE
)

rp_table <- data.table::rbindlist(
  lapply(site_results, `[[`, "return_periods"),
  fill = TRUE
)

data.table::fwrite(
  dist_table,
  file.path(cfg$output_dir, "selected_marginal_distributions_and_copulas.csv")
)

data.table::fwrite(
  rp_table,
  file.path(cfg$output_dir, "multivariate_flood_frequency_results.csv")
)

saveRDS(
  list(
    distributions = dist_table,
    return_periods = rp_table
  ),
  file.path(cfg$output_dir, "multivariate_flood_frequency_results.rds")
)
