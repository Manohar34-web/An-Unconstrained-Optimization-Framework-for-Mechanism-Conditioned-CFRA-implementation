# ============================================================
# Utility functions
# ============================================================

safe_numeric_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x
}

remove_feb29 <- function(x, dates) {
  keep <- !(lubridate::month(dates) == 2 & lubridate::day(dates) == 29)
  list(data = x[keep, , drop = FALSE], dates = dates[keep])
}

trim_initial_rows <- function(x, n = 0L) {
  if (n <= 0L) return(x)
  if (n >= nrow(x)) stop("initial_rows_to_remove is >= number of rows.")
  x[(n + 1L):nrow(x), , drop = FALSE]
}

make_365_blocks <- function(x) {
  n_complete <- floor(nrow(x) / 365)
  if (n_complete < 1L) stop("Not enough observations for a complete 365-day block.")
  x[seq_len(n_complete * 365L), , drop = FALSE]
}

annual_max_and_date <- function(x, dates) {
  x <- safe_numeric_matrix(x)
  n_years <- floor(nrow(x) / 365)
  out_value <- matrix(NA_real_, n_years, ncol(x))
  out_date  <- matrix(NA_integer_, n_years, ncol(x))

  for (j in seq_len(ncol(x))) {
    z <- x[seq_len(n_years * 365L), j]
    blocks <- split(z, rep(seq_len(n_years), each = 365L))
    for (yy in seq_len(n_years)) {
      b <- blocks[[yy]]
      if (all(is.na(b))) next
      ii <- which.max(replace(b, is.na(b), -Inf))
      if (length(ii) == 0L || !is.finite(b[ii])) next
      out_value[yy, j] <- b[ii]
      out_date[yy, j] <- ii
    }
  }
  list(value = out_value, doy = out_date)
}

# Circular statistics following the formulation used in the manuscript.
circular_components <- function(doy, m = 365) {
  z <- doy[is.finite(doy)]
  if (!length(z)) {
    return(c(xbar = NA_real_, ybar = NA_real_,
             R = NA_real_, mean_doy = NA_real_, n = 0L))
  }

  theta <- 2 * pi * z / m
  xbar <- mean(cos(theta))
  ybar <- mean(sin(theta))
  R <- sqrt(xbar^2 + ybar^2)

  ang <- atan2(ybar, xbar)
  if (ang < 0) ang <- ang + 2 * pi
  mean_doy <- ang * m / (2 * pi)

  c(xbar = xbar, ybar = ybar, R = R,
    mean_doy = mean_doy, n = length(z))
}

# Unconstrained quadratic minimization:
# min_alpha ||v_f - V %*% alpha||^2
#
# This is mathematically equivalent to the unconstrained QP used in the
# methodology:
# 1/2 alpha' (2 V'V) alpha - (2 v_f'V) alpha.
unconstrained_qp <- function(V, vf) {
  Dmat <- 2 * crossprod(V)
  dvec <- 2 * crossprod(V, vf)

  alpha <- tryCatch(
    solve(Dmat, dvec),
    error = function(e) qr.solve(Dmat, dvec)
  )

  as.numeric(alpha)
}

# ------------------------------------------------------------
# L-moment marginal distributions
# ------------------------------------------------------------

fit_lmom_candidates <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 10L) stop("Too few observations for marginal fitting.")

  lm <- lmom::samlmu(x)

  fits <- list(
    Gumbel = lmom::pelgum(lm),
    GEV = lmom::pelgev(lm),
    LogPearsonIII = lmom::pelpe3(lm),
    Weibull = lmom::pelwei(lm),
    Exponential = lmom::pelexp(lm),
    Gamma = lmom::pelgam(lm),
    GeneralizedLogistic = lmom::pelglo(lm),
    GeneralizedNormal = lmom::pelgno(lm),
    Kappa = lmom::pelkap(lm),
    LogNormalIII = lmom::pelln3(lm),
    Normal = lmom::pelnor(lm)
  )

  cdfs <- list(
    Gumbel = function(z) lmom::cdfgum(z, fits$Gumbel),
    GEV = function(z) lmom::cdfgev(z, fits$GEV),
    LogPearsonIII = function(z) lmom::cdfpe3(z, fits$LogPearsonIII),
    Weibull = function(z) lmom::cdfwei(z, fits$Weibull),
    Exponential = function(z) lmom::cdfexp(z, fits$Exponential),
    Gamma = function(z) lmom::cdfgam(z, fits$Gamma),
    GeneralizedLogistic = function(z) lmom::cdfglo(z, fits$GeneralizedLogistic),
    GeneralizedNormal = function(z) lmom::cdfgno(z, fits$GeneralizedNormal),
    Kappa = function(z) lmom::cdfkap(z, fits$Kappa),
    LogNormalIII = function(z) lmom::cdfln3(z, fits$LogNormalIII),
    Normal = function(z) lmom::cdfnor(z, fits$Normal)
  )

  quantiles <- list(
    Gumbel = function(p) lmom::quagum(p, fits$Gumbel),
    GEV = function(p) lmom::quagev(p, fits$GEV),
    LogPearsonIII = function(p) lmom::quape3(p, fits$LogPearsonIII),
    Weibull = function(p) lmom::quawei(p, fits$Weibull),
    Exponential = function(p) lmom::quaexp(p, fits$Exponential),
    Gamma = function(p) lmom::quagam(p, fits$Gamma),
    GeneralizedLogistic = function(p) lmom::quaglo(p, fits$GeneralizedLogistic),
    GeneralizedNormal = function(p) lmom::quagno(p, fits$GeneralizedNormal),
    Kappa = function(p) lmom::quakap(p, fits$Kappa),
    LogNormalIII = function(p) lmom::qualn3(p, fits$LogNormalIII),
    Normal = function(p) lmom::quanor(p, fits$Normal)
  )

  cvm <- vapply(names(fits), function(nm) {
    ans <- tryCatch(
      goftest::cvm.test(x, null = cdfs[[nm]]),
      error = function(e) NULL
    )
    if (is.null(ans)) Inf else as.numeric(ans$statistic)
  }, numeric(1))

  selected <- names(which.min(cvm))

  list(
    selected = selected,
    fit = fits[[selected]],
    cdf = cdfs[[selected]],
    quantile = quantiles[[selected]],
    statistics = cvm,
    fits = fits
  )
}

# Pseudo-observations based on fitted marginal CDFs.
to_pobs <- function(x) {
  u <- as.numeric(x)
  u <- pmin(pmax(u, 1e-8), 1 - 1e-8)
  u
}

# Kendall distribution for the Archimedean copulas used in the analysis.
kendall_distribution_archimedean <- function(t, family, theta) {
  t <- pmin(pmax(t, 1e-10), 1 - 1e-10)

  if (family == 3) { # Clayton
    phi <- function(p) (p^(-theta) - 1) / theta
    phi_prime <- function(p) -p^(-theta - 1)
  } else if (family == 4) { # Gumbel
    phi <- function(p) (-log(p))^theta
    phi_prime <- function(p) -theta * (-log(p)^(theta - 1)) / p
  } else if (family == 5) { # Frank
    phi <- function(p) -log((exp(-theta * p) - 1) /
                              (exp(-theta) - 1))
    phi_prime <- function(p) theta * exp(-theta * p) /
                              (exp(-theta * p) - 1)
  } else {
    stop("Kendall formulation implemented for Clayton, Gumbel and Frank only.")
  }

  out <- t - phi(t) / phi_prime(t)
  pmin(pmax(out, 0), 1)
}

fit_bicopula <- function(u, v) {
  # Restrict to copulas for which the analytical Kendall distribution
  # used below is available.
  VineCopula::BiCopSelect(
    u, v,
    familyset = c(3, 4, 5),
    rotations = FALSE,
    selectioncrit = "AIC"
  )
}

marginal_probability <- function(data, q) {
  sum(data <= q, na.rm = TRUE) / (sum(is.finite(data)) + 1)
}

return_from_probability <- function(p) {
  p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
  1 / (1 - p)
}
