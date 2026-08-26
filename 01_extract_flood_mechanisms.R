# ============================================================
# 01 - Flood-generating mechanism event-date extraction
# ============================================================

source("R/utils.R")
source("config/analysis_config.R")

dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

read_daily <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  data.table::fread(path, header = TRUE)
}

# ---- Dates ----
dates <- seq(cfg$start_date, cfg$end_date, by = "day")

if (length(dates) != nrow(read_daily(cfg$precipitation_file))) {
  stop("Configured date range does not match precipitation row count.")
}

# ---- Read variables ----
P  <- safe_numeric_matrix(read_daily(cfg$precipitation_file))
SM <- safe_numeric_matrix(read_daily(cfg$soil_moisture_file))
WS <- safe_numeric_matrix(read_daily(cfg$wind_speed_file))

if (!all(dim(P) == dim(SM)) || !all(dim(P) == dim(WS))) {
  stop("Precipitation, soil moisture and wind-speed matrices must have identical dimensions.")
}

# Remove February 29 so all annual blocks contain 365 days.
if (cfg$remove_feb29) {
  tmp <- remove_feb29(P, dates); P <- tmp$data; dates <- tmp$dates
  SM <- SM[!(lubridate::month(seq(cfg$start_date, cfg$end_date, by="day")) == 2 &
             lubridate::day(seq(cfg$start_date, cfg$end_date, by="day")) == 29), , drop=FALSE]
  WS <- WS[!(lubridate::month(seq(cfg$start_date, cfg$end_date, by="day")) == 2 &
             lubridate::day(seq(cfg$start_date, cfg$start_date + length(dates)-1, by="day")) == 29), , drop=FALSE]
}

# Re-align after leap-day removal using the final date vector.
all_dates <- seq(cfg$start_date, cfg$end_date, by = "day")
keep_dates <- !(month(all_dates) == 2 & day(all_dates) == 29)
P  <- P[keep_dates, , drop = FALSE]
SM <- SM[keep_dates, , drop = FALSE]
WS <- WS[keep_dates, , drop = FALSE]
dates <- all_dates[keep_dates]

# Preserve the original 151-record trimming before hydrological-year blocks.
P  <- trim_initial_rows(P, cfg$initial_rows_to_remove)
SM <- trim_initial_rows(SM, cfg$initial_rows_to_remove)
WS <- trim_initial_rows(WS, cfg$initial_rows_to_remove)
dates <- dates[(cfg$initial_rows_to_remove + 1L):length(dates)]

P  <- make_365_blocks(P)
SM <- make_365_blocks(SM)
WS <- make_365_blocks(WS)
dates <- dates[seq_len(nrow(P))]

n_years <- nrow(P) / 365L
n_sites <- ncol(P)

# ------------------------------------------------------------
# D1: Annual maximum precipitation
# ------------------------------------------------------------
D1_value <- D1_doy <- matrix(NA_real_, n_years, n_sites)

# ------------------------------------------------------------
# D2: 2-day precipitation maximum
# ------------------------------------------------------------
D2_value <- D2_doy <- matrix(NA_real_, n_years, n_sites)

# ------------------------------------------------------------
# D3: Compound soil-moisture / precipitation
#
# select the highest soil-moisture day for which precipitation exceeds the annual
# positive-precipitation 90th percentile within the prescribed
# forward concurrence window.
# ------------------------------------------------------------
D3_doy <- D3_sm_value <- D3_pr_value <- matrix(NA_real_, n_years, n_sites)

# ------------------------------------------------------------
# D4: Compound wind-speed / precipitation
# ------------------------------------------------------------
D4_doy <- D4_ws_value <- D4_pr_value <- matrix(NA_real_, n_years, n_sites)

for (j in seq_len(n_sites)) {
  cat("Extracting mechanisms for catchment", j, "of", n_sites, "\n")

  p <- P[, j]
  sm <- SM[, j]
  ws <- WS[, j]

  for (yy in seq_len(n_years)) {
    idx <- ((yy - 1L) * 365L + 1L):(yy * 365L)

    pb  <- p[idx]
    smb <- sm[idx]
    wsb <- ws[idx]

    # D1
    if (any(is.finite(pb))) {
      i1 <- which.max(replace(pb, !is.finite(pb), -Inf))
      D1_value[yy, j] <- pb[i1]
      D1_doy[yy, j] <- i1
    }

    # D2
    if (sum(is.finite(pb)) >= 2L) {
      two_day <- zoo::rollapply(pb, 2, sum, fill = NA, align = "right",
                                na.rm = FALSE)
      i2 <- which.max(replace(two_day, !is.finite(two_day), -Inf))
      if (is.finite(two_day[i2])) {
        D2_value[yy, j] <- two_day[i2]
        D2_doy[yy, j] <- i2
      }
    }

    # Annual positive-precipitation threshold
    ppos <- pb[is.finite(pb) & pb > 0]
    if (!length(ppos)) next
    pthr <- as.numeric(stats::quantile(ppos, cfg$precipitation_quantile,
                                        names = FALSE, na.rm = TRUE))

    # D3: high antecedent soil moisture + subsequent high precipitation
    sm_order <- order(smb, decreasing = TRUE, na.last = NA)
    for (d in sm_order) {
      if (!is.finite(smb[d])) next
      e <- min(d + cfg$concurrence_window_days - 1L, 365L)
      prw <- pb[d:e]
      hit <- which(is.finite(prw) & prw >= pthr)
      if (length(hit)) {
        dd <- d + hit[1L] - 1L
        D3_doy[yy, j] <- d
        D3_sm_value[yy, j] <- smb[d]
        D3_pr_value[yy, j] <- pb[dd]
        break
      }
    }

    # D4: highest wind extreme with concurrent/subsequent high precipitation
    ws_order <- order(wsb, decreasing = TRUE, na.last = NA)
    for (d in ws_order) {
      if (!is.finite(wsb[d])) next
      e <- min(d + cfg$concurrence_window_days - 1L, 365L)
      prw <- pb[d:e]
      hit <- which(is.finite(prw) & prw >= pthr)
      if (length(hit)) {
        dd <- d + hit[1L] - 1L
        D4_doy[yy, j] <- d
        D4_ws_value[yy, j] <- wsb[d]
        D4_pr_value[yy, j] <- pb[dd]
        break
      }
    }
  }
}

colnames(D1_value) <- colnames(D1_doy) <- colnames(P)
colnames(D2_value) <- colnames(D2_doy) <- colnames(P)
colnames(D3_doy) <- colnames(D3_sm_value) <- colnames(D3_pr_value) <- colnames(P)
colnames(D4_doy) <- colnames(D4_ws_value) <- colnames(D4_pr_value) <- colnames(P)

write_matrix <- function(x, name) {
  data.table::fwrite(as.data.table(x), file.path(cfg$output_dir, name))
}

write_matrix(D1_value, "D1_AMX_Pr_value.csv")
write_matrix(D1_doy, "D1_AMX_Pr_date.csv")
write_matrix(D2_value, "D2_AMX_MD_Pr_value.csv")
write_matrix(D2_doy, "D2_AMX_MD_Pr_date.csv")
write_matrix(D3_doy, "D3_CSMP_date.csv")
write_matrix(D3_sm_value, "D3_CSMP_soil_moisture.csv")
write_matrix(D3_pr_value, "D3_CSMP_precipitation.csv")
write_matrix(D4_doy, "D4_CWSP_date.csv")
write_matrix(D4_ws_value, "D4_CWSP_wind_speed.csv")
write_matrix(D4_pr_value, "D4_CWSP_precipitation.csv")

saveRDS(list(
  D1_value = D1_value, D1_doy = D1_doy,
  D2_value = D2_value, D2_doy = D2_doy,
  D3_doy = D3_doy, D3_sm_value = D3_sm_value, D3_pr_value = D3_pr_value,
  D4_doy = D4_doy, D4_ws_value = D4_ws_value, D4_pr_value = D4_pr_value
), file.path(cfg$output_dir, "mechanism_events.rds"))
