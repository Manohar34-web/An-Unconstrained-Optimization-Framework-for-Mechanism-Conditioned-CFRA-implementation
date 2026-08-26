# ============================================================
# 02 - Circular seasonality statistics
# ============================================================

source("R/utils.R")
source("config/analysis_config.R")

events <- readRDS(file.path(cfg$output_dir, "mechanism_events.rds"))

read_streamflow_dates <- function(path) {
  x <- data.table::fread(path, header = TRUE)
  as.matrix(x)
}

# Flood dates are supplied as annual maximum flood dates.
flood_date <- read_streamflow_dates(cfg$streamflow_dates_file)

n_sites <- ncol(events$D1_doy)

get_stats_matrix <- function(x) {
  out <- matrix(NA_real_, n_sites, 5)
  colnames(out) <- c("xbar", "ybar", "R", "mean_doy", "n")
  for (j in seq_len(n_sites)) {
    out[j, ] <- circular_components(x[, j], cfg$mean_days_per_year)
  }
  out
}

stats_D1 <- get_stats_matrix(events$D1_doy)
stats_D2 <- get_stats_matrix(events$D2_doy)
stats_D3 <- get_stats_matrix(events$D3_doy)
stats_D4 <- get_stats_matrix(events$D4_doy)
stats_F  <- get_stats_matrix(flood_date)

write_stats <- function(x, name) {
  data.table::fwrite(
    data.table::as.data.table(x, keep.rownames = "catchment"),
    file.path(cfg$output_dir, name)
  )
}

write_stats(stats_D1, "circular_D1_AMX_Pr.csv")
write_stats(stats_D2, "circular_D2_AMX_MD_Pr.csv")
write_stats(stats_D3, "circular_D3_CSMP.csv")
write_stats(stats_D4, "circular_D4_CWSP.csv")
write_stats(stats_F,  "circular_flood.csv")

saveRDS(
  list(D1 = stats_D1, D2 = stats_D2, D3 = stats_D3,
       D4 = stats_D4, Flood = stats_F),
  file.path(cfg$output_dir, "circular_statistics.rds")
)
