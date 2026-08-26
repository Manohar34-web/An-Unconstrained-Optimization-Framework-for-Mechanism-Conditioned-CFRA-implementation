# ============================================================
# Configuration
# ============================================================

set.seed(101)

cfg <- list(

  # ---- Daily period ----
  start_date = as.Date("1965-01-01"),
  end_date   = as.Date("2018-12-31"),

  # The analysis removes February 29 and subsequently works with
  # 365-day hydrological-year blocks.
  remove_feb29 = TRUE,

  # Number of initial daily records removed before the 365-day blocks.
  # This preserves the original workflow supplied for the analysis.
  initial_rows_to_remove = 151,

  # ---- Compound-event definitions ----
  # Precipitation threshold for CSMP/CWSP:
  precipitation_quantile = 0.90,

  # Number of days including the selected extreme day over which
  # precipitation concurrence is checked.
  concurrence_window_days = 3,

  # ---- Flood-seasonality ----
  # Average number of days per year used in the circular transformation.
  mean_days_per_year = 365,

  # ---- Return periods ----
  return_periods = c(10, 50, 100),

  # ---- Input files ----
  precipitation_file = "data/input/precipitation_daily.csv",
  soil_moisture_file = "data/input/soil_moisture_daily.csv",
  wind_speed_file = "data/input/wind_speed_daily.csv",
  streamflow_file = "data/input/annual_max_streamflow.csv",
  streamflow_dates_file = "data/input/annual_max_streamflow_dates.csv",

  # ---- Output ----
  output_dir = "data/output"
)
