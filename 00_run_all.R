# ============================================================
# Master workflow
# ============================================================

set.seed(101)

required <- c(
  "data.table", "lubridate", "zoo", "hydroTSM",
  "hydromad", "lmom", "goftest", "VineCopula"
)

missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing)) {
  stop(
    "Install the following packages before running the analysis: ",
    paste(missing, collapse = ", ")
  )
}

# Run sequentially because each stage produces inputs for the next.
source("R/01_extract_flood_mechanisms.R")
source("R/02_circular_statistics.R")
source("R/03_dominant_mechanism_optimization.R")
source("R/04_multivariate_flood_frequency.R")

message("Complete workflow finished.")
