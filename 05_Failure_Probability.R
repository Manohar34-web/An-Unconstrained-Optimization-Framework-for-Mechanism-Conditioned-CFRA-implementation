# ============================================================
# FAILURE PROBABILITY ANALYSIS
# ============================================================
#
# Purpose:
# Calculate the probability that at least one design-level
# flood event occurs during a specified infrastructure
# service life.
#
# Methodology:
#
# FP_T = 1 - (1 - P)^T
#
# where:
#   FP_T = failure probability over T years
#   P    = annual exceedance probability
#   T    = design life in years
#
# Three multivariate hazard formulations are considered:
#   1. AND
#   2. OR
#   3. Kendall/Copula
#
# ============================================================


# ------------------------------------------------------------
# 1. INITIALIZATION
# ------------------------------------------------------------

rm(list = ls())

library(data.table)


# ------------------------------------------------------------
# 2. USER SETTINGS
# ------------------------------------------------------------

# Design lives considered in the analysis
design_lives <- c(10, 25, 50, 100)


# ------------------------------------------------------------
# 3. READ MULTIVARIATE RETURN PERIOD RESULTS
# ------------------------------------------------------------

AND_JRP <- fread("AND_JRP.csv")
OR_JRP  <- fread("OR_JRP.csv")
Ken_JRP <- fread("Ken_JRP.csv")


# Convert to matrices
AND_JRP <- as.matrix(AND_JRP)
OR_JRP  <- as.matrix(OR_JRP)
Ken_JRP <- as.matrix(Ken_JRP)


# ------------------------------------------------------------
# 4. FAILURE-PROBABILITY FUNCTION
# ------------------------------------------------------------
#
# Equation S5:
#
# FP_T = 1 - (1 - P)^T
#
# ------------------------------------------------------------

failure_probability <- function(P, T) {
  
  # Avoid numerical problems
  P <- pmin(pmax(P, 0), 1)
  
  FP <- 1 - (1 - P)^T
  
  return(FP)
}


# ------------------------------------------------------------
# 5. NUMBER OF BASINS / STATIONS
# ------------------------------------------------------------

n_sites <- nrow(AND_JRP)


# ------------------------------------------------------------
# 6. INITIALIZE OUTPUT MATRICES
# ------------------------------------------------------------

FP_AND <- matrix(
  NA_real_,
  nrow = n_sites,
  ncol = length(design_lives)
)

FP_OR <- matrix(
  NA_real_,
  nrow = n_sites,
  ncol = length(design_lives)
)

FP_Ken <- matrix(
  NA_real_,
  nrow = n_sites,
  ncol = length(design_lives)
)


# Column names
colnames(FP_AND) <- paste0("FP_", design_lives, "yr")
colnames(FP_OR)  <- paste0("FP_", design_lives, "yr")
colnames(FP_Ken) <- paste0("FP_", design_lives, "yr")


# ------------------------------------------------------------
# 7. FAILURE PROBABILITY CALCULATION
# ------------------------------------------------------------

for (k in 1:n_sites) {
  
  cat("Processing site:", k, "of", n_sites, "\n")
  
  
  # ----------------------------------------------------------
  # 7.1 AND SCENARIO
  # ----------------------------------------------------------
  #
  # Annual exceedance probability:
  #
  # P_AND = 1 / JRP_AND
  #
  # Equation S8
  #
  
  
  P_AND <- 1 / AND_JRP[k, ]
  
  
  # Failure probability for each design life
  #
  # Equation S9
  #
  
  FP_AND[k, ] <- sapply(
    design_lives,
    function(T) failure_probability(P_AND, T)
  )
  
  
  # ----------------------------------------------------------
  # 7.2 OR SCENARIO
  # ----------------------------------------------------------
  #
  # Annual exceedance probability:
  #
  # P_OR = 1 / JRP_OR
  #
  # Equation S6
  #
  
  P_OR <- 1 / OR_JRP[k, ]
  
  
  # Failure probability
  #
  # Equation S7
  #
  
  FP_OR[k, ] <- sapply(
    design_lives,
    function(T) failure_probability(P_OR, T)
  )
  
  
  # ----------------------------------------------------------
  # 7.3 KENDALL / COPULA SCENARIO
  # ----------------------------------------------------------
  #
  # Annual exceedance probability:
  #
  # P_C = 1 / JRP_Ken
  #
  # Failure probability:
  #
  # Equation S10
  #
  
  P_C <- 1 / Ken_JRP[k, ]
  
  
  FP_Ken[k, ] <- sapply(
    design_lives,
    function(T) failure_probability(P_C, T)
  )
  
}


# ------------------------------------------------------------
# 8. ROUND RESULTS FOR OUTPUT
# ------------------------------------------------------------

FP_AND_Output <- round(FP_AND, 6)
FP_OR_Output  <- round(FP_OR, 6)
FP_Ken_Output <- round(FP_Ken, 6)


# ------------------------------------------------------------
# 9. WRITE OUTPUT FILES
# ------------------------------------------------------------

fwrite(
  as.data.table(FP_AND_Output),
  "Failure_Probability_AND.csv"
)

fwrite(
  as.data.table(FP_OR_Output),
  "Failure_Probability_OR.csv"
)

fwrite(
  as.data.table(FP_Ken_Output),
  "Failure_Probability_Kendall.csv"
)


# ------------------------------------------------------------
# 10. COMBINED FAILURE-PROBABILITY OUTPUT
# ------------------------------------------------------------

FP_Combined <- data.table(
  Site = 1:n_sites
)


# AND
for (i in seq_along(design_lives)) {
  
  FP_Combined[
    ,
    paste0("AND_FP_", design_lives[i], "yr") :=
      FP_AND_Output[, i]
  ]
  
}


# OR
for (i in seq_along(design_lives)) {
  
  FP_Combined[
    ,
    paste0("OR_FP_", design_lives[i], "yr") :=
      FP_OR_Output[, i]
  ]
  
}


# Kendall / Copula
for (i in seq_along(design_lives)) {
  
  FP_Combined[
    ,
    paste0("Kendall_FP_", design_lives[i], "yr") :=
      FP_Ken_Output[, i]
  ]
  
}


# Save combined results
fwrite(
  FP_Combined,
  "Failure_Probability_All_Scenarios.csv"
)


# ------------------------------------------------------------
# 11. SUMMARY
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("FAILURE PROBABILITY ANALYSIS COMPLETED\n")
cat("============================================================\n")
cat("Design lives:", paste(design_lives, collapse = ", "), "years\n")
cat("Number of sites:", n_sites, "\n")
cat("\nOutput files:\n")
cat("  - Failure_Probability_AND.csv\n")
cat("  - Failure_Probability_OR.csv\n")
cat("  - Failure_Probability_Kendall.csv\n")
cat("  - Failure_Probability_All_Scenarios.csv\n")
cat("============================================================\n")