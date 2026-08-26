# ============================================================
# 03 - Dominant flood-generating mechanism
#      Unconstrained quadratic optimization
# ============================================================

source("R/utils.R")
source("config/analysis_config.R")

cs <- readRDS(file.path(cfg$output_dir, "circular_statistics.rds"))

mechanism_names <- c("D1_AMX_Pr", "D2_AMX_MD_Pr", "D3_CSMP", "D4_CWSP")
n_sites <- nrow(cs$Flood)

results <- vector("list", n_sites)

for (j in seq_len(n_sites)) {
  vf <- c(cs$Flood[j, "xbar"], cs$Flood[j, "ybar"])

  V <- cbind(
    c(cs$D1[j, "xbar"], cs$D1[j, "ybar"]),
    c(cs$D2[j, "xbar"], cs$D2[j, "ybar"]),
    c(cs$D3[j, "xbar"], cs$D3[j, "ybar"]),
    c(cs$D4[j, "xbar"], cs$D4[j, "ybar"])
  )

  ok <- is.finite(vf) & all(is.finite(V))
  if (!ok) {
    results[[j]] <- data.frame(
      catchment = j,
      alpha_D1 = NA_real_, alpha_D2 = NA_real_,
      alpha_D3 = NA_real_, alpha_D4 = NA_real_,
      dominant_mechanism = NA_character_,
      objective = NA_real_
    )
    next
  }

  alpha <- unconstrained_qp(V, vf)
  reconstructed <- as.numeric(V %*% alpha)
  objective <- sum((vf - reconstructed)^2)

  # The manuscript interprets the largest optimized alpha as dominant.
  dom <- mechanism_names[which.max(alpha)]

  results[[j]] <- data.frame(
    catchment = j,
    alpha_D1 = alpha[1],
    alpha_D2 = alpha[2],
    alpha_D3 = alpha[3],
    alpha_D4 = alpha[4],
    dominant_mechanism = dom,
    objective = objective
  )
}

dominant <- data.table::rbindlist(results, fill = TRUE)

data.table::fwrite(
  dominant,
  file.path(cfg$output_dir, "dominant_mechanisms.csv")
)

saveRDS(dominant, file.path(cfg$output_dir, "dominant_mechanisms.rds"))
