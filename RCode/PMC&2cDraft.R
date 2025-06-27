# ---- LIBRARIES ----
library(data.table)
library(dplyr)
library(readxl)
library(tibble)
library(ggplot2)
library(Deriv)
library(Ryacas)
library(numDeriv)
library(rootSolve)

rm(list=ls())

# ---- PATHS ----
dir_to_script <- getwd()

# For PMC/ANOVA
pmc_input_file <- file.path(dir_to_script, "..", "GC", "cosall.xlsx")
if (!file.exists(pmc_input_file)) stop("cosall.xlsx not found at: ", pmc_input_file)

# For Cosinor-2c
cosinor_input_folder <- file.path(dir_to_script, "..", "01-input", "AllTwin2Day")
if (!dir.exists(cosinor_input_folder)) stop("AllTwin2Day folder not found: ", cosinor_input_folder)

# Output
output_path <- file.path(dir_to_script, "..", "02-output")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

# ---- PMC ANALYSIS (cosall.xlsx) ----
twinsdata <- read_excel(pmc_input_file)
groupedtwinsdata <- twinsdata %>% group_by(PID)

runanova <- function(subpopulation, rhythmcharacteristics, data, formula_str) {
  model_formula <- as.formula(formula_str)
  a <- aov(model_formula, data=data)
  a_summary <- summary(a)
  MSBetweenGroups <- a_summary[[1]][1,3]
  MSWithinGroups <- a_summary[[1]][2,3]
  rl <- (MSBetweenGroups - MSWithinGroups) / (MSBetweenGroups + MSWithinGroups)
  F_stat <- MSBetweenGroups/MSWithinGroups
  p_value <- a_summary[[1]][1,5]
  k <- length(unique(data$PID))
  results <- list(
    subpopulation = subpopulation,
    metric = rhythmcharacteristics,
    ms_between = MSBetweenGroups,
    ms_within = MSWithinGroups,
    icc = rl,
    f_stat = F_stat,
    p_value = p_value,
    n_pairs = k,
    residuals = a$residuals
  )
  return(results)
}

runAmplitudeAndMesorAnovas <- function(subpopulation, data) {
  amplitude_results <- runanova(subpopulation, "Amplitude", data, "A24h ~ PID")
  mesor_results <- runanova(subpopulation, "MESOR", data, "M24h ~ PID")
  return(list(amplitude = amplitude_results, mesor = mesor_results))
}

all_data_results <- runAmplitudeAndMesorAnovas("AllData", groupedtwinsdata)
low_bw_data <- groupedtwinsdata[groupedtwinsdata$Gp=="1L",]
low_bw_results <- runAmplitudeAndMesorAnovas("LowBodyWeight", low_bw_data)
medium_bw_data <- groupedtwinsdata[groupedtwinsdata$Gp=="2M",]
medium_bw_results <- runAmplitudeAndMesorAnovas("MediumBodyWeight", medium_bw_data)
high_bw_data <- groupedtwinsdata[groupedtwinsdata$Gp=="3H",]
high_bw_results <- runAmplitudeAndMesorAnovas("HighBodyWeight", high_bw_data)

consolidated_results <- data.frame(
  Subpopulation = c(
    "AllData", "AllData",
    "LowBodyWeight", "LowBodyWeight",
    "MediumBodyWeight", "MediumBodyWeight",
    "HighBodyWeight", "HighBodyWeight"
  ),
  Metric = c(
    "Amplitude", "MESOR",
    "Amplitude", "MESOR",
    "Amplitude", "MESOR",
    "Amplitude", "MESOR"
  ),
  MS_Between = c(
    all_data_results$amplitude$ms_between, all_data_results$mesor$ms_between,
    low_bw_results$amplitude$ms_between, low_bw_results$mesor$ms_between,
    medium_bw_results$amplitude$ms_between, medium_bw_results$mesor$ms_between,
    high_bw_results$amplitude$ms_between, high_bw_results$mesor$ms_between
  ),
  MS_Within = c(
    all_data_results$amplitude$ms_within, all_data_results$mesor$ms_within,
    low_bw_results$amplitude$ms_within, low_bw_results$mesor$ms_within,
    medium_bw_results$amplitude$ms_within, medium_bw_results$mesor$ms_within,
    high_bw_results$amplitude$ms_within, high_bw_results$mesor$ms_within
  ),
  ICC = c(
    all_data_results$amplitude$icc, all_data_results$mesor$icc,
    low_bw_results$amplitude$icc, low_bw_results$mesor$icc,
    medium_bw_results$amplitude$icc, medium_bw_results$mesor$icc,
    high_bw_results$amplitude$icc, high_bw_results$mesor$icc
  ),
  F_statistic = c(
    all_data_results$amplitude$f_stat, all_data_results$mesor$f_stat,
    low_bw_results$amplitude$f_stat, low_bw_results$mesor$f_stat,
    medium_bw_results$amplitude$f_stat, medium_bw_results$mesor$f_stat,
    high_bw_results$amplitude$f_stat, high_bw_results$mesor$f_stat
  ),
  P_value = c(
    all_data_results$amplitude$p_value, all_data_results$mesor$p_value,
    low_bw_results$amplitude$p_value, low_bw_results$mesor$p_value,
    medium_bw_results$amplitude$p_value, medium_bw_results$mesor$p_value,
    high_bw_results$amplitude$p_value, high_bw_results$mesor$p_value
  ),
  N_Pairs = c(
    all_data_results$amplitude$n_pairs, all_data_results$mesor$n_pairs,
    low_bw_results$amplitude$n_pairs, low_bw_results$mesor$n_pairs,
    medium_bw_results$amplitude$n_pairs, medium_bw_results$mesor$n_pairs,
    high_bw_results$amplitude$n_pairs, high_bw_results$mesor$n_pairs
  )
)

write.csv(consolidated_results, file.path(output_path, "consolidated_analysis_results.csv"), row.names = FALSE)

# ---- PMC PLOTS ----
twinsdata$GroupLabel <- factor(twinsdata$Gp, levels = c("1L", "2M", "3H"), labels = c("Low BW", "Medium BW", "High BW"))

png(file.path(output_path, "combined_a24h_boxplot.png"), width=800, height=600)
boxplot(A24h ~ GroupLabel, data=twinsdata,
        main="Amplitude (A24h) by Weight Group",
        xlab="Weight Group", ylab="Amplitude (A24h)",
        col=c("lightcoral", "lightgreen", "lightblue"))
dev.off()

png(file.path(output_path, "combined_m24h_boxplot.png"), width=800, height=600)
boxplot(M24h ~ GroupLabel, data=twinsdata,
        main="MESOR (M24h) by Weight Group",
        xlab="Weight Group", ylab="MESOR (M24h)",
        col=c("lightcoral", "lightgreen", "lightblue"))
dev.off()

# Residual histograms
png(file.path(output_path, "combined_a24h_residuals.png"), width=1000, height=800)
par(mfrow=c(2,2))
hist(all_data_results$amplitude$residuals, main="A24h Residuals - All Data", xlab="Residuals", col="lightgray", border="black")
hist(low_bw_results$amplitude$residuals, main="A24h Residuals - Low BW", xlab="Residuals", col="lightcoral", border="black")
hist(medium_bw_results$amplitude$residuals, main="A24h Residuals - Medium BW", xlab="Residuals", col="lightgreen", border="black")
hist(high_bw_results$amplitude$residuals, main="A24h Residuals - High BW", xlab="Residuals", col="lightblue", border="black")
par(mfrow=c(1,1))
dev.off()

png(file.path(output_path, "combined_m24h_residuals.png"), width=1000, height=800)
par(mfrow=c(2,2))
hist(all_data_results$mesor$residuals, main="M24h Residuals - All Data", xlab="Residuals", col="lightgray", border="black")
hist(low_bw_results$mesor$residuals, main="M24h Residuals - Low BW", xlab="Residuals", col="lightcoral", border="black")
hist(medium_bw_results$mesor$residuals, main="M24h Residuals - Medium BW", xlab="Residuals", col="lightgreen", border="black")
hist(high_bw_results$mesor$residuals, main="M24h Residuals - High BW", xlab="Residuals", col="lightblue", border="black")
par(mfrow=c(1,1))
dev.off()

cat("PMC/ANOVA analysis complete. Results in:", output_path, "\n")

# ---- COSINOR-2C ANALYSIS (AllTwin2Day/*.xlsx) ----

# Helper functions from Cosinor-2c.R
ConvertToHccDegrees <- function(acrophase) { if(acrophase > 0) acrophase - 360 else acrophase }
kittAcrophaseQuadrantCorrection <- function (cosCoeff, sinCoeff) {
  acrophase <- atan(abs(sinCoeff/cosCoeff))
  g <- 1
  if(cosCoeff * sinCoeff > 0) { g <- -1 }
  dphi <- 0
  if(cosCoeff < 0) { dphi <- -pi }
  acrophase <- dphi + g*acrophase
  return(acrophase)
}
ExtremaOf <- function(f, lowerT, upperT, method = "numeric", tol = .Machine$double.eps^0.5, ...) {
  df_fun <- function(x) numDeriv::grad(func = f, x = x)
  roots <- rootSolve::uniroot.all(f = df_fun, lower = lowerT, upper = upperT, tol = tol, ...)
  roots <- sort(unique(roots))
  roots <- roots[roots >= lowerT & roots <= upperT]
  vals <- vapply(roots, f, numeric(1))
  data.frame(t = roots, y = vals)
}
compute_two_component_cosinor <- function(local_data) {
  if (nrow(local_data) < 10) return(NULL)
  local_data$time_hours <- local_data$time * 24
  local_data$cos_24 <- cos(2 * pi * local_data$time_hours / 24)
  local_data$sin_24 <- sin(2 * pi * local_data$time_hours / 24)
  local_data$cos_12 <- cos(2 * pi * local_data$time_hours / 12)
  local_data$sin_12 <- sin(2 * pi * local_data$time_hours / 12)
  model <- lm(HR ~ cos_24 + sin_24 + cos_12 + sin_12, data = local_data)
  coefs <- coef(model)
  summary_model <- summary(model)
  A1 <- sqrt(coefs[2]^2 + coefs[3]^2)
  A2 <- sqrt(coefs[4]^2 + coefs[5]^2)
  acro24AsRadians <- kittAcrophaseQuadrantCorrection(coefs[2], coefs[3])
  if (acro24AsRadians < 0 ) acro24AsRadians <- acro24AsRadians + 2*pi
  acro24 <- acro24AsRadians * 180 / pi
  if (acro24 < 0) acro24 <- acro24 + 360
  acro12AsRadians <- kittAcrophaseQuadrantCorrection(coefs[4], coefs[5])
  if(acro12AsRadians < 0) acro12AsRadians <- acro12AsRadians + 2*pi
  acro12 <- acro12AsRadians * 180 / pi
  if (acro12 < 0) acro12 <- acro12 + 360
  percent_rhythm <- summary_model$r.squared * 100
  fstat <- summary_model$fstatistic
  pval <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
  MESOR <- coefs[1]
  cosinor2C <- function(t) {
    MESOR +
      A1 * cos((2 * pi / 24) * t + acro24AsRadians ) +
      A2 * cos((2 * pi / 12) * t + acro12AsRadians )
  }
  minT <- min(local_data$time_hours, na.rm=TRUE)
  range <- minT + 24
  extremaOf <- ExtremaOf(cosinor2C, minT, range, "numeric")
  Ortho_y <- max(extremaOf$y, na.rm=TRUE)
  Bathy_y <- min(extremaOf$y, na.rm=TRUE)
  Magnitude <- (Ortho_y - Bathy_y)/2
  imax <- which.max(extremaOf$y)
  Ortho_t_hour <- ((extremaOf[imax, ]$t))%%24
  Ortho_t <- Ortho_t_hour * -360/24
  imin <- which.min(extremaOf$y)
  Bathy_t_hour <- ((extremaOf[imin, ]$t))%%24
  Bathy_t <- Bathy_t_hour * -360/24
  data.frame(
    ID = unique(local_data$PID),
    N_observations = nrow(local_data),
    PercentRhythm = percent_rhythm,
    MESOR = MESOR,
    Amplitude_24h = A1,
    Acrophase_24h = ConvertToHccDegrees(acro24),
    Amplitude_12h = A2,
    Acrophase_12h = ConvertToHccDegrees(acro12),
    MAGNITUDE = Magnitude,
    Orthophase = ConvertToHccDegrees(Ortho_t),
    Bathyphase = ConvertToHccDegrees(Bathy_t),
    Pvalue = pval,
    F_statistic = fstat[1]
  )
}

# Process all xlsx files in AllTwin2Day
cosinor_files <- list.files(cosinor_input_folder, pattern = "\\.xlsx$", full.names = TRUE)
if (length(cosinor_files) == 0) stop("No xlsx files found in: ", cosinor_input_folder)

cosinor_results_list <- list()
for (f in cosinor_files) {
  cat("Processing cosinor file:", basename(f), "\n")
  d <- read_excel(f)
  # Expect columns: PID, time, HR (and possibly others)
  if (!all(c("PID","time","HR") %in% colnames(d))) {
    warning("File ", basename(f), " does not have required columns. Skipping.")
    next
  }
  # Split by PID if multiple PIDs per file
  for (this_pid in unique(d$PID)) {
    local_data <- d[d$PID == this_pid, ]
    cres <- compute_two_component_cosinor(local_data)
    if (!is.null(cres)) cosinor_results_list[[paste0(basename(f), "_", this_pid)]] <- cres
  }
}
cosinor_results <- do.call(rbind, cosinor_results_list)
write.csv(cosinor_results, file.path(output_path, "cosinor_results.csv"), row.names = FALSE)

cat("Cosinor-2c analysis complete. Results in:", output_path, "\n")

