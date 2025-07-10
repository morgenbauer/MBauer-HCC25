library(readxl)
library(this.path) # For determining the directory where the script is running

# Directory Setup
dir_to_script <- dirname(this.path::this.path())
input_path <- file.path(dir_to_script, "..", "GC")
if (!dir.exists(input_path)) stop("directory not found: ", input_path)
output_path <- file.path(dir_to_script, "..", "02-output")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

# Obtain the input file
xlsxFileToProcess <- list.files(
  path = input_path,
  pattern = "cosall.xlsx",
  full.names = TRUE
)
if (length(xlsxFileToProcess) == 0) stop("No cosall.xlsx file found in input_path!")

# Function to load data
process_sheet_data <- function(file_path, sheet_name) {
  read_excel(file_path, sheet = sheet_name)
}

input_data <- process_sheet_data(xlsxFileToProcess[1], "cosall")
if (is.null(input_data)) stop("Failed to load input_data!")

# Generalized PMC function
pmc_general <- function(data, PR_col, M_col, A_col, phi_col, phi_unit = "deg") {
  gdata <- data
  
  PR <- data[[PR_col]]
  M <- data[[M_col]]
  A <- data[[A_col]]
  phi <- data[[phi_col]]
  
  # Remove NA values
  valid_idx <- !is.na(PR) & !is.na(M) & !is.na(A) & !is.na(phi)
  PR <- PR[valid_idx]
  M <- M[valid_idx]
  A <- A[valid_idx]
  phi <- phi[valid_idx]
  n <- length(PR)
  
  phi_rad <- if (phi_unit == "deg") phi * pi / 180 else phi
  beta <- A * cos(phi_rad) 
  gamma <- A * -sin(phi_rad)
  
  gbeta <- beta
  ggamma <- gamma
  
  # Means
  mean_PR <- mean(PR, na.rm = TRUE)
  mean_M <- mean(M, na.rm = TRUE)
  mean_beta <- mean(beta, na.rm = TRUE)
  mean_gamma <- mean(gamma, na.rm = TRUE)
  mean_A <- sqrt(mean_beta^2 + mean_gamma^2)
  mean_phi_rad <- atan2(-mean_gamma, mean_beta)
  mean_phi_deg <- mean_phi_rad * 180 / pi
  if (mean_phi_deg > 0) { mean_phi_deg <- mean_phi_deg - 360 }
  
  SS_beta <- (sum((beta - mean_beta)^2)) / (n-1)
  SS_gamma <- (sum((gamma - mean_gamma)^2)) / (n-1)
  SS_beta_gamma <- (sum((beta - mean_beta) * (gamma - mean_gamma))) / (n-1)
  matrix2x2 <- matrix(c(SS_beta, SS_beta_gamma, SS_beta_gamma, SS_gamma), nrow=2, ncol=2, byrow=TRUE)
  invert2x2 <- solve(matrix2x2)
  
  #See kits equation 58 and 59
  # Compute variance and covariance matrix
  c22 <- (SS_beta*mean_beta^2 + 2*(SS_beta_gamma)*(mean_beta)*(mean_gamma) + SS_gamma * mean_gamma^2) / (n*(mean_A^2))
  c23 <- (-(SS_beta - SS_gamma) * mean_beta * mean_gamma + SS_beta_gamma * (mean_beta^2 - mean_gamma^2)) / (n*(mean_A^2))
  c33 <- (SS_beta * mean_gamma^2 - 2*(SS_beta_gamma)*(mean_beta)*(mean_gamma) + SS_gamma *mean_beta^2) / (n*(mean_A^2))
  
  #Kit's formula 61, for variable r (named ratio here)
  ratio <- SS_beta_gamma / (sqrt(SS_beta * SS_gamma))
  # Cell Q113 in PMCr excel (intermediate equation)
  GCparen <- (mean_beta^2 / SS_beta) - 2*ratio*mean_beta*mean_gamma / sqrt(SS_beta * SS_gamma) + (mean_gamma^2 / SS_gamma)
  # Cell Q114 
  GCmult <- (n*(n-2))/((2*(n-1))*(1-ratio^2))
  # F Value from Excel sheet (Q115)
  fval <- GCmult*GCparen
  #Now we have a the value of the f distribution
  # P Value for overall model, 2 df in numerator (elipse for A and phi)
  pval <- df(fval,2,n-2)

    # Confidence intervals for amplitude and acrophase
  den <- mean_A^2
  num1 <- (mean_beta^2 * c22) + (2 * mean_beta * mean_gamma * c23) + (mean_gamma^2 * c33)
  num2 <- (mean_gamma^2 * c22) - (2 * mean_beta * mean_gamma * c23) + (mean_beta^2 * c33)
  num3 <- if (den > 0) num2 / den^2 else NA
  sqrt_num3 <- if (!is.na(num3) && num3 > 0) sqrt(num3) else NA
  
  # 95% confidence intervals
  t_crit <- qt(0.975, df = n - 1)
  F_crit <- qf(0.95, df1 = 2, df2 = n - 2)
  
  # Amplitude CI
  CI_A <- if (!is.na(sqrt_num3)) sqrt(2 * F_crit * sqrt_num3) else NA
  A_CI_lower <- mean_A - CI_A
  A_CI_upper <- mean_A + CI_A
  
  # Acrophase CI
  SE_phi_rad <- if (mean_A > 0) sqrt(num1) / mean_A else NA
  SE_phi_deg <- if (!is.na(SE_phi_rad)) SE_phi_rad * 180 / pi else NA
  CI_phi_deg <- if (!is.na(SE_phi_deg)) t_crit * SE_phi_deg else NA
  phi_CI_lower <- mean_phi_deg - CI_phi_deg
  phi_CI_upper <- mean_phi_deg + CI_phi_deg
  
  # Output all three unique entries for 2x2 matrix and its inverse
  list(
    n = n,
    mean_PR = mean_PR,
    mean_M = mean_M,
    mean_A = mean_A,
    mean_phi_deg = mean_phi_deg,
    pval = pval,
    # Amplitude confidence intervals
    A_CI_lower = A_CI_lower,
    A_CI_upper = A_CI_upper,
    # Acrophase confidence intervals
    SE_phi_deg = SE_phi_deg,
    phi_CI_lower = phi_CI_lower,
    phi_CI_upper = phi_CI_upper,
    # Matrix entries
    SS_beta = SS_beta,
    SS_beta_gamma = SS_beta_gamma,
    SS_gamma = SS_gamma,
    # Variance-covariance entries
    c22 = c22,
    c23 = c23,
    c33 = c33,
    # 2x2 matrix and its inverse (all entries)
    matrix2x2_11 = matrix2x2[1,1],
    matrix2x2_12 = matrix2x2[1,2],
    matrix2x2_21 = matrix2x2[2,1],
    matrix2x2_22 = matrix2x2[2,2],
    invert2x2_11 = invert2x2[1,1],
    invert2x2_12 = invert2x2[1,2],
    invert2x2_21 = invert2x2[2,1],
    invert2x2_22 = invert2x2[2,2]
  )
  }



# Define analyses: list of lists with column names and a label
# Changed column names in excel sheet PR(Tau1) changed to PRTau1. Otherwise RStudio, interpretes it as a function
# Fix me: Chase wonders what is pmc_analysis, is it about column names from spreadsheet or something else?
pmc_analyses <- list(
  PMC_24h = list(PR = "PRTau1", MESOR = "MESOR", A = "A1", Phi = "Phi1", label = "24h"),
  PMC_12h = list(PR = "PRTau2", MESOR = "MESOR", A = "A2", Phi = "Phi2", label = "12h"),
  PMC_ortho = list(PR = "PR", MESOR = "MESOR", A = "Aoverall", Phi = "ϕOrthophase", label = "Composite-Orthophase"),
  PMC_bathy = list(PR = "PR", MESOR = "MESOR", A = "Aoverall", Phi = "ϕBathyphase", label = "Composite-Bathyphase")
)

# Run for entire population
results_population <- list()
for (analysis in names(pmc_analyses)) {
  cols <- pmc_analyses[[analysis]]
  res <- pmc_general(
    data = input_data,
    PR_col = cols$PR,
    M_col = cols$MESOR,
    A_col = cols$A,
    phi_col = cols$Phi,
    phi_unit = "deg"
  )
  results_population[[analysis]] <- res
  cat(sprintf("\nPopulation Mean Cosinor (%s):\n", cols$label))
  print(res)
}

# Run for each BW group (col BWgrp)
results_by_bw <- list()
if ("BWgrp" %in% names(input_data)) {
  bw_groups <- unique(input_data$BWgrp)
  for (bw in bw_groups) {
    cat(sprintf("\n==== BW Group: %s ====\n", bw))
    group_data <- input_data[input_data$BWgrp == bw, ]
    group_results <- list()
    for (analysis in names(pmc_analyses)) {
      cols <- pmc_analyses[[analysis]]
      res <- pmc_general(
        data = group_data,
        PR_col = cols$PR,
        M_col = cols$MESOR,
        A_col = cols$A,
        phi_col = cols$Phi,
        phi_unit = "deg"
      )
      group_results[[analysis]] <- res
      cat(sprintf("\nBW Group %s - PMC (%s):\n", bw, cols$label))
      print(res)
    }
    results_by_bw[[as.character(bw)]] <- group_results
  }
}

# --- Save ALL BW group results in ONE CSV ---
if (length(results_by_bw) > 0) {
  # Create a single consolidated data frame for all BW groups
  all_bw_df <- do.call(rbind, lapply(names(results_by_bw), function(bw) {
    # For each BW group, create rows for each analysis
    bw_group_df <- do.call(rbind, lapply(names(results_by_bw[[bw]]), function(analysis) {
      cbind(
        BW_group = bw,
        analysis = analysis,
        as.data.frame(results_by_bw[[bw]][[analysis]])
      )
    }))
    return(bw_group_df)
  }))
  
  # Write the consolidated BW results to a single CSV
  write.csv(all_bw_df, 
            file = file.path(output_path, "all_BWgroups_pmc_results.csv"), 
            row.names = FALSE)
  
  cat("All BW group results saved to: all_BWgroups_pmc_results.csv\n")
}

# --- Save population results as CSV ---
# Convert results_population (a list) to a data frame
population_df <- do.call(rbind, lapply(names(results_population), function(name) {
  cbind(analysis = name, as.data.frame(results_population[[name]]))
}))
write.csv(population_df, file = file.path(output_path, "population_pmc_results.csv"), row.names = FALSE)


