
library(readxl)

# Directory Setup
dir_to_script <- getwd()
input_path <- file.path(dir_to_script, "..", "GC")
if (!dir.exists(input_path)) stop("directory not found: ", input_path)
output_path <- file.path(dir_to_script, "..", "02-output")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

# Obtain the input file
xlsxFileToProcess <- list.files(
  path       = input_path, 
  pattern    = "cosall.xlsx", 
  full.names = TRUE
)
if (length(xlsxFileToProcess) == 0) stop("No cosall.xlsx file found in input_path!")

# Function to load data
process_sheet_data <- function(file_path, sheet_name) {
  cat("Processing sheet:", sheet_name, "\n")
  cat("Processing path:", file_path, "\n")
  data <- tryCatch({
    read_excel(file_path, sheet = sheet_name)
  }, error = function(e) {
    warning(paste("Error reading sheet:", e$message))
    return(NULL)
  })
  return(data)
}

# Load the data
input_data <- process_sheet_data(xlsxFileToProcess[1], "cosall")
if (is.null(input_data)) stop("Failed to load input_data!")

# Generalized PMC function
pmc_general <- function(data, PR_col, M_col, A_col, phi_col, phi_unit = "deg") {
  PR <- data[[PR_col]]
  M <- data[[M_col]]
  A <- data[[A_col]]
  phi <- data[[phi_col]]
  phi_rad <- if (phi_unit == "deg") phi * pi / 180 else phi
  beta <- A * cos(phi_rad)
  gamma <- A * sin(phi_rad)
  mean_PR <- mean(PR, na.rm = TRUE)
  mean_M <- mean(M, na.rm = TRUE)
  mean_beta <<- mean(beta, na.rm = TRUE)
  mean_gamma <<- mean(gamma, na.rm = TRUE)
  mean_A <- sqrt(mean_beta^2 + mean_gamma^2)
  mean_phi_rad <- atan2(mean_gamma, mean_beta)
  mean_phi_deg <- mean_phi_rad * 180 / pi
  # Ensure acrophases are always reported as negative values
  if (mean_phi_deg > 0) { mean_phi_deg <- mean_phi_deg - 360}
  
  beta_sq <<- sum((beta - mean_beta)^2)
  gamma_sq <<- sum((gamma - mean_gamma)^2)
  beta_gamma_crs <<- sum((beta - mean_beta) * (gamma - mean_gamma))
  matrix2x2 <<- matrix(c(beta_sq, beta_gamma_crs, beta_gamma_crs, gamma_sq), nrow=2, ncol=2, byrow=TRUE)
  invert2x2 <<- solve(matrix2x2)
  
  list(
    mean_PR = mean_PR,
    mean_M = mean_M,
    mean_A = mean_A,
    mean_phi_deg = mean_phi_deg
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

# --- Save population results as CSV ---
# Convert results_population (a list) to a data frame
population_df <- do.call(rbind, lapply(names(results_population), function(name) {
  cbind(analysis = name, as.data.frame(results_population[[name]]))
}))
write.csv(population_df, file = file.path(output_path, "population_pmc_results.csv"), row.names = FALSE)

# --- Save BW group results as CSV ---
# For each BW group, create a CSV
if (length(results_by_bw) > 0) {
  for (bw in names(results_by_bw)) {
    bw_df <- do.call(rbind, lapply(names(results_by_bw[[bw]]), function(name) {
      cbind(analysis = name, as.data.frame(results_by_bw[[bw]][[name]]))
    }))
    write.csv(bw_df,
              file = file.path(output_path, paste0("BWgroup_", bw, "_pmc_results.csv")),
              row.names = FALSE)
  }
}




