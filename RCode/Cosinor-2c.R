#2c Cosinor
#Purpose here is to create an R studio verson for two component cosinor model
#Goal 1: Create path to read in GC spreadsheet
#Goal 2: Find RStudio Libraries that know how to import .xlsx for specific tabs in them
#Goal 3: Process imported data to remove or ignore columms E and H (HR051a)
#Goal 4: Process data for 2c

# Two-Component Cosinor Analysis for Heart Rate Data
# HR(t) = M + A1*cos(2π(t-t0)/24 + φ1) + A2*cos(2π(t-t0)/12 + φ2) + e(t)

library(dplyr)
library(ggplot2)
library(readxl)
# libs for ExtremaOf
library(Deriv)
library(Ryacas)
library(numDeriv)
library(rootSolve)

#ConvertToHccDegrees ensures results are consistent with HCC philosophy that
#the zero degree is at the 12 o'clock positon.
ConvertToHccDegrees <- function(acrophase){
  if(acrophase >0) {acrophase - 360} else{acrophase}
}

# ExtremaOf returns the list of x,y pairs where the slope is zero.
# From this list, select the MinMax to return the Orthophase and Bathyphase
# ExtremaOf(f, lowerT, upperT, "numeric")
ExtremaOf <- function(f,
                      lowerT, upperT,
                      method = c("symbolic", "numeric"),
                      tol    = .Machine$double.eps^0.5,
                      ...) 
{
  method <- match.arg(method)
  
  if (method == "symbolic") {
    # --- Symbolic: f must be given as an *expression* in `t`, e.g. quote(sin(t)+t^2)
    if (!requireNamespace("Deriv", quietly = TRUE) ||
        !requireNamespace("Ryacas", quietly = TRUE)) {
      stop("For symbolic method, install Deriv and Ryacas")
    }
    f_expr  <- substitute(f)
    df_expr <- Deriv::Deriv(f_expr, "t")
    
    # Build a Yacas solve call
    sol_txt <- Ryacas::y_fn("Solve",
                            Ryacas::ysym(df_expr) == 0,
                            "t")
    # parse the returned list of rules
    raw    <- Ryacas::yac_str(sol_txt)
    # e.g. raw == "{ {t -> 1/2}, {t -> 3/2} }"
    # strip braces & split
    nums   <- gsub("^\\{+|\\}+","", raw)
    parts  <- unlist(strsplit(nums, "\\},\\{", perl=TRUE))
    roots  <- vapply(parts, function(p) {
      # each p like "t -> 1/2"
      as.numeric(eval(parse(text=sub(".*->","", p))))
    }, numeric(1))
    
  } else {
    # --- Numeric: f is a numeric function f(t)
    if (!requireNamespace("numDeriv", quietly = TRUE) ||
        !requireNamespace("rootSolve", quietly = TRUE)) {
      stop("For numeric method, install numDeriv and rootSolve")
    }
    # numeric derivative
    df_fun <- function(x) numDeriv::grad(func = f, x = x)
    # find *all* zeros of df_fun in [lowerT,upperT]
    roots <- rootSolve::uniroot.all(f = df_fun,
                                    lower = lowerT,
                                    upper = upperT,
                                    tol   = tol,
                                    ...)
  }
  
  # Filter and evaluate f at each root
  roots <- sort(unique(roots))
  roots <- roots[roots >= lowerT & roots <= upperT]
  
  # Return a two‐column data.frame: t and f(t)
  vals <- vapply(roots, f, numeric(1))
  data.frame(t = roots, y = vals)
}


# Function to perform cosinor analysis (full + individual components)
compute_two_component_cosinor <- function(local_data) {
  if (nrow(local_data) < 10) {
    warning("Insufficient data points for cosinor analysis")
    return(NULL)
  }
  
  tryCatch({
    local_data$time_hours <- local_data$time * 24
    local_data$cos_24 <- cos(2 * pi * local_data$time_hours / 24)
    local_data$sin_24 <- sin(2 * pi * local_data$time_hours / 24)
    local_data$cos_12 <- cos(2 * pi * local_data$time_hours / 12)
    local_data$sin_12 <- sin(2 * pi * local_data$time_hours / 12)
    
    rows <- list()
    
    # Full model
    model <- lm(HR ~ cos_24 + sin_24 + cos_12 + sin_12, data = local_data)
    coefs <- coef(model)
    summary_model <- summary(model)
    A1 <- sqrt(coefs[2]^2 + coefs[3]^2)
    A2 <- sqrt(coefs[4]^2 + coefs[5]^2)
    acro24 <- -atan2(coefs[3], coefs[2]) * 180 / pi
    acro12 <- -atan2(coefs[5], coefs[4]) * 180 / pi
    if (acro24 < 0) acro24 <- acro24 + 360
    if (acro12 < 0) acro12 <- acro12 + 360
    percent_rhythm <- summary_model$r.squared * 100
    fstat <- summary_model$fstatistic
    pval <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
    MESOR <- coefs[1]
    
    # BEGIN COMPUTE Orthophase and Bathyphase
    # Create 2-component cosinor function
    # NOTE: t must be expressed in hours
    cosinor2C <- function(t) {
      MESOR +
        A1 * cos((2 * pi / 24) * t * acro24) +
        A2 * cos((2 * pi / 12) * t * acro12)
    }
    # specify the time range to use
    minT <- min(local_data$time_hours, na.rm=TRUE)
    maxT <- max(local_data$time_hours, na.rm=TRUE)
    # Ensure one complete cycle using the largest period of the harmonic
    range <- minT + 24
    # TODO : If tau components are not harmonic, return NA

    # Compute where zero slope occurs
    extremaOf <- ExtremaOf(cosinor2C, minT, range, "numeric")
 
    # Compute the MAGNITUDE by subtracting Amplitude of Orthophase and Bathyphase
    Ortho_y  <- max(extremaOf$y, na.rm=TRUE)
    Bathy_y  <- min(extremaOf$y, na.rm=TRUE)
    Magnitude <- Ortho_y - Bathy_y
    
    # Determine what times are associated with thos maximums
    imax  <- which.max(extremaOf$y)
    Ortho_t <- (extremaOf[imax, ]$t / 24) * 180 / pi
    if (Ortho_t < 0) Ortho_t <- Ortho_t + 360
    
    # Which row has the smallest y?
    imin  <- which.min(extremaOf$y)
    Bathy_t <- ( extremaOf[imin, ]$t / 24) * 180 / pi
    if (Bathy_t < 0) Bathy_t <- Bathy_t + 360
    
    
    # END COMPUTE Orthophase and Bathyphase
    
    rows[["Full"]] <- data.frame(
      Model_Type = "Full",
      N_observations = nrow(local_data),
      PercentRhythm = percent_rhythm,
      MESOR = MESOR,
      Amplitude = NA,
      Acrophase_deg = NA,
      Pvalue = pval,
      F_statistic = fstat[1],
      Amplitude_24h_full = A1,
      Acrophase_24h_full = ConvertToHccDegrees(acro24),
      Amplitude_12h_full = A2,
      Acrophase_12h_full = ConvertToHccDegrees(acro12),
      MAGNITUDE = Magnitude,
      Orthophase = ConvertToHccDegrees(Ortho_t),
      Bathyphase = ConvertToHccDegrees(Bathy_t)
    )
    
    # 24h only
    model_24 <- lm(HR ~ cos_24 + sin_24, data = local_data)
    coefs24 <- coef(model_24)
    A24 <- sqrt(coefs24[2]^2 + coefs24[3]^2)
    acro24_indiv <- -atan2(coefs24[3], coefs24[2]) * 180 / pi
    if (acro24_indiv < 0) acro24_indiv <- acro24_indiv + 360
    summary24 <- summary(model_24)
    PR24 <- summary24$r.squared * 100
    f24 <- summary24$fstatistic
    p24 <- pf(f24[1], f24[2], f24[3], lower.tail = FALSE)
    
    rows[["24h"]] <- data.frame(
      Model_Type = "24h_only",
      N_observations = nrow(local_data),
      PercentRhythm = PR24,
      MESOR = coefs24[1],
      Amplitude = A24,
      Acrophase_deg = ConvertToHccDegrees(acro24_indiv),
      Pvalue = p24,
      F_statistic = f24[1],
      Amplitude_24h_full = NA,
      Acrophase_24h_full = NA,
      Amplitude_12h_full = NA,
      Acrophase_12h_full = NA,
      MAGNITUDE = NA,
      Orthophase = NA,
      Bathyphase = NA
    )
    
    # 12h only
    model_12 <- lm(HR ~ cos_12 + sin_12, data = local_data)
    coefs12 <- coef(model_12)
    A12 <- sqrt(coefs12[2]^2 + coefs12[3]^2)
    acro12_indiv <- -atan2(coefs12[3], coefs12[2]) * 180 / pi
    if (acro12_indiv < 0) acro12_indiv <- acro12_indiv + 360
    summary12 <- summary(model_12)
    PR12 <- summary12$r.squared * 100
    f12 <- summary12$fstatistic
    p12 <- pf(f12[1], f12[2], f12[3], lower.tail = FALSE)
    
    rows[["12h"]] <- data.frame(
      Model_Type = "12h_only",
      N_observations = nrow(local_data),
      PercentRhythm = PR12,
      MESOR = coefs12[1],
      Amplitude = A12,
      Acrophase_deg = ConvertToHccDegrees(acro12_indiv),
      Pvalue = p12,
      F_statistic = f12[1],
      Amplitude_24h_full = NA,
      Acrophase_24h_full = NA,
      Amplitude_12h_full = NA,
      Acrophase_12h_full = NA,
      MAGNITUDE = NA,
      Orthophase = NA,
      Bathyphase = NA
    )
    
    # Add sheet ID to all rows
    for (i in names(rows)) {
      rows[[i]]$Sheet_ID <- local_data$Sheet_ID[1]
    }
    
    # Bind all model results together
    return(do.call(rbind, rows))
    
  }, error = function(e) {
    warning(paste("Cosinor analysis error:", e$message))
    return(NULL)
  })
}

# Process a sheet
process_sheet_data <- function(file_path, sheet_name) {
  cat("\nProcessing sheet:", sheet_name, "\n")
  data <- tryCatch({
    read_excel(file_path, sheet = sheet_name)
  }, error = function(e) {
    warning(paste("Error reading sheet:", e$message))
    return(NULL)
  })
  
  if (is.null(data) || nrow(data) < 10) return(NULL)
  
  if (sheet_name == "HR051a") {
    time_col <- 3
    hr_col <- 4
  } else if (sheet_name == "HR051b") {
    time_col <- 3
    hr_col <- 4
  } else {
    warning("Unknown sheet format")
    return(NULL)
  }
  
  simplified_data <- data.frame(
    Sheet_ID = sheet_name,
    time = as.numeric(data[[time_col]]),
    HR = as.numeric(data[[hr_col]])
  )
  
  simplified_data <- simplified_data[!is.na(simplified_data$time) & !is.na(simplified_data$HR), ]
  if (nrow(simplified_data) < 10) return(NULL)
  
  compute_two_component_cosinor(simplified_data)
}

# Directory Setup
dir_to_script <- getwd()
input_path <- file.path(dir_to_script, "..", "GC")
output_path <- file.path(dir_to_script, "..", "02-output")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

target_file <- "MB-Tw051-2c.xlsx"
excel_file <- file.path(input_path, target_file)
if (!file.exists(excel_file)) stop("File not found: ", excel_file)

sheets_to_process <- c("HR051a", "HR051b")
all_results <- list()

for (sheet_name in sheets_to_process) {
  result <- process_sheet_data(excel_file, sheet_name)
  if (!is.null(result)) {
    all_results[[sheet_name]] <- result
    cat("✔ Processed:", sheet_name, "\n")
  } else {
    cat("✖ Failed:", sheet_name, "\n")
  }
}

# Output in long format
if (length(all_results) > 0) {
  results_df <- do.call(rbind, all_results)
  results_df <- results_df[, c(
    "Sheet_ID", "Model_Type", "N_observations", "PercentRhythm", "MESOR",
    "Amplitude", "Acrophase_deg",
    "Amplitude_24h_full", "Acrophase_24h_full",
    "Amplitude_12h_full", "Acrophase_12h_full",
    "Pvalue", "F_statistic", "MAGNITUDE", "Orthophase", "Bathyphase"
  )]
  
  output_file <- file.path(output_path, "2c_cosionor_results.csv")
  write.csv(results_df, output_file, row.names = FALSE)
  cat("\n Output saved to:", output_file, "\n")
  print(results_df)
} else {
  cat("\n No valid results generated.\n")
}


