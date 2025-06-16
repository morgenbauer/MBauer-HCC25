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
library(readxl)  # For reading Excel files

# Function to perform two-component cosinor analysis
compute_two_component_cosinor <- function(local_data) {
  # Ensure we have enough data points
  if(nrow(local_data) < 10) {
    warning("Insufficient data points for two-component cosinor analysis")
    return(NULL)
  }
  
  tryCatch({
    # Convert time from days to hours for the analysis
    #TODO: might be incorrect
    local_data$time_hours <- local_data$time * 24
    
    # Create the two-component cosinor model
    # 24-hour component (circadian)
    local_data$cos_24 <- cos(2 * pi * local_data$time_hours / 24)
    local_data$sin_24 <- sin(2 * pi * local_data$time_hours / 24)
    
    # 12-hour component (ultradian)  
    local_data$cos_12 <- cos(2 * pi * local_data$time_hours / 12)
    local_data$sin_12 <- sin(2 * pi * local_data$time_hours / 12)
    
    # Fit the two-component model
    model <- lm(HR ~ cos_24 + sin_24 + cos_12 + sin_12, data = local_data)
    
    # Extract coefficients
    coefs <- coef(model)
    M <- coefs[1]  # MESOR (mean level)
    
    # 24-hour component parameters
    coef_cos <- coefs[2]
    coef_sin <- coefs[3]
    A1 <- sqrt(coef_cos^2 + coef_sin^2)  # Amplitude for 24h component
    
    # Calculate acrophase for 24-hour component (in hours)
    acrophase_24_rad <- atan2(coef_sin, coef_cos)
    acrophase_24_hours <- -acrophase_24_rad * 24 / (2 * pi)
    if(acrophase_24_hours < 0) acrophase_24_hours <- acrophase_24_hours + 24
    
    # 12-hour component parameters
    coef_cos_12 <- coefs[4]
    coef_sin_12 <- coefs[5]
    A2 <- sqrt(coef_cos_12^2 + coef_sin_12^2)  # Amplitude for 12h component
    
    # Calculate acrophase for 12-hour component (in hours)
    acrophase_12_rad <- atan2(coef_sin_12, coef_cos_12)
    acrophase_12_hours <- -acrophase_12_rad * 12 / (2 * pi)
    if(acrophase_12_hours < 0) acrophase_12_hours <- acrophase_12_hours + 12
    
    # Compute model statistics
    model_summary <- summary(model)
    
    # Percent rhythm calculation (from R-squared)
    percent_rhythm <- model_summary$r.squared * 100
    
    # P-value calculation from overall F-statistic of the model
    f_value <- model_summary$fstatistic
    p_value <- pf(f_value[1], f_value[2], f_value[3], lower.tail = FALSE)
    
    return(list(
      Sheet_ID = as.character(local_data$Sheet_ID[1]),
      PercentRhythm = percent_rhythm,
      MESOR = M,
      Amplitude_24h = A1,
      Acrophase_24h_hours = acrophase_24_hours,
      Amplitude_12h = A2,
      Acrophase_12h_hours = acrophase_12_hours,
      Pvalue = p_value,
      F_statistic = f_value[1],
      N_observations = nrow(local_data)
    ))
  }, error = function(e) {
    warning(paste("Error in two-component cosinor analysis:", e$message))
    return(NULL)
  })
}

# Function to process Excel sheet data for two-component analysis
process_sheet_data <- function(file_path, sheet_name) {
  # Add debug flag
  verbose <- TRUE
  
  if(verbose) cat("\nProcessing sheet", sheet_name, "from", basename(file_path), "\n")
  
  # Read the Excel sheet
  data <- tryCatch({
    readxl::read_excel(file_path, sheet = sheet_name)
  }, error = function(e) {
    warning(paste("Error reading sheet", sheet_name, "from file", file_path, ":", e$message))
    return(NULL)
  })
  
  # Check if data was loaded successfully
  if(is.null(data) || nrow(data) < 10) {
    warning(paste("Insufficient data in sheet:", sheet_name))
    return(NULL)
  }
  
  if(verbose) cat("Raw data rows:", nrow(data), "columns:", ncol(data), "\n")
  
  # Extract columns C and D (3rd and 4th columns)
  if(ncol(data) < 4) {
    warning(paste("Sheet", sheet_name, "does not have enough columns (need at least 4)"))
    return(NULL)
  }
  
  # Create simplified dataset with columns C (time) and D (HR)
  simplified_data <- data.frame(
    Sheet_ID = sheet_name,
    time = as.numeric(data[[3]]),  # Column C
    HR = as.numeric(data[[4]])     # Column D
  )
  
  if(verbose) cat("Extracted time and HR data from columns C and D\n")
  
  # Remove rows with NAs in critical columns
  orig_rows <- nrow(simplified_data)
  simplified_data <- simplified_data[!is.na(simplified_data$time) & !is.na(simplified_data$HR), ]
  if(verbose) cat("Removed", orig_rows - nrow(simplified_data), "rows with NAs. Remaining:", nrow(simplified_data), "\n")
  
  # Check if we still have enough data
  if(nrow(simplified_data) < 10) {
    warning(paste("Insufficient valid data points in sheet:", sheet_name))
    return(NULL)
  }
  
  # Fit two-component cosinor model
  if(verbose) cat("Running two-component cosinor analysis with", nrow(simplified_data), "data points\n")
  cosinor_results <- compute_two_component_cosinor(simplified_data)
  
  # If results are NULL, return NULL
  if(is.null(cosinor_results)) {
    if(verbose) cat("Two-component cosinor analysis failed for sheet", sheet_name, "\n")
    return(NULL)
  } else {
    if(verbose) cat("Two-component cosinor analysis successful for sheet", sheet_name, "\n")
  }
  
  return(cosinor_results)
}

# Set the directory path where this script is running
dir_to_RStudio_script <- getwd()

# Set paths relative to project root
input_path <- file.path(dir_to_RStudio_script, "..", "GC")
output_path <- file.path(dir_to_RStudio_script, "..", "02-output")

# Create output directory if it doesn't exist
if(!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
  cat("Created output directory:", output_path, "\n")
}

cat("RStudio Script:", dir_to_RStudio_script, "\n")
cat("Input Path:", input_path, "\n")
cat("Output Path:", output_path, "\n")

# Verify the GC directory exists
if(!dir.exists(input_path)) {
  stop(paste("GC directory not found at:", input_path))
}

# Look for Excel files in the input directory
excel_files <- list.files(input_path, pattern = "\\.xlsx$", full.names = TRUE)
#Now we filter out temporary excel files
excel_files <- excel_files[!grepl("^~", basename(excel_files))]

cat("Found", length(excel_files), "Excel files to process\n")

# Print the Excel files found for debugging
if(length(excel_files) > 0) {
  cat("Excel files found:\n")
  for(file in excel_files) {
    cat(" -", basename(file), "\n")
  }
} else {
  cat("No Excel files found in:", input_path, "\n")
  cat("Please verify that .xlsx files exist in the GC directory\n")
}

# Check if any Excel files were found
if(length(excel_files) == 0) {
  stop("No Excel files found in input directory. Please check the path.")
}

# Process the specific sheets HR051a and HR051b from Excel file
all_results <- list()
sheets_to_process <- c("HR051a", "HR051b")

for(excel_file in excel_files) {
  cat("\nProcessing Excel file:", basename(excel_file), "\n")
  
  for(sheet_name in sheets_to_process) {
    cat("Processing sheet:", sheet_name, "\n")
    
    result <- process_sheet_data(excel_file, sheet_name)
    
    if(!is.null(result)) {
      # Create a unique identifier for this result
      result_id <- paste(basename(excel_file), sheet_name, sep = "_")
      all_results[[result_id]] <- as.data.frame(result)
      cat("Successfully processed:", result_id, "\n")
    } else {
      cat("Failed to process sheet:", sheet_name, "from file:", basename(excel_file), "\n")
    }
  }
}

# Print results summary
cat("\nResults summary:\n")
cat("Total results collected:", length(all_results), "\n")

# Check if any results exist
if(length(all_results) > 0) {
  # Convert to a dataframe
  results_df <- do.call(rbind, all_results)
  
  # Reorder columns for clarity
  results_df <- results_df[, c("Sheet_ID", "N_observations", "PercentRhythm", "MESOR", 
                               "Amplitude_24h", "Acrophase_24h_hours", 
                               "Amplitude_12h", "Acrophase_12h_hours", 
                               "Pvalue", "F_statistic")]
  
  # Write results to CSV
  output_file <- file.path(output_path, "two_component_cosinor_results.csv")
  write.csv(results_df, output_file, row.names = FALSE)
  
  # Print success message
  cat("\nSuccessfully created", output_file, "with two-component cosinor analysis results\n")
  cat("Results contain data for", nrow(results_df), "sheet analyses\n")
  
  # Print summary of results
  cat("\nSummary of Results:\n")
  cat("==================\n")
  print(results_df)
  
  cat("\nInterpretation Guide:\n")
  cat("- Sheet_ID: Which sheet was analyzed (HR051a or HR051b)\n")
  cat("- N_observations: Number of data points used in analysis\n")
  cat("- PercentRhythm: Percentage of variance explained by the model\n")
  cat("- MESOR: Mean heart rate level\n")
  cat("- Amplitude_24h: Magnitude of 24-hour (circadian) rhythm\n")
  cat("- Acrophase_24h_hours: Peak time of 24-hour rhythm (hours from midnight)\n")
  cat("- Amplitude_12h: Magnitude of 12-hour (ultradian) rhythm\n")
  cat("- Acrophase_12h_hours: Peak time of 12-hour rhythm (hours within 12h cycle)\n")
  cat("- Pvalue: Statistical significance of the overall model\n")
  cat("- F_statistic: F-statistic value for model significance\n")
  
} else {
  cat("\nERROR: No valid results found. Please check your input files.\n")
  cat("Possible issues:\n")
  cat("1. Excel files don't contain sheets named 'HR051a' or 'HR051b'\n")
  cat("2. Columns C and D don't contain valid time and HR data\n")
  cat("3. Not enough valid data points after removing NA values (need at least 10)\n")
  cat("4. Two-component cosinor model failed to fit\n")
}

# Print final diagnostic information
cat("\nFinal Summary:\n")
cat("=============\n")
cat("Excel files found:", length(excel_files), "\n")
cat("Sheets processed per file:", length(sheets_to_process), "\n")
cat("Expected total analyses:", length(excel_files) * length(sheets_to_process), "\n")
cat("Successful analyses:", length(all_results), "\n")
cat("Failed analyses:", (length(excel_files) * length(sheets_to_process)) - length(all_results), "\n")

