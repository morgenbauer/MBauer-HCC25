library(dplyr)
library(ggplot2)
library(readxl)  # For reading Excel files
library(this.path) # For determining the directory where the script is running

# Function to compute cosinor parameters with the correct acrophase calculation
compute_cosinor <- function(local_data) {
  # Ensure we have enough data points
  if(nrow(local_data) < 3) {
    warning("Insufficient data points for cosinor analysis")
    return(NULL)
  }
  
  # Fit the model
  tryCatch({
    model <- lm(HR ~ cos(2*pi*time) + sin(2*pi*time), data = local_data)
    
    # Extracting coefficients
    M <- coef(model)[1]  # MESOR
    beta_cos <- coef(model)[2]
    beta_sin <- coef(model)[3]
    
    # Calculating amplitude
    A <- sqrt(beta_cos^2 + beta_sin^2)
    
    # CORRECT ACROPHASE CALCULATION
    acrophase <- 0
    acrophase <- atan(abs(beta_sin/beta_cos))
    g <- 1
    if(beta_cos*beta_sin > 0) {
      g <- -1
    }
    dphi <- 0
    if(beta_cos < 0) {
      dphi <- -pi
    }
    acrophase <- dphi + g*acrophase
    acrophase_degrees <- acrophase * 180/pi
    
    # Compute model statistics
    model_summary <- summary(model)
    
    # Percent rhythm calculation (from R-squared)
    percent_rhythm <- model_summary$r.squared * 100
    
    # P-value calculation from overall F-statistic of the model
    f_value <- model_summary$fstatistic
    p_value <- pf(f_value[1], f_value[2], f_value[3], lower.tail = FALSE)
    
    return(list(
      Twin_ID = as.character(local_data$Twin_ID[1]),
      PercentRhythm = percent_rhythm,
      MESOR = M, 
      Amplitude = A, 
      Acrophase_degrees = acrophase_degrees,
      Pvalue = p_value,
      F_statistic = f_value[1]
    ))
  }, error = function(e) {
    warning(paste("Error in cosinor analysis:", e$message))
    return(NULL)
  })
}

# Function to convert Excel files to CSV
convert_excel_to_csv <- function(dir_path) {
  # List all Excel files in the directory
  excel_files <- list.files(dir_path, pattern = "\\.xlsx$", full.names = TRUE)
  
  cat("Found", length(excel_files), "Excel files to convert\n")
  
  # Create a converted_files vector to track successful conversions
  converted_files <- character(0)
  
  # Process each Excel file
  for(excel_file in excel_files) {
    # Extract the base filename without extension
    base_name <- sub("\\.xlsx$", "", basename(excel_file))
    
    # Define the output CSV file path
    csv_file <- file.path(dir_path, paste0(base_name, ".csv"))
    
    cat("Converting:", basename(excel_file), "to CSV\n")
    
    # Try to read the Excel file
    tryCatch({
      # Read the Excel file
      excel_data <- readxl::read_excel(excel_file)
      
      # Write to CSV
      write.csv(excel_data, file = csv_file, row.names = FALSE)
      
      # Add to the list of successfully converted files
      converted_files <- c(converted_files, csv_file)
      
      cat("Successfully converted:", basename(excel_file), "to", basename(csv_file), "\n")
    }, error = function(e) {
      warning(paste("Error converting", basename(excel_file), ":", e$message))
    })
  }
  
  cat("Converted", length(converted_files), "Excel files to CSV format\n")
  return(converted_files)
}

# Function to process twin data
process_twin_data <- function(file_path) {
  # Add debug flag
  verbose <- TRUE
  
  # Extract twin ID from filename
  twin_id <- basename(file_path)
  twin_id <- sub("\\.csv$", "", twin_id)
  
  if(verbose) cat("\nProcessing twin", twin_id, "\n")
  
  # Read the CSV file 
  data <- tryCatch({
    read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
  }, error = function(e) {
    warning(paste("Error reading file", file_path, ":", e$message))
    return(NULL)
  })
  
  # Check if data was loaded successfully
  if(is.null(data) || nrow(data) < 3) {
    warning(paste("Insufficient data in file:", file_path))
    return(NULL)
  }
  
  if(verbose) cat("Raw data rows:", nrow(data), "\n")
  
  # Print column names to help with debugging
  if(verbose) cat("Original column names:", paste(colnames(data), collapse=", "), "\n")
  
  # Check if we have the expected columns
  required_cols <- c("ID", "time", "HR")
  # Make column names case-insensitive by converting to uppercase for checking
  data_cols_upper <- toupper(colnames(data))
  required_cols_upper <- toupper(required_cols)
  
  missing_cols <- required_cols[!required_cols_upper %in% data_cols_upper]
  if(length(missing_cols) > 0) {
    warning(paste("Missing required columns in file", file_path, ":", paste(missing_cols, collapse=", ")))
    if(verbose) cat("Available columns:", paste(colnames(data), collapse=", "), "\n")
    return(NULL)
  }
  
  # Get the actual column names that match our requirements (preserving original case)
  id_col <- colnames(data)[toupper(colnames(data)) == "ID"]
  time_col <- colnames(data)[toupper(colnames(data)) == "TIME"]
  hr_col <- colnames(data)[toupper(colnames(data)) == "HR"]
  
  # Create a simplified dataset with just the columns we need
  simplified_data <- data.frame(
    Twin_ID = data[[id_col]],
    time = as.numeric(data[[time_col]]),
    HR = as.numeric(data[[hr_col]])
  )
  
  if(verbose) cat("Created simplified dataset with columns: Twin_ID, time, HR\n")
  
  # Check for NAs in critical columns
  na_counts <- c(
    Twin_ID = sum(is.na(simplified_data$Twin_ID)),
    time = sum(is.na(simplified_data$time)),
    HR = sum(is.na(simplified_data$HR))
  )
  
  if(verbose) cat("NA counts in critical columns:", paste(names(na_counts), na_counts, sep="=", collapse=", "), "\n")
  
  # Remove rows with NAs in critical columns
  orig_rows <- nrow(simplified_data)
  simplified_data <- simplified_data[!is.na(simplified_data$time) & !is.na(simplified_data$HR), ]
  if(verbose) cat("Removed", orig_rows - nrow(simplified_data), "rows with NAs. Remaining:", nrow(simplified_data), "\n")
  
  # Fill in missing Twin_ID values with the filename
  if(sum(is.na(simplified_data$Twin_ID)) > 0 || all(simplified_data$Twin_ID == "")) {
    if(verbose) cat("Filling in missing Twin_ID values with", twin_id, "\n")
    simplified_data$Twin_ID <- twin_id
  }
  
  # Check if we still have enough data
  if(nrow(simplified_data) < 3) {
    warning(paste("Insufficient valid data points in file:", file_path))
    return(NULL)
  }
  
  # Fit cosinor model
  if(verbose) cat("Running cosinor analysis with", nrow(simplified_data), "data points\n")
  cosinor_results <- compute_cosinor(simplified_data)
  
  # If results are NULL, return NULL
  if(is.null(cosinor_results)) {
    if(verbose) cat("Cosinor analysis failed for twin", twin_id, "\n")
    return(NULL)
  } else {
    if(verbose) cat("Cosinor analysis successful for twin", twin_id, "\n")
  }
  
  return(cosinor_results)
}

# Main script execution

# Set the directory path where this script is running
dir_to_RStudio_script <- dirname(this.path::this.path())
# The RStudio folder is at same level as "01-input".  Hence the following ".."
input_path <- file.path(dir_to_RStudio_script, "..", "01-input", "AllTwinData")
# This is where output files will go
output_path <- file.path(dir_to_RStudio_script, "..", "02-output")

cat("Input Path:", input_path, "\n")
cat("Output Path:", output_path, "\n")

# First, convert all Excel files to CSV
cat("\nStep 1: Converting Excel files to CSV\n")
convert_excel_to_csv(input_path)

# Now list all CSV files in the directory (including newly converted ones)
cat("\nStep 2: Processing CSV files\n")
csv_files <- list.files(input_path, pattern = "\\.csv$", full.names = TRUE)
cat("Found", length(csv_files), "CSV files to process\n")

# Check if any CSV files were found
if(length(csv_files) == 0) {
  stop("No CSV files found in directory. Please check that Excel conversion worked correctly.")
}

# Process all files and combine results
all_results <- lapply(csv_files, function(file) {
  cat("\nProcessing file:", basename(file), "\n")
  result <- process_twin_data(file)
  if(!is.null(result)) {
    as.data.frame(result)
  } else {
    cat("Failed to process file:", basename(file), "\n")
    NULL
  }
})

# Print results summary before removing NULL values
cat("\nResults summary:\n")
for(i in 1:length(csv_files)) {
  file_name <- basename(csv_files[i])
  result_status <- if(is.null(all_results[[i]])) "Failed" else "Success"
  cat(file_name, ": ", result_status, "\n", sep="")
}

# Remove NULL results and combine
all_results <- all_results[!sapply(all_results, is.null)]

# Print the count of successful analyses
cat("\nNumber of successful analyses:", length(all_results), "\n")

# Check if any results exist
if(length(all_results) > 0) {
  # Convert to a dataframe
  results_df <- do.call(rbind, all_results)
  
  # Reorder columns as requested: Twin_ID, PercentRhythm, MESOR, Amplitude, Acrophase_degrees, Pvalue
  results_df <- results_df[, c("Twin_ID", "PercentRhythm", "MESOR", "Amplitude", 
                               "Acrophase_degrees", "Pvalue", "F_statistic")]
  
  # Write results to CSV
  output_file <- file.path(output_path, "twin_hr_cosinor_results.csv")
  write.csv(results_df, output_file, row.names = FALSE)
  
  # Print success message
  cat("\nSuccessfully created", output_file, "with cosinor analysis results\n")
  cat("Results contain data for", nrow(results_df), "twins\n")
} else {
  cat("\nERROR: No valid results found. Please check your input files.\n")
  cat("Possible issues:\n")
  cat("1. Required columns (ID, time, HR) not found\n")
  cat("2. Not enough valid data points after removing NA values\n")
  cat("3. Cosinor model failed to fit\n")
}

# Print final diagnostic information
cat("\nSummary:\n")
cat("Total files processed:", length(csv_files), "\n")
cat("Successful analyses:", length(all_results), "\n")
cat("Failed analyses:", length(csv_files) - length(all_results), "\n")
#Confirming I made a change to this line to test Git
