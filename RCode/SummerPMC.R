
library(readxl)
# Directory Setup
dir_to_script <- getwd()
# dir_to_script <- file.path("C:", "Users", "chase", "Documents", "HCC", "_git", "MBauer-HCC25", "RCode" )
if (!dir.exists(dir_to_script)) stop("Directory not found: ", dir_to_script)

input_path <- file.path(dir_to_script, "..", "GC")
if (!dir.exists(input_path)) stop("directory not found: ", input_path)

output_path <- file.path(dir_to_script, "..", "02-output")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

# Obtain a list of files in the input_path
xlsxFileToProcess <- list.files(
  path       = input_path, 
  pattern    = "cosall.xlsx", 
  full.names = TRUE
)

sheetNameFromPath <- function(aPath){
  fs::path_ext_remove(basename(aPath))
}

all_results <- list()

# Process a sheet
process_sheet_data <- function(file_path, sheet_name) {
  cat("\nProcessing sheet:", sheet_name, "\n")
  data <- tryCatch({
    read_excel(file_path, sheet = sheet_name)
  }, error = function(e) {
    warning(paste("Error reading sheet:", e$message))
    return(NULL)
  })
}

input_data <- process_sheet_data(xlsxFileToProcess, "cosall")



