# library imports
library(data.table)
library(dplyr)
library(readxl)
library(tibble)

rm(list=ls())
code_directory <- "/Users/morgenbauer/Desktop/PHSL4094"

file_path <- paste(code_directory, "TwinsGroupAllData.xlsx", sep="/")

twinsdata <- read_excel(file_path)
groupedtwinsdata <- twinsdata %>% group_by(ID)

# Function to run ANOVA analysis on a specified subpopulation using specified model
runanova <- function(subpopulation, rhythmcharacteristics, data, formula_str) {
  # Convert formula string to actual formula
  model_formula <- as.formula(formula_str)
  
  # Run ANOVA
  a <- aov(model_formula, data=data)
  a_summary <- summary(a)
  
  # Print the summary
  print(paste("ANOVA Summary for", rhythmcharacteristics, "in", subpopulation))
  print(a_summary)
  
  # Calculate ICC
  MSBetweenGroups <- a_summary[[1]][1,3]
  MSWithinGroups <- a_summary[[1]][2,3]
  rl <- (MSBetweenGroups - MSWithinGroups) / (MSBetweenGroups + MSWithinGroups)
  
  # Calculate F-statistic and p-value
  F_stat <- MSBetweenGroups/MSWithinGroups
  p_value <- a_summary[[1]][1,5]
  
  # Get number of twin pairs
  k <- length(unique(data$ID))
  
  # Create results dataframe
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
  
  # Print the key results
  print(paste(rhythmcharacteristics, "in", subpopulation, "- ICC (rl):", rl))
  print(paste(rhythmcharacteristics, "in", subpopulation, "- p-value:", p_value))
  
  return(results)
}

# Function to run both amplitude and mesor analyses for a given subpopulation
runAmplitudeAndMesorAnovas <- function(subpopulation, data) {
  print(paste("Running analysis for", subpopulation))
  
  # Run amplitude analysis
  amplitude_results <- runanova(subpopulation, "Amplitude", data, "A24h ~ ID")
  
  # Run mesor analysis
  mesor_results <- runanova(subpopulation, "MESOR", data, "M24h ~ ID")
  
  # Return both results
  return(list(amplitude = amplitude_results, mesor = mesor_results))
}

# Execute the analyses for all data
all_data_results <- runAmplitudeAndMesorAnovas("AllData", groupedtwinsdata)

# Low body weight group analysis - using the exact subsetting syntax 
low_bw_data <- groupedtwinsdata[groupedtwinsdata$Gp=="1L",]
low_bw_results <- runAmplitudeAndMesorAnovas("LowBodyWeight", low_bw_data)

# Medium body weight group analysis - using the exact subsetting syntax 
medium_bw_data <- groupedtwinsdata[groupedtwinsdata$Gp=="2M",]
medium_bw_results <- runAmplitudeAndMesorAnovas("MediumBodyWeight", medium_bw_data)

# High body weight group analysis - using the exact subsetting syntax 
high_bw_data <- groupedtwinsdata[groupedtwinsdata$Gp=="3H",]
high_bw_results <- runAmplitudeAndMesorAnovas("HighBodyWeight", high_bw_data)

# Create a consolidated results dataframe
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

# Save the consolidated results to a CSV file
consolidated_results_file <- paste(code_directory, "consolidated_analysis_results.csv", sep="/")
write.csv(consolidated_results, consolidated_results_file, row.names = FALSE)

# Print summary of results
print("Summary of ICC Analysis Results:")
print(consolidated_results)

# Create a more readable group factor for plotting
twinsdata$GroupLabel <- factor(twinsdata$Gp, 
                               levels = c("1L", "2M", "3H"), 
                               labels = c("Low BW", "Medium BW", "High BW"))

# Create combined boxplot for A24h by weight group
png(paste(code_directory, "combined_a24h_boxplot.png", sep="/"), width=800, height=600)
boxplot(A24h ~ GroupLabel, data=twinsdata, 
        main="Amplitude (A24h) by Weight Group", 
        xlab="Weight Group", 
        ylab="Amplitude (A24h)",
        col=c("lightcoral", "lightgreen", "lightblue"))
dev.off()

# Create combined boxplot for M24h by weight group
png(paste(code_directory, "combined_m24h_boxplot.png", sep="/"), width=800, height=600)
boxplot(M24h ~ GroupLabel, data=twinsdata, 
        main="MESOR (M24h) by Weight Group", 
        xlab="Weight Group", 
        ylab="MESOR (M24h)",
        col=c("lightcoral", "lightgreen", "lightblue"))
dev.off()

# Create combined histograms of residuals for A24h
png(paste(code_directory, "combined_a24h_residuals.png", sep="/"), width=1000, height=800)
par(mfrow=c(2,2))

# All data residuals
hist(all_data_results$amplitude$residuals, 
     main="A24h Residuals - All Data", 
     xlab="Residuals", 
     col="lightgray", 
     border="black")

# Low BW residuals
hist(low_bw_results$amplitude$residuals, 
     main="A24h Residuals - Low BW", 
     xlab="Residuals", 
     col="lightcoral", 
     border="black")

# Medium BW residuals
hist(medium_bw_results$amplitude$residuals, 
     main="A24h Residuals - Medium BW", 
     xlab="Residuals", 
     col="lightgreen", 
     border="black")

# High BW residuals
hist(high_bw_results$amplitude$residuals, 
     main="A24h Residuals - High BW", 
     xlab="Residuals", 
     col="lightblue", 
     border="black")

par(mfrow=c(1,1))
dev.off()

# Create combined histograms of residuals for M24h
png(paste(code_directory, "combined_m24h_residuals.png", sep="/"), width=1000, height=800)
par(mfrow=c(2,2))

# All data residuals
hist(all_data_results$mesor$residuals, 
     main="M24h Residuals - All Data", 
     xlab="Residuals", 
     col="lightgray", 
     border="black")

# Low BW residuals
hist(low_bw_results$mesor$residuals, 
     main="M24h Residuals - Low BW", 
     xlab="Residuals", 
     col="lightcoral", 
     border="black")

# Medium BW residuals
hist(medium_bw_results$mesor$residuals, 
     main="M24h Residuals - Medium BW", 
     xlab="Residuals", 
     col="lightgreen", 
     border="black")

# High BW residuals
hist(high_bw_results$mesor$residuals, 
     main="M24h Residuals - High BW", 
     xlab="Residuals", 
     col="lightblue", 
     border="black")

par(mfrow=c(1,1))
dev.off()

print("Analysis completed. Files saved:")
print(paste("1. Consolidated Results CSV:", consolidated_results_file))
print(paste("2. Combined A24h Boxplot:", code_directory, "/combined_a24h_boxplot.png", sep=""))
print(paste("3. Combined M24h Boxplot:", code_directory, "/combined_m24h_boxplot.png", sep=""))
print(paste("4. Combined A24h Residuals Histograms:", code_directory, "/combined_a24h_residuals.png", sep=""))
print(paste("5. Combined M24h Residuals Histograms:", code_directory, "/combined_m24h_residuals.png", sep=""))