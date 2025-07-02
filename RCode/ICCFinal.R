# install.packages("readxl")
library(data.table)
library(dplyr)
library(readxl)
library(tibble)

rm(list=ls())
code_directory <- "/Users/morgenbauer/Desktop/PHSL4094"

file_path <- paste(code_directory, "TwinsGroupAllData.xlsx", sep="/")

twinsdata <- read_excel(file_path)
groupedtwinsdata <- twinsdata %>% group_by(ID)

# Run ANOVA for A24h
a <- aov(A24h ~ ID, data=groupedtwinsdata)
a_summary <- summary(a)

# Run ANOVA for M(24h)
b <- aov(`M(24h)` ~ ID, data=groupedtwinsdata)
b_summary <- summary(b)

# Removed the ANOVA for Phi24h (acrophases)

# Histograms and boxplots for A24h
hist(a$residuals)
boxplot(A24h ~ ID, data=groupedtwinsdata)

# Calculate ICC for A24h
MSBetweenGroups_A <- a_summary[[1]][1,3]
MSWithinGroups_A <- a_summary[[1]][2,3]
rl_A <- (MSBetweenGroups_A - MSWithinGroups_A) / (MSBetweenGroups_A + MSWithinGroups_A)

# Calculate ICC for M(24h)
MSBetweenGroups_M <- b_summary[[1]][1,3]
MSWithinGroups_M <- b_summary[[1]][2,3]
rl_M <- (MSBetweenGroups_M - MSWithinGroups_M) / (MSBetweenGroups_M + MSWithinGroups_M)

# Print full ANOVA summaries
print("ANOVA Summary for A24h:")
print(a_summary)

print("ANOVA Summary for M(24h):")
print(b_summary)

print(sprintf("A24h - MS Between Groups: %f", MSBetweenGroups_A))
print(sprintf("A24h - MS Within Groups: %f", MSWithinGroups_A))
print(sprintf("A24h - ICC (rl): %f", rl_A))

print(sprintf("M(24h) - MS Between Groups: %f", MSBetweenGroups_M))
print(sprintf("M(24h) - MS Within Groups: %f", MSWithinGroups_M))
print(sprintf("M(24h) - ICC (rl): %f", rl_M))

# Group the twins by body weight groups
weight_groups <- twinsdata %>%
  group_by(Gp) %>%
  summarise(
    count = n(),
    mean_A24h = mean(A24h, na.rm = TRUE),
    sd_A24h = sd(A24h, na.rm = TRUE),
    mean_M24h = mean(`M(24h)`, na.rm = TRUE),
    sd_M24h = sd(`M(24h)`, na.rm = TRUE)
  )

print("Summary by Body Weight Groups:")
print(weight_groups)

# Create a data frame with your ICC results for both metrics
icc_results <- data.frame(
  Metric = c("A24h - MS Between Groups", "A24h - MS Within Groups", "A24h - ICC (rl)",
             "M(24h) - MS Between Groups", "M(24h) - MS Within Groups", "M(24h) - ICC (rl)"),
  Value = c(MSBetweenGroups_A, MSWithinGroups_A, rl_A,
            MSBetweenGroups_M, MSWithinGroups_M, rl_M)
)

# Save the ANOVA summaries to CSV files
write.csv(as.data.frame(a_summary[[1]]), paste(code_directory, "a24h_anova_summary.csv", sep="/"), row.names = TRUE)
write.csv(as.data.frame(b_summary[[1]]), paste(code_directory, "m24h_anova_summary.csv", sep="/"), row.names = TRUE)

# Save the ICC results to a CSV file
write.csv(icc_results, paste(code_directory, "icc_results.csv", sep="/"), row.names = FALSE)

# Save the weight group summary
write.csv(weight_groups, paste(code_directory, "weight_group_summary.csv", sep="/"), row.names = FALSE)

# Save the histogram as a PNG file for A24h
png(paste(code_directory, "a24h_residuals_histogram.png", sep="/"), width=800, height=600)
hist(a$residuals, main="Histogram of A24h Residuals", xlab="Residuals", col="lightblue", border="black")
dev.off()

# Save the histogram as a PNG file for M(24h)
png(paste(code_directory, "m24h_residuals_histogram.png", sep="/"), width=800, height=600)
hist(b$residuals, main="Histogram of M(24h) Residuals", xlab="Residuals", col="lightblue", border="black")
dev.off()

# Save the boxplot as a PNG file for A24h
png(paste(code_directory, "a24h_boxplot.png", sep="/"), width=800, height=600)
boxplot(A24h ~ ID, data=groupedtwinsdata, main="Boxplot of A24h by ID", xlab="ID", ylab="A24h", col="lightgreen")
dev.off()

# Save the boxplot as a PNG file for M(24h)
png(paste(code_directory, "m24h_boxplot.png", sep="/"), width=800, height=600)
boxplot(`M(24h)` ~ ID, data=groupedtwinsdata, main="Boxplot of M(24h) by ID", xlab="ID", ylab="M(24h)", col="lightgreen")
dev.off()

# Create boxplots by weight group for both metrics
png(paste(code_directory, "a24h_by_weight_group.png", sep="/"), width=800, height=600)
boxplot(A24h ~ Gp, data=twinsdata, main="A24h by Weight Group", xlab="Weight Group", ylab="A24h", col=c("red", "green", "blue"))
dev.off()

png(paste(code_directory, "m24h_by_weight_group.png", sep="/"), width=800, height=600)
boxplot(`M(24h)` ~ Gp, data=twinsdata, main="M(24h) by Weight Group", xlab="Weight Group", ylab="M(24h)", col=c("red", "green", "blue"))
dev.off()

print("Files saved to your directory:")
print(paste("1. A24h ANOVA Summary: ", code_directory, "/a24h_anova_summary.csv", sep=""))
print(paste("2. M(24h) ANOVA Summary: ", code_directory, "/m24h_anova_summary.csv", sep=""))
print(paste("3. ICC Results: ", code_directory, "/icc_results.csv", sep=""))
print(paste("4. Weight Group Summary: ", code_directory, "/weight_group_summary.csv", sep=""))
print(paste("5. A24h Residuals Histogram: ", code_directory, "/a24h_residuals_histogram.png", sep=""))
print(paste("6. M(24h) Residuals Histogram: ", code_directory, "/m24h_residuals_histogram.png", sep=""))
print(paste("7. A24h Boxplot: ", code_directory, "/a24h_boxplot.png", sep=""))
print(paste("8. M(24h) Boxplot: ", code_directory, "/m24h_boxplot.png", sep=""))
print(paste("9. A24h by Weight Group: ", code_directory, "/a24h_by_weight_group.png", sep=""))
print(paste("10. M(24h) by Weight Group: ", code_directory, "/m24h_by_weight_group.png", sep=""))

