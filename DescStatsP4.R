#### DESCRIPTIVE STATISTICS ###

# Load necessary libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)

# Ensure gridExtra is installed
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")

# Define the methods based on the column names
methods <- c("ilu_ori_sph", "ilu_alt_RKI", "ilu_alt2_man", 
             "ont_ori_mbi", "ont_alt_ART", "ont_alt2_man", 
             "ilu_dnv_spa", "hybrid")

# List of the performance metrics
metrics <- c("genomic_fraction_", "AA_subs_", "N_subs_", "nextclade_qval_", 
             "AA_dels_", "N_dels_", "AA_unique_", "N_unique_subs_", 
             "N_del_unique_", "AA_del_unique_")

desc_stats_list <- list()

# Compute the descriptive statistics
for (metric in metrics) {
  metric_cols <- grep(metric, names(df), value = TRUE)
  metric_stats <- data.frame()
  for (col in metric_cols) {
    stats <- c(mean(df[[col]], na.rm = TRUE), 
               median(df[[col]], na.rm = TRUE), 
               sd(df[[col]], na.rm = TRUE), 
               min(df[[col]], na.rm = TRUE), 
               max(df[[col]], na.rm = TRUE))
    metric_stats <- rbind(metric_stats, stats)
  }
  colnames(metric_stats) <- c("Mean", "Median", "Standard Deviation", "Min", "Max")
  rownames(metric_stats) <- metric_cols
  desc_stats_list[[metric]] <- metric_stats
}

# Combine the list into a dataframe
desc_stats <- do.call(rbind, desc_stats_list)

# Display the descriptive statistics
print(desc_stats)

## VISUALIZATION OF DESC STATISTICS

# Let's plot for the AA_subs_ metric as an example
# Extract the rows related to the AA_subs_ metric
subset_metric <- desc_stats[grepl("AA_subs_", rownames(desc_stats)), ]

# Reset row names for visualization
subset_metric <- subset_metric %>%
  mutate(Method = rownames(subset_metric))

# Convert to long format
long_data <- subset_metric %>%
  pivot_longer(cols = c("Mean", "Median", "Standard Deviation", "Min", "Max"), 
               names_to = "Statistic", 
               values_to = "Value")

# Create the bar plot
plot <- ggplot(long_data, aes(x = Method, y = Value, fill = Statistic)) +
  geom_bar(stat = "identity", position = "dodge") + 
  theme_minimal() +
  labs(title = "Descriptive Statistics for AA_subs Metric by Method",
       x = "Method",
       y = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot)

# Generate boxplots for each metric
generate_plot <- function(metric) {
  subset_data <- df[, grep(metric, colnames(df))]
  
  # Convert to long format
  long_subset_data <- pivot_longer(subset_data, cols = everything(), names_to = "Method", values_to = "Value")
  
  # Extract only the method part of the column names
  methods <- gsub(metric, "", colnames(subset_data))
  long_subset_data$Method <- factor(long_subset_data$Method, levels = colnames(subset_data), labels = methods)
  
  plot <- ggplot(long_subset_data, aes(x = Method, y = Value, fill = Method)) +
    geom_boxplot() + 
    theme_minimal() +
    labs(title = paste("Distribution of", metric, "by Method"),
         x = "Method",
         y = "Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  return(plot)
}

# Generate list of plots
plots_list <- lapply(metrics, generate_plot)

# Arrange plots into a grid
grid.arrange(grobs = plots_list, ncol = 5)
