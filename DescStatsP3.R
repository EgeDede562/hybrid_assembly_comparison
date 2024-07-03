###DESCRIPTIVE STATISTICS III - STATISTICAL ANALYSIS AND VISUALIZATION###

#PART 3.1 - Load necessary libraries
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)
library(FSA)
## END OF PART 3.1 ##


## PART 3.2 - Normality testing ##
# Subset the data for the metric
genomic_fraction_data <- df[, grepl("genomic_fraction_", colnames(df))]

# Shapiro-Wilk test for normality
shapiro_results <- apply(genomic_fraction_data, 2, shapiro.test)

# Extract p-values from the test results
shapiro_pvalues <- sapply(shapiro_results, function(x) x$p.value)

print(shapiro_pvalues)
## END OF PART 3.2 ##


## PART 3.3 - Kruskal-Wallis Test for Multiple Metrics ##
metrics_names <- c("genomic_fraction_", "AA_subs_", "N_subs_", "nextclade_qval_", 
                   "AA_dels_", "N_dels_", "AA_unique_", "N_unique_subs_", "N_del_unique", "AA_del_unique")

kruskal_results <- lapply(metrics_names, function(metric) {
  # Subset data
  metric_data <- df[, grepl(pattern = metric, colnames(df))]
  
  # Convert to long format
  long_data <- pivot_longer(metric_data, cols = everything(), names_to = "Method", values_to = "Value")
  
  # Apply Kruskal-Wallis test
  kruskal.test(Value ~ Method, data = long_data)
})

# Store the p-values for each metric
p_values <- sapply(kruskal_results, function(result) result$p.value)

# Display p-values for each metric
data.frame(Metric = metrics_names, P_Value = p_values)
##END OF PART 3.3##


##PART 3.4 - Dunn's Test##
# Function to apply Dunn's test and handle potential issues
run_dunns_test <- function(metric) {
  metric_data <- df[, grepl(pattern = metric, colnames(df))]
  long_data <- pivot_longer(metric_data, cols = everything(), names_to = "Method", values_to = "Value")
  
  # Remove rows with NA values in 'Value'
  long_data <- na.omit(long_data)
  
  # Apply Dunn's test if there are enough non-NA values
  if (nrow(long_data) > 0) {
    dunnTest(Value ~ Method, data = long_data, method = "bonferroni")
  } else {
    NULL
  }
}

# Apply Dunn's test to each metric
dunns_results <- lapply(metrics_names, run_dunns_test)
 
dunns_results <- lapply(metrics_names, function(metric) {
  metric_data <- df[, grepl(pattern = metric, colnames(df))]
  long_data <- pivot_longer(metric_data, cols = everything(), names_to = "Method", values_to = "Value")
  dunnTest(Value ~ Method, data = long_data, method = "bonferroni")
})

# Process Dunn's test results
dunns_comparison <- lapply(dunns_results, function(result) {
  result$res
})

# Combine results for all metrics into a single dataframe
all_dunns_results <- do.call(rbind, lapply(1:length(metrics_names), function(i) {
  cbind(Metric = metrics_names[i], dunns_comparison[[i]])
}))

# Filter significant results
significant_results <- subset(all_dunns_results, P.adj < 0.05)
print(head(significant_results))
## END OF PART 3.4 ##


## PART 3.5 - Visualizing Dunn's Test Results ##
# Separate the comparison into two columns: Method1 and Method2
heatmap_data <- significant_results %>% 
  separate(Comparison, into = c("Method1", "Method2"), sep = "-") %>%
  filter(P.adj < 0.05)

#map the methods
method_map <- c(
  "1" = "ilu_ori_sph",
  "2" = "ilu_alt_RKI",
  "3" = "ilu_alt2_man",
  "4" = "ont_ori_mbi",
  "5" = "ont_alt_ART",
  "6" = "ont_alt2_man",
  "7" = "ilu_dnv_spa",
  "8" = "hybrid"
)

# Simplify the method names
heatmap_data$Method1 <- method_map[gsub(".*_([0-9]+).*", "\\1", heatmap_data$Method1)]
heatmap_data$Method2 <- method_map[gsub(".*_([0-9]+).*", "\\1", heatmap_data$Method2)]

# Define significance markers based on adjusted p-values
heatmap_data$Significance <- cut(heatmap_data$P.adj,
                                 breaks = c(0, 0.001, 0.01, 0.05, 1),
                                 labels = c("***", "**", "*", "ns"),
                                 include.lowest = TRUE, right = FALSE)

# Generate heatmaps
heatmap_plot <- ggplot(heatmap_data, aes(x = Method1, y = Metric, fill = Z, label = Significance)) + 
  geom_tile() + 
  geom_text(size = 3) + 
  facet_wrap(~ Method2, ncol = 3, scales = "free_y") + 
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "white", colour = "grey", linewidth = 1), 
    strip.text = element_text(size = 12)
  ) + 
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red", name = "Z Value") + 
  labs(title = "Dunn's Test Results")

print(heatmap_plot)
##END OF PART 3.5##