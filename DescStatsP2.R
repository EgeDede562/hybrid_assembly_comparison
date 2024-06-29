###DESCRIPTIVE STATISTICS - PART II###
##This code snippet is to calculate and graph descriptive statistics for each metric-method combination##
. 
# Load the necessary library. 
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

# Get the column names for the metric_method combinations from the original dataframe
metric_method_names <- colnames(df)[13:92]
metric_method_names_adjusted <- str_replace_all(metric_method_names, "_\\d+_", "_")

# Generate descriptive statistics for each metric-method combination
desc_stats_df <- df %>%
  pivot_longer(
    cols = matches("^((genomic_fraction_)|(AA_subs_)|(N_subs_)|(nextclade_qval_)|(AA_dels_)|(N_dels_)|(AA_unique_)|(N_unique_subs_)|(N_del_unique)|(AA_del_unique_)).*"),
    names_to = "metric_method",
    values_to = "value"
  ) %>%
  # Adjusted separate call
  separate(metric_method, into = c("metric", "method"), sep = "_\\d+_", extra = "merge") %>%
  mutate(
    metric = str_replace(metric, "_\\d+$", ""), # Remove trailing number if any from metric
    method = str_replace(method, "^\\d+_", "") # Remove leading number if present in method
  ) %>%
  group_by(metric, method) %>%
  summarize(
    Mean = mean(value, na.rm = TRUE),
    Median = median(value, na.rm = TRUE),
    SD = sd(value, na.rm = TRUE),
    Min = min(value, na.rm = TRUE),
    Max = max(value, na.rm = TRUE),
    .groups = "drop" # Ensures the result is ungrouped
  )

# Back-up Descriptive Statistics Dataframe
desc_stats_dfbu <- desc_stats_df

# EXAMPLE FOR A SINGLE METRIC VISUALIZATION
# This plot is not ordered. 
single_metric_plot <- ggplot(filter(desc_stats_df, metric == "AA_subs"), aes(x = method, y = Mean, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Mean AA Substitutions by Method", x = "Method", y = "Mean Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Explicitly print the plot
print(single_metric_plot)
