install.packages("randomForest")
library(randomForest)
 #Change the dataformat
colnames(df) <- gsub("-", "_", colnames(df))
colnames(df) <- gsub("\\(|\\)", "", colnames(df))
colnames(df) <- gsub("%", "per", colnames(df))
colnames(df) <- gsub("#", "n", colnames(df))

#start with regression
# Convert detect_hy to a binary format
df$detect_hy_binary <- ifelse(df$detect_hy == "YES", 1, 0)

# Fit a logistic regression model
logit_model <- glm(detect_hy_binary ~ n_input_reads_ilu + n_of_SARS_CoV_2_reads_ilu + per_SARS_CoV_2_reads_ilu + 
                     per_Trimmed_reads_ilu + General_Mean_Quality_FastQC_ilu + n_input_reads_ont + 
                     n_of_SARS_CoV_2_reads_ont + per_SARS_CoV_2_reads_ont + per_Trimmed_reads_ont + 
                     General_Mean_Quality_FastQC_ont, 
                   family = binomial, data = df)
summary(logit_model)

library(rpart)
tree_model <- rpart(detect_hy_binary ~ n_input_reads_ilu + n_of_SARS_CoV_2_reads_ilu + per_SARS_CoV_2_reads_ilu + 
                      per_Trimmed_reads_ilu + General_Mean_Quality_FastQC_ilu + n_input_reads_ont + 
                      n_of_SARS_CoV_2_reads_ont + per_SARS_CoV_2_reads_ont + per_Trimmed_reads_ont + 
                      General_Mean_Quality_FastQC_ont, 
                    data = df, method="class")
plot(tree_model, uniform=TRUE, margin=0.1)
text(tree_model, use.n=TRUE, all=TRUE, cex=.8)

# Determine the best cp (complexity parameter)
printcp(tree_model)

# Plot the cross-validated error
plotcp(tree_model)

# Choose a complexity parameter and prune the tree
pruned_tree <- prune(tree_model, cp = tree_model$cptable[which.min(tree_model$cptable[,"xerror"]),"CP"])

# Visualize the pruned tree
plot(pruned_tree, uniform=TRUE, margin=0.1)
text(pruned_tree, use.n=TRUE, all=TRUE, cex=.8)

#####################
#random_forest_part2
# Convert the variable to factor
df$detect_hy_binary <- as.factor(df$detect_hy_binary)

# Fit the Random Forest model again
rf_model <- randomForest(detect_hy_binary ~ n_input_reads_ilu + n_of_SARS_CoV_2_reads_ilu + per_SARS_CoV_2_reads_ilu + 
                           per_Trimmed_reads_ilu + General_Mean_Quality_FastQC_ilu + n_input_reads_ont + 
                           n_of_SARS_CoV_2_reads_ont + per_SARS_CoV_2_reads_ont + per_Trimmed_reads_ont + 
                           General_Mean_Quality_FastQC_ont, 
                         data = df, ntree=100)

# Print the model summary
print(rf_model)

####END OF IMPORTANT SECTION FROM NOW ON IT IS EXTRA####

#random_forest_start
df_clean_hy <- df[df$detect_hy == "YES",]
dnv_cols <- grep("_dnv", names(df_clean_hy), value = TRUE)
df_clean_no_dnv <- df_clean_hy[ , !(names(df_clean_hy) %in% dnv_cols)]

metric_columns <- grep("_hybrid$", names(df_clean_no_dnv), value = TRUE)

run_rf_model <- function(metric_col_name, data) {
  # Define the formula dynamically based on the metric column
  formula <- as.formula(paste(metric_col_name, "~", paste(colnames(data)[2:13], collapse=" + ")))
  
  # Determine if the metric should be treated as classification
  if (metric_col_name == "AA_del_unique_8_hybrid" && length(unique(data[[metric_col_name]])) <= 5) {
    model <- randomForest(formula, data=data, ntree=500, regression=FALSE)
  } else {
    model <- randomForest(formula, data=data, ntree=500)
  }
  
  return(model)
}

models_list <- list()

for(metric in metric_columns) {
  cat("Training model for", metric, "\n")
  models_list[[metric]] <- run_rf_model(metric, df_clean_no_dnv)
}

importance(rf_model)


##start with other methods

#data prep

method_properties <- list(
  "ilu_ori_sph" = list(suffix = "ilu_ori_sph", properties = 2:6),
  "ilu_alt_RKI" = list(suffix = "ilu_alt_RKI", properties = 2:6),
  "ilu_alt2_man" = list(suffix = "ilu_alt2_man", properties = 2:6),
  "ont_ori_mbi" = list(suffix = "ont_ori_mbi", properties = 7:11),
  "ont_alt_ART" = list(suffix = "ont_alt_ART", properties = 7:11),
  "ont_alt2_man" = list(suffix = "ont_alt2_man", properties = 7:11),
  "hybrid" = list(suffix = "hybrid", properties = 2:11)
)

method_properties <- list(
  ilu_ori_sph = 2:6,
  ilu_alt_RKI = 2:6,
  ilu_alt2_man = 2:6,
  ont_ori_mbi = 7:11,
  ont_alt_ART = 7:11,
  ont_alt2_man = 7:11,
  hybrid = 2:11
)

# Function to determine which metric columns to use based on the method's suffix
get_metric_columns <- function(suffix) {
  grep(paste0("_", suffix, "$"), names(df_clean), value = TRUE)
}

# Function to filter samples based on those detected by the hybrid method
filter_samples <- function(df_clean, method) {
  hybrid_samples <- unique(df_clean_hy$Sample)
  df_filtered <- df_clean[df_clean$Sample %in% hybrid_samples, ]
  return(df_filtered)
}

# Updated run_rf_model to accommodate the method's properties and handle NA values
run_rf_model <- function(metric_col_name, data, properties) {
  # Define the formula dynamically based on the metric column and method's properties
  formula <- as.formula(paste(metric_col_name, "~", paste(colnames(data)[properties], collapse=" + ")))
  
  # Check if the metric should be treated as classification
  if (metric_col_name == "AA_del_unique_8_hybrid" && length(unique(data[[metric_col_name]])) <= 5) {
    model <- randomForest(formula, data=data, ntree=500, regression=FALSE, na.action=na.omit)
  } else {
    model <- randomForest(formula, data=data, ntree=500, na.action=na.omit)
  }
  
  return(model)
}

# Updated run_models_for_method to use properties and filter samples
run_models_for_method <- function(df_clean, method, method_suffix, properties) {
  data <- filter_samples(df_clean, method)
  metric_columns <- get_metric_columns(method_suffix)
  
  models_list <- list()
  
  for(metric in metric_columns) {
    cat("Training model for", metric, "using", method_suffix, "data\n")
    models_list[[metric]] <- run_rf_model(metric, data, properties)
  }
  
  return(models_list)
}

# Master function to manage the entire process
run_all_models <- function(df_clean) {
  all_models <- list()
  for(method in names(method_properties)) {
    all_models[[method]] <- run_models_for_method(df_clean, method, method_properties[[method]]$suffix, method_properties[[method]]$properties)
  }
  return(all_models)
}

# Execute
all_models <- run_all_models(df_clean)

#Model Evaluation:

# Loop through each method and metric to print the metrics
for(method in names(all_models)) {
  for(metric in names(all_models[[method]])) {
    model <- all_models[[method]][[metric]]
    
    if(model$type == "classification") {
      oob_error <- round(model$err.rate[nrow(model$err.rate), 1], 4)
      cat(paste("Method:", method, "| Metric:", metric, "| OOB Error Rate:", oob_error, "\n"))
    } else if(model$type == "regression") {
      mse_final <- round(model$mse[length(model$mse)], 4)
      rsq_final <- round(model$rsq[length(model$rsq)], 4)
      cat(paste("Method:", method, "| Metric:", metric, "| MSE:", mse_final, "| R-squared:", rsq_final, "\n"))
    }
  }
}

plot_importance <- function(model) {
  importance_data <- importance(model)
  importance_df <- data.frame(Variables = rownames(importance_data), Importance = importance_data[, "MeanDecreaseGini"])
  ggplot(importance_df, aes(x = reorder(Variables, Importance), y = Importance)) + 
    geom_bar(stat = "identity") +
    coord_flip() + 
    labs(title = "Feature Importance")
}


library(ggplot2)

# Initialize empty vectors to store results for visualization
methods_vec <- c()
metrics_vec <- c()
rsq_values <- c()
mse_values <- c()

results_list <- list()

# Loop through each method and metric to extract metrics
for(method in names(all_models)) {
  for(metric in names(all_models[[method]])) {
    model <- all_models[[method]][[metric]]
    
    temp_data <- list()
    
    temp_data$Method <- method
    temp_data$Metric <- metric
    
    if(model$type == "classification") {
      temp_data$R2 <- NA
      temp_data$MSE <- NA
    } else if(model$type == "regression") {
      temp_data$MSE <- model$mse[length(model$mse)]
      temp_data$R2 <- model$rsq[length(model$rsq)]
    }
    
    results_list[[length(results_list) + 1]] <- temp_data
  }
}
#
# Convert the results list to a dataframe
results_df <- do.call(rbind.data.frame, results_list)
# Extract the metric name, assuming the metric names are before the first underscore.
results_df$MetricName <- gsub("_.+$", "", results_df$Metric)
results_df <- na.omit(results_df)

plot_R2_by_method <- function(results_df) {
  # Correcting the MetricName
  results_df$MetricName <- gsub("_.+$", "", results_df$Metric)
  
  unique_metrics <- unique(results_df$MetricName)
  
  pdf("R2_by_Metric.pdf")
  for(metric in unique_metrics) {
    metric_data <- subset(results_df, MetricName == metric)
    
    p <- ggplot(metric_data, aes(x = Method, y = R2, fill = Method)) + 
      geom_bar(stat = "identity") +
      labs(title = paste("R-squared Values by Method for", metric), y = "R-squared Value") +
      theme_minimal() +
      theme(legend.position = "none")
    
    print(p)
  }
  dev.off()
}

plot_R2_by_method(results_df)

########
plot_R2_by_method <- function(results_df) {
  # Correcting the MetricName
  results_df$MetricName <- gsub("_.+$", "", results_df$Metric)
  
  unique_metrics <- unique(results_df$MetricName)
  
  # Define fixed y-axis limits
  y_min <- min(results_df$R2) - 0.1 * abs(min(results_df$R2))
  y_max <- max(results_df$R2) + 0.1 * abs(max(results_df$R2))
  
  pdf("R2_by_Metric.pdf")
  for(metric in unique_metrics) {
    metric_data <- subset(results_df, MetricName == metric)
    
    p <- ggplot(metric_data, aes(x = Method, y = R2, fill = Method)) + 
      geom_bar(stat = "identity") +
      labs(title = paste("R-squared Values by Method for", metric), y = "R-squared Value") +
      theme_minimal() +
      theme(legend.position = "none") +
      coord_cartesian(ylim = c(y_min, y_max))  # set fixed y-axis range here
    
    print(p)
  }
  dev.off()
}

plot_R2_by_method(results_df)



****************

library(ggplot2)

unique_metrics <- unique(results_df$MetricName)

pdf("R2_by_Metric.pdf")
for(metric in unique_metrics) {
  metric_data <- subset(results_df, MetricName == metric)
  
  p <- ggplot(metric_data, aes(x = Method, y = R2, fill = Method)) + 
    geom_bar(stat = "identity") +
    labs(title = paste("R-squared Values by Method for", metric), y = "R-squared Value") +
    theme_minimal() +
    theme(legend.position = "none")
  
  print(p)
}
dev.off()


# Assuming "Genomic_Fraction" is one of the metrics
metric_data <- subset(results_df, grepl("genomic_fraction", Metric))

# Plot R^2 for Genomic_Fraction
ggplot(metric_data, aes(x = Method, y = R2, fill = Method)) + 
  geom_bar(stat = "identity") +
  labs(title = "R-squared Values by Method for Genomic_Fraction", y = "R-squared Value") +
  theme_minimal() +
  theme(legend.position = "none")


# List of unique metric names
metric_names = unique(results_df$MetricName)

# For each metric name, generate a bar plot
for(metric_name in metric_names){
  metric_data <- subset(results_df, MetricName == metric_name)
  
  ggplot(metric_data, aes(x = Method, y = R2, fill = Method)) + 
    geom_bar(stat = "identity") +
    labs(title = paste("R-squared Values by Method for", metric_name), y = "R-squared Value") +
    theme_minimal() +
    theme(legend.position = "none") +
}

