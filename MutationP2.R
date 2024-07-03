library(dplyr)
library(tidyr)
library(stringr)
library(readxl)

# Function to prepare the final dataframe
prepare_final_df <- function(data) {
  # Convert the tool columns to character
  data <- data %>%
    mutate(across(c("ilu_ori", "ilu_alt", "ont_ori", "ont_alt", "ilu_dnv", "hybrid", "ilu_alt_2", "ont_alt_2"), 
                  as.character))
  
  # Splitting the mutations
  data <- data %>%
    mutate(across(c("ilu_ori", "ilu_alt", "ont_ori", "ont_alt", "ilu_dnv", "hybrid", "ilu_alt_2", "ont_alt_2"), 
                  ~strsplit(., ","), 
                  .names = "{.col}_split"))
  
  data <- data %>%
    pivot_longer(cols = ends_with("_split"), 
                 names_to = "tool",
                 values_to = "mutation") %>%
    unnest(cols = mutation)
  
  # Clean up tool names
  data$tool <- gsub("_split", "", data$tool)
  
  mutation_count <- data %>%
    group_by(tool, mutation) %>%
    summarise(frequency = n(), .groups = "drop")
  
  # Pivot the data to wide format
  final_df <- mutation_count %>%
    pivot_wider(names_from = tool, values_from = frequency, values_fill = list(frequency = 0))
  
  # Ensure "mutation" is the first column and order the remaining columns
  tool_cols <- c("ilu_ori", "ilu_alt", "ont_ori", "ont_alt", "ilu_dnv", "hybrid", "ilu_alt_2", "ont_alt_2")
  final_df <- final_df %>% select(mutation, all_of(tool_cols))
  
  # Standardize mutation format to "Gene:Mutation"
  final_df$mutation <- str_replace_all(final_df$mutation, "_", ":")
  
  # Now extract the gene names from the standardized mutation strings
  final_df <- final_df %>%
    mutate(Gene = ifelse(str_detect(mutation, ":|_"), str_extract(mutation, "^[^:]+"), Gene),
           Mutation = ifelse(str_detect(mutation, ":|_"), str_extract(mutation, "(?<=:|_)(.*)$"), mutation))
  
  # Removing the trailing ":" or "_"
  final_df$Gene <- gsub(":|_", "", final_df$Gene)
  
  # Filter out rows with '0' combination
  final_df <- final_df %>% filter(!(Gene == "0" & Mutation == "0"))
  
  return(final_df)
}

# Function to calculate metrics
calculate_metrics <- function(df, threshold) {
  # Determine consensus mutations based on the threshold
  consensus_mutations <- df$mutation[rowSums(df[,-1]) >= threshold]
  
  # Initialize a metrics dataframe
  metrics <- data.frame(Method = names(df)[-1], TP = NA, FP = NA, FN = NA)
  
  # Calculate TP, FP, and FN for each method
  for (method in names(df)[-1]) {
    detected_mutations <- df$mutation[df[[method]] > 0]
    
    TP <- length(intersect(detected_mutations, consensus_mutations))
    FP <- length(setdiff(detected_mutations, consensus_mutations))
    FN <- length(setdiff(consensus_mutations, detected_mutations))
    
    metrics[metrics$Method == method, "TP"] <- TP
    metrics[metrics$Method == method, "FP"] <- FP
    metrics[metrics$Method == method, "FN"] <- FN
  }
  
  # Calculate Sensitivity and Precision
  metrics$Sensitivity <- metrics$TP / (metrics$TP + metrics$FN)
  metrics$Precision <- metrics$TP / (metrics$TP + metrics$FP)
  
  return(metrics)
}

# Example usage:
# Load the datasets
mutation_data_highConf <- read_excel("hybrid_subsitutions_real_highConf.xlsx", na = "NA")
mutation_data_midConf <- read_excel("hybrid_subsitutions_real_midConf.xlsx", na = "NA")

# Prepare the final_df for each dataset
M_final_df <- prepare_final_df(mutation_data_highConf)
H_final_df <- prepare_final_df(mutation_data_midConf)

# Set threshold and calculate metrics
threshold <- 4
metrics_highConf <- calculate_metrics(M_final_df, threshold)
metrics_midConf <- calculate_metrics(H_final_df, threshold)

print(metrics_highConf)
print(metrics_midConf)
