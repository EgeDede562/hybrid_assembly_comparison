library(readxl)
library(ggplot2)

# Load the datasets
mutation_data_highConf <- read_excel("hybrid_subsitutions_real_highConf.xlsx", na = "NA")
mutation_data_midConf <- read_excel("hybrid_subsitutions_real_midConf.xlsx", na = "NA")
row.names(mutation_data_highConf)<-mutation_data_highConf$...1
row.names(mutation_data_midConf)<-mutation_data_midConf$...1
mutation_data_highConf$...1<-NULL
mutation_data_midConf$...1<-NULL

generate_plot <- function(mutation_data) {
  
  # Extract method names from the dataset
  method_names <- colnames(mutation_data)
  mutations_all <- list()
  
  # Extract mutations for each method
  for (method in method_names) {
    # Split mutations by comma and unlist them
    mutations_all[[method]] <- unlist(strsplit(unlist(mutation_data[, method]), split = ","))
  }
  
  # Filter out NA and empty values from the mutations
  mutations_all <- lapply(mutations_all, function(x) x[!is.na(x) & x != ""])
  
  # Count the occurrence of each mutation for every method
  mutations_counts <- lapply(mutations_all, function(x) table(x))
  
  # Extract unique mutations for each method (mutations that appear only once)
  unique_mutations <- lapply(mutations_counts, function(x) names(x[x == 1]))
  
  # Combine mutations from all methods
  all_mutations <- unlist(mutations_all)
  
  # Filter out empty strings
  all_mutations <- all_mutations[all_mutations != ""]
  
  # Compute the global occurrence of each mutation across all methods
  all_mutations_count <- table(all_mutations)
  
  # Categorize mutations based on their occurrence across methods
  mutation_categories <- list(
    Unique = names(all_mutations_count[all_mutations_count == 1]),
    Universal = names(all_mutations_count[all_mutations_count == 8]),
    High_commonality = names(all_mutations_count[all_mutations_count %in% c(6,7)]),
    Medium_commonality = names(all_mutations_count[all_mutations_count %in% c(4,5)]),
    Low_commonality = names(all_mutations_count[all_mutations_count %in% c(2,3)])
  )
  
  # Count the number of mutations for each category and method
  category_counts <- matrix(nrow = length(method_names), 
                            ncol = length(names(mutation_categories)), 
                            dimnames = list(method_names, names(mutation_categories)))
  
  for (method in method_names) {
    for (category in names(mutation_categories)) {
      category_counts[method, category] <- sum(mutations_all[[method]] %in% mutation_categories[[category]])
    }
  }
  
  # Convert matrix to dataframe for easier visualization
  df <- as.data.frame(category_counts)
  
  # Convert dataframe to long format for ggplot2
  df$Method <- row.names(df)
  df_long <- tidyr::pivot_longer(df, cols = -Method, names_to = "Category", values_to = "Count")
  
  # 1. Reorder Methods
  desired_order <- c("ilu_ori", "ilu_alt", "ilu_alt_2", "ont_ori", "ont_alt", "ont_alt_2", "ilu_dnv", "hybrid")
  df_long$Method <- factor(df_long$Method, levels = desired_order)
  
  # 3. Reorder Stacking
  category_order <- c("Universal", "High_commonality", "Medium_commonality", "Low_commonality", "Unique")
  df_long$Category <- factor(df_long$Category, levels = category_order)
  
  # Plot with custom colors
  custom_colors <- c("Universal" = "#6c757d", "Low_commonality" = "#f0e68c", "Medium_commonality" = "#98fb10", "High_commonality" = "#dda0dd", "Unique" = "#ff7f50")
  
  ggplot(df_long, aes(x = Method, y = Count, fill = Category)) + 
    geom_bar(stat = "identity", position = "stack") + 
    scale_fill_manual(values = custom_colors) +
    labs(title = "Mutation Categories per Method", x = "Method", y = "Number of Mutations") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
 return(ggplot(df_long, aes(x = Method, y = Count, fill = Category)) + 
           geom_bar(stat = "identity", position = "stack") + 
           labs(title = "Mutation Categories per Method", x = "Method", y = "Number of Mutations") +
           theme_minimal() +
           theme(axis.text.x = element_text(angle = 45, hjust = 1)))
}

unique_counts_highConf <- generate_plot(mutation_data_highConf)
unique_counts_midConf <- generate_plot(mutation_data_midConf)

print(unique_counts_highConf)
print(unique_counts_midConf)
##GRAPH PART 1 FINISHED
