###DESCRIPTIVE STATISTICS PART I###
##This code snippet is used to calculate and graph the success rate for the eight different methods##

#LOAD THE NECESSARY LIBRARIES
# Ensure necessary packages are installed and loaded
packages <- c("readxl", "ggplot2")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(packages, library, character.only = TRUE)


#IMPORT THE DATASET AND PREPARE THE DATASET
# Check if the Excel file exists to avoid read errors
file_path <- "01_analysis_table_real.xlsx"
if(!file.exists(file_path)) {
  stop("File does not exist: ", file_path)
}

# Import the dataset
df <- read_excel(file_path, na = "NA", col_names = TRUE)

# Convert the tibble to a standard data frame
df <- as.data.frame(df)

# Set the first column as row names and remove it from the dataframe
rownames(df) <- df$Sample
df[,1] <- NULL

# Identify columns tht need to be converted and convert them
df <- df %>%
  mutate(across(contains("nextclade_qval_"), ~as.numeric(as.character(.x))))

# Backup the original dataframe
dfbu <- df

# CALCULATE THE SUCCESS RATE FOR METHODS
# Define the columns that contain the genome fraction columns
gf_cols <- colnames(df)[15:22]

# Ensure the specified columns exist in the dataframe
missing_cols <- setdiff(gf_cols, colnames(df))
if(length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

# Count of non-NA values in genome fraction columns
row1 <- sapply(df[gf_cols], function(x) sum(!is.na(x)))

# Count of genomic fractions greater than 0.90
row2 <- sapply(df[gf_cols], function(x) sum(x > 0.9, na.rm = TRUE))

# Calculate the percentage of genomes > 0.90
row3 <- row2 / row1 * 100

# Combines the calculated rows into a results  matrix
results_matrix <- rbind(row1, row2, row3)

# Output the results matrix
print(results_matrix)

# GENERATE A BAR PLOT FOR SUCCESS RATES
# Converts the results matrix to a dataframe for ggplot
df_plot <- data.frame(
  Method = colnames(results_matrix),
  Total = results_matrix[1, ],
  Over_90 = results_matrix[2, ]
)

# Creates a bar plot using ggplot2
ggplot(df_plot, aes(x = Method, y = Total, fill = "Total")) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_bar(aes(y = Over_90, fill = ">90% Genomic Fraction"), stat = "identity", position = "dodge") +
  labs(title = "Number of Consensus Sequences by Method", y = "Number of Sequences") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Total" = "#a6cee3", ">90% Genomic Fraction" = "#1f78b4"))

#END OF DESCRIPTIVE STATISTICS PART I#



