file_names = list.files(pattern="*.csv$")
file_names



for (i in 1:length(file_names)){
  name=file_names[i]
  y=as.matrix(read.table(file_names[i],header=TRUE))
  assign(name,y)
}

all_genes <- c()
all_pvals <- c()

# Loop through each file and extract the gene names and adjusted p-values
for (file in file_names) {
  # Read the CSV file into a matrix
  data <- as.matrix(read.table(file, header = TRUE, stringsAsFactors = FALSE))
  
  # Extract gene names and adjusted p-values (removing potential commas in p-values)
  genes <- data[, 1]  # First column contains gene names
  pvals <- as.numeric(sub(",", "", data[, 2]))  # Second column contains adjusted p-values
  
  # Append to the overall vectors
  all_genes <- c(all_genes, genes)
  all_pvals <- c(all_pvals, pvals)
}

# Combine into a matrix
result_matrix <- cbind(all_genes, all_pvals)



# Convert to a data frame for easier handling
result_df <- data.frame(Gene = all_genes, Adjusted_PValue = all_pvals, stringsAsFactors = FALSE)



result_df$Adjusted_PValue <- as.numeric(result_df$Adjusted_PValue)


#Sort the data frame by Adjusted_PValue in increasing order
sorted_result_df <- result_df[order(result_df$Adjusted_PValue), ]

# Print the first few rows to check
head(sorted_result_df)

write.csv(sorted_result_df, "output_expression.csv", row.names = FALSE)



library(ggplot2)
library(dplyr)
deg_results = sorted_result_df

# Assuming you have a data frame `deg_results` with columns:
# - adj.P.Val (adjusted p-value)
# - Gene (gene name)

# Select the top genes by lowest adjusted p-value
top_n <- 10  # Number of top genes to highlight
top_genes <- deg_results %>% 
  arrange(Adjusted_PValue) %>% 
  slice_head(n = top_n)

print(top_genes)

# Create a bar plot
bar_plot <- ggplot(top_genes, aes(x = reorder(Gene, Adjusted_PValue), y = -log10(Adjusted_PValue))) +
  geom_col(fill = "steelblue") +
  coord_flip() +  # Flip for better readability
  theme_minimal() +
  labs(title = "Top Differentially Expressed Genes",
       x = "Gene",
       y = "-log10 Adjusted P-Value")

# Print the plot
print(bar_plot)
