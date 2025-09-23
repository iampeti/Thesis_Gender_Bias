file_names = list.files(pattern="*.txt$")
file_names



for (i in 1:length(file_names)){
  name=file_names[i]
  y=read.table(file_names[i])
  assign(name,y)
}

all_results = list()

for (i in 1:length(file_names)){
  y=unlist(get(file_names[i]))
  all_results = c(all_results, y)
}

length(all_results)


counts = list()


for (item in all_results) {
  if (is.null(counts[[item]])) {
    counts[[item]] <- 1  # Initialize if first occurrence
  } else {
    counts[[item]] <- counts[[item]] + 1  # Increment if already exists
  }
}
counts

# Convert to a named vector for sorting
counts = unlist(counts)

sorted_counts <- sort(counts, decreasing = TRUE)

head(sorted_counts,578)


# Convert to a two-column matrix
result_matrix <- cbind(names(sorted_counts), as.numeric(sorted_counts))

# Convert to a data frame for better handling
result_df <- data.frame(Item = names(sorted_counts), Frequency = as.numeric(sorted_counts), row.names = NULL)

print(result_df)

write.csv(result_df, "output_frequency_correct.csv", row.names = FALSE)

