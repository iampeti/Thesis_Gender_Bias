inter=read.csv("intersection_all_DEGs.csv",header=FALSE, stringsAsFactors = FALSE)
inter
colnames(inter) = c("Gene")
inter

gene_counts <- as.data.frame(table(inter$Gene))

colnames(gene_counts) <- c("Gene", "Count")

gene_counts <- gene_counts[order(-gene_counts$Count), ]

write.csv(gene_counts, "intersection_all_DEGs_counted.csv", row.names = FALSE)

head(gene_counts)
