library (limma)
library (matrixStats)


files = list.files(pattern="*.csv$") # Creates a list of strings containing all the names of the .csv files in our woring directory (factors).
a = read.table(list.files(pattern = "*.clean$"), header=TRUE, sep="\t") # Reads our dataset (should have .clean as extension).

a[a == "null"] = NA  # transform possible null values to NA


geneNames = a[,2] # create a variable containing only the GeneNames (IDENTIFIER)
bst = as.matrix(apply(a[,-c(1,2)], 2, as.numeric)) # create a matrix containing our dataset minus the first two columns (ID_REF and IDENTIFIER) applying
# as numeric to each column of a

#### Normalization

# Plotting the Histograms of b and log2(b) we can decide if our data is normalized or not.
hist(bst)
hist(log2(bst))
boxplot(bst[sample(1:nrow(bst),1000),])
boxplot(log2(bst[sample(1:nrow(bst),1000),]))



### !! Only do this step if the data is not normalized beforehand !!!


#bst = log2(bst)





#Create variables named as the factor categories, which contain the GSMs.

for (i in 1:length(files)){
  name=sub(".csv","", files[i])
  y=as.character(read.table(files[i], sep=","))
  assign(name,y)
}

files=sub(".csv","", files)   #ommiting the .csv
files


######create factors

individuals = colnames(a[,-c(1,2)])


split_names = strsplit(files, "\\.") # split each element of files at .

factor_names = unique(sapply(split_names, `[`, 1))   # '[' is a subsetting operator and 1 specifies which part toextract. together: extract the first element

factors = list() # initialize factors as an empty list


for (factor in factor_names) {                                               # iterate through each factor                                  
  levels = sapply(split_names,'[',2)[sapply(split_names,'[',1) == factor]    ### extract levels for each factor (only apply the first sapply if the second returns true)
  
  
  factor_values = rep(NA, length(individuals))                              # Initialize a vector to store factor levels for each individual
  
  
  for (level in levels) {                                                   # iterate through each level of the factor 
    var_name = paste0(factor, ".", level)                                   # paste the factor and level we are currently iterating through
    gsm_list = get(var_name)                                                # get the value of pasted factor and level (GSM's)
    factor_values[individuals %in% gsm_list] = level                        # change from NA to factor level if the individual is in the gsm list for this level
  }
  
  
  
  factors[[factor]] = factor(factor_values, levels = levels)                # create each factor and store them in a list called factors (plural)
  
}





factors



#######  with intercept, female control omental as reference level:

#factors$disease = relevel(factors$disease,"control")

design2 = model.matrix(~factors$gender*factors$disease_state*factors$tissue)
design2

colnames(design2)=c("intercept","Males","Resistant","Subcutaneous","Males_Resistant","Males_Subcutaneous","Resistant_Subcutaneous","Males_Resistant_Subcutaneous")
design2

# intercept = female control omental
# Males =  effect of being male in control omental
# resistant = effect of resistant in omental females
# Subcutaneous = effect of subcutaneous in control females
# Males_resistant = Male vs Female in resistant omental * coef 1
# Males_Subcutaneous = Male vs Female in control in subcutaneous
# Resistant_Subcutaneous = control vs Resistant in females subcutaneous
# Male_Resistant_Resistant = male vs female in resistant and subcutaneous * coef 2


#resistant omental females:  intercept + resistant 
#resistant omental males: intercept + males + resistant + males_resistant
#resistant omental females - resistant omental males =  males + males_resistant

#############################
#resistant subcutaneous females: intercept + resistant + subcutaneous + resistant_subcutaneous
#resistant subcutaneous males: intercept + males + resistant + subcutaneous + males_resistant + males_subcutaneous + resistant_subcutaneos + males_resistant_subcutaneous
#afairesh: Males + Males_Resistant + Males_Subcutaneous + Males_Resistant_Subcutaneous



contrasts2 = makeContrasts(Males+Males_Resistant,Males + Males_Resistant + Males_Subcutaneous + Males_Resistant_Subcutaneous,  levels = design2)
contrasts2


fit1 = lmFit(bst,design=design2)
fit2 = contrasts.fit(fit1, contrasts = contrasts2)
fit2 = eBayes(fit2)



########################################

res2_coef1 = topTable(fit2, coef=1,n=Inf)

l2_1= length(which(res2_coef1$adj.P.Val < 0.05))

l2_1

res2_coef1
if (l2_1 > 0) {
  write.table(geneNames[as.numeric(rownames(res2_coef1)[1:l2_1])], "with intercept female vs male omental", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "with intercept female vs male omental", quote=F, row.names=F, col.names=F)
}

expr_levels = data.frame(geneNames[as.numeric(rownames(res2_coef1))],res2_coef1$logFC,res2_coef1$adj.P.Val)
write.csv(expr_levels, "expr_levels23_1.csv", row.names = FALSE)


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

# Rename columns for easier reference
expr_levels <- expr_levels %>%
  rename(Gene = geneNames.as.numeric.rownames.res2_coef1..., 
         logFC = res2_coef1.logFC, 
         adj.P.Val = res2_coef1.adj.P.Val)

# Identify the top DEG based on the lowest adjusted p-value
top_gene <- expr_levels %>% 
  arrange(adj.P.Val) %>% 
  slice_head(n = 10)  # Select the most significant gene

# Create volcano plot highlighting the top DEG
volcano_plot <- ggplot(expr_levels, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = adj.P.Val < 0.05), alpha = 0.5) +  # Color significant genes
  scale_color_manual(values = c("gray", "red")) +
  geom_point(data = top_gene, aes(x = logFC, y = -log10(adj.P.Val)), 
             color = "blue", size = 4) +  # Highlight top gene
  geom_text_repel(data = top_gene, aes(label = Gene), size = 5, box.padding = 0.5) +
  theme_minimal() +
  labs(title = "Volcano Plot Highlighting the Top DEG",
       x = "Log Fold Change",
       y = "-log10 Adjusted P-Value",
       color = "Significant")

# Print the plot
print(volcano_plot)



####################################

res2_coef2 = topTable(fit2, coef=2, n= 10000)

l2_2= length(which(res2_coef2$adj.P.Val < 0.05))

l2_2


if (l2_2 > 0) {
  write.table(geneNames[as.numeric(rownames(res2_coef2)[1:l2_2])], "with intercept female vs male subcutaneous", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "with intercept female vs male subcutaneous", quote=F, row.names=F, col.names=F)
}

expr_levels = data.frame(geneNames[as.numeric(rownames(res2_coef2)[1:l2_2])],res2_coef2$adj.P.Val[1:l2_2])
write.csv(expr_levels, "expr_levels23_2.csv", row.names = FALSE)



