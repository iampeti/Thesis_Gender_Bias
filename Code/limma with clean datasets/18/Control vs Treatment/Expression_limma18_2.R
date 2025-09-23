library (limma)
library (matrixStats)


files = list.files(pattern="*.csv$") # Creates a list of strings containing all the names of the .csv files in our woring directory (factors).
a = read.table(list.files(pattern = "*.clean$"), header=TRUE, sep="\t") # Reads our dataset (should have .clean as extension).

a[a == "null"] = NA  # transform possible null values to NA


geneNames = a[,2] # create a variable containing only the GeneNames (IDENTIFIER)
bst = as.matrix(apply(a[,-c(1,2)], 2, as.numeric)) # create a matrix containing our dataset minus the first two columns (ID_REF and IDENTIFIER) applying
# as numeric to each column of a



#### Normalization

# Plotting the Histograms or boxplots of bst and log2(bst) we can decide if our data is normalized or not.
hist(bst)
hist(log2(bst))

boxplot(bst[sample(1:nrow(bst),1000),])
boxplot(log2(bst[sample(1:nrow(bst),1000),]))


### !! Only do this step if the data is not normalized beforehand !!!


bst = log2(bst)





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





group=paste(factors$disease_state,factors$gender,factors$age)



design=model.matrix(~0+group)
design
colnames(design)=c("con_f_a","con_f_mid","con_f_old","con_f_y","con_m_a","con_m_mid","con_m_old","con_m_y","sch_f_a","sch_f_mid","sch_f_old","sch_f_y","sch_m_a","sch_m_mid","sch_m_old","sch_m_y")
design


contrasts=makeContrasts( (con_f_a+con_f_mid+con_f_old+con_f_y+con_m_a+con_m_mid+con_m_old+con_m_y)-(sch_f_a+sch_f_mid+sch_f_old+sch_f_y+sch_m_a+sch_m_mid+sch_m_old+sch_m_y)
  ,                        levels = colnames(design))
contrasts



fit1 = lmFit(bst,design=design)
fit2 = contrasts.fit(fit1, contrasts = contrasts)
fit2 = eBayes(fit2)




res_coef1 = topTable(fit2, n= 10000)


l1= length(which(res_coef1$adj.P.Val < 0.05))
l1

if (l1 > 0) {
  write.table(geneNames[as.numeric(rownames(res_coef1)[1:l1])], "control vs treatment", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "contol vs treatment", quote=F, row.names=F, col.names=F)
}

expr_levels = data.frame(geneNames[as.numeric(rownames(res_coef1)[1:l1])],res_coef1$adj.P.Val[1:l1])
write.csv(expr_levels, "expr_levels18_2.csv", row.names = FALSE)
