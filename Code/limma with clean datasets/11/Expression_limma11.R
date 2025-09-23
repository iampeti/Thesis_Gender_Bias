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

# OR if the .soft file has under value !dataset_value_type = transformed count it means it is normalized in some way? But ist probably
# safer to check every dataset seperatly with hist

# soft = read.table(list.files(pattern="*.soft$"), nrows=500, sep="\t") # reading the soft file

# for (i in 1:50) {                                                    # find the line that contains the dataset_value_type and save it
#   if (grepl("!dataset_value_type = ", as.character(soft[i,]))){
#    value_type = soft[i,]
#   }
# }


# if (value_type == "!dataset_value_type = count"){                  #normalize the data if its not normalized
#   bst = log2(bst)
# }




### !! Only do this step if the data is not normalized beforehand !!!

# Pavlidis: kapou diavasa kapou oti gia raw data count kai limma kalutera to voom normalization pou ginetai sto telos. ? na to doume?


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





factors


### Without intercept

#group = paste(factors$gender,factors$disease_state, sep='_')
#group

design=model.matrix(~0+factors$gender)
design
colnames(design)=c("female","male")





contrasts=makeContrasts(female-male,
                        levels = colnames(design))
contrasts



fit1 = lmFit(bst,design=design)
fit2 = contrasts.fit(fit1, contrasts = contrasts)
fit2 = eBayes(fit2)




res1 = topTable(fit2, n= 10000)
res1


l1= length(which(res1$adj.P.Val < 0.05))
l1

if (l1 > 0) {
  write.table(geneNames[as.numeric(rownames(res1)[1:l1])], "female_vs_male", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "female_vs_male", quote=F, row.names=F, col.names=F)
}

expr_levels = data.frame(geneNames[as.numeric(rownames(res1)[1:l1])],res1$adj.P.Val[1:l1])
write.csv(expr_levels, "expr_levels11.csv", row.names = FALSE)



#######  with intercept, female as reference level:



design2 = model.matrix(~factors$gender)
design2

colnames(design2)=c("intercept","Males")
design2

# intercept = female
# Males =  difference between Males and Females


contrasts2 = makeContrasts(Males,  levels = design2)
contrasts2


fit1 = lmFit(bst,design=design2)
fit2 = contrasts.fit(fit1, contrasts = contrasts2)
fit2 = eBayes(fit2)



 
res2 = topTable(fit2, n= 10000)

l2= length(which(res2$adj.P.Val < 0.05))



if (l2 > 0) {
  write.table(geneNames[as.numeric(rownames(res2)[1:l2])], "with_intercept_female_vs_male", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "with_intercept_female_vs_male", quote=F, row.names=F, col.names=F)
}



