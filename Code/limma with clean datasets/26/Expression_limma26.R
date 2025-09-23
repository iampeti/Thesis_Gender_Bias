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


### Without intercept

group = paste(factors$gender,factors$individual, sep='_')
group

design=model.matrix(~0+group)
design



contrasts=makeContrasts(groupfemale_nomadic-groupmale_nomadic, groupfemale_rural-groupmale_rural, groupfemale_urban-groupmale_urban,
                        levels = colnames(design))
contrasts



fit1 = lmFit(bst,design=design)
fit2 = contrasts.fit(fit1, contrasts = contrasts)
fit2 = eBayes(fit2)



######## coef = 1 : female vs male in nomadic
######## coef = 2 : female vs male in rural
######## coef = 3 : female vs male in urban

res_coef1 = topTable(fit2, coef = 1, n= 10000)



l1= length(which(res_coef1$adj.P.Val < 0.05))
l1

if (l1 > 0) {
  write.table(geneNames[as.numeric(rownames(res_coef1)[1:l1])], "female vs male nomadic", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "female vs male nomadic", quote=F, row.names=F, col.names=F)
}

expr_levels = data.frame(geneNames[as.numeric(rownames(res_coef1)[1:l1])],res_coef1$adj.P.Val[1:l1])
write.csv(expr_levels, "expr_levels26_1.csv", row.names = FALSE)

#######################################

res_coef2 = topTable(fit2, coef = 2, n= 10000)



l2= length(which(res_coef2$adj.P.Val < 0.05))
l2



if (l2 > 0) {
  write.table(geneNames[as.numeric(rownames(res_coef2)[1:l2])], "female vs male rural",, quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "female vs male rural", quote=F, row.names=F, col.names=F)
}

expr_levels = data.frame(geneNames[as.numeric(rownames(res_coef2)[1:l2])],res_coef2$adj.P.Val[1:l2])
write.csv(expr_levels, "expr_levels26_2.csv", row.names = FALSE)

#######################################

res_coef3 = topTable(fit2, coef = 3, n= 10000)



l3= length(which(res_coef3$adj.P.Val < 0.05))
l3



if (l3 > 0) {
  write.table(geneNames[as.numeric(rownames(res_coef3)[1:l3])], "female vs male urban",, quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "female vs male urban", quote=F, row.names=F, col.names=F)
}

expr_levels = data.frame(geneNames[as.numeric(rownames(res_coef3)[1:l3])],res_coef3$adj.P.Val[1:l3])
write.csv(expr_levels, "expr_levels26_3.csv", row.names = FALSE)


#######  with intercept, female nomadic as reference level:

#factors$disease = relevel(factors$disease,"control")

design2 = model.matrix(~factors$gender*factors$individual)
design2

colnames(design2)=c("intercept","Males","Rural","Urban","Males_Rural","Males_Urban")
design2

# intercept = female nomadic
# Males =  effect of being male in nomadic
# Rural = effect of rural vs nomadic in females
# Urban = effect of urban vs nomadic in females
# Males_Rural = male in rural vs female rural
# Males_Urban = male in urban vs female urban


contrasts2 = makeContrasts(Males, Males_Rural+Males, Males_Urban+Males,  levels = design2)
contrasts2


fit1 = lmFit(bst,design=design2)
fit2 = contrasts.fit(fit1, contrasts = contrasts2)
fit2 = eBayes(fit2)



########################################

res2_coef1 = topTable(fit2, coef=1, n= 10000)

l2_1= length(which(res2_coef1$adj.P.Val < 0.05))

l2_1


if (l2_1 > 0) {
  write.table(geneNames[as.numeric(rownames(res2_coef1)[1:l2_1])], "with intercept female vs male nomadic", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "with intercept female vs male nomadic", quote=F, row.names=F, col.names=F)
}


####################################

res2_coef2 = topTable(fit2, coef=2, n= 10000)

l2_2= length(which(res2_coef2$adj.P.Val < 0.05))

l2_2


if (l2_2 > 0) {
  write.table(geneNames[as.numeric(rownames(res2_coef2)[1:l2_2])], "with intercept female vs male rural", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "with intercept female vs male rural", quote=F, row.names=F, col.names=F)
}

###############################

res2_coef3 = topTable(fit2, coef=3, n= 10000)

l2_3= length(which(res2_coef3$adj.P.Val < 0.05))

l2_3


if (l2_3 > 0) {
  write.table(geneNames[as.numeric(rownames(res2_coef3)[1:l2_3])], "with intercept female vs male urban", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "with intercept female vs male urban", quote=F, row.names=F, col.names=F)
}


