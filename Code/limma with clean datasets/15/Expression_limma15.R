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

group = paste(factors$gender,factors$age, sep='_')
group

design=model.matrix(~0+group)
design



contrasts=makeContrasts( groupfemale_young - groupmale_young , 
                         groupfemale_middle - groupmale_middle ,
                         groupfemale_old - groupmale_old ,
                         groupfemale_ancient - groupmale_ancient,
                         levels = colnames(design))
contrasts


# coef = 1: young
# coef = 2: middle
# coef = 3: old
# coef = 4: ancient

fit1 = lmFit(bst,design=design)
fit2 = contrasts.fit(fit1, contrasts = contrasts)
fit2 = eBayes(fit2)




res_coef1 = topTable(fit2, coef = 1, n= 10000)


l1= length(which(res_coef1$adj.P.Val < 0.05))
l1

if (l1 > 0) {
  write.table(geneNames[as.numeric(rownames(res_coef1)[1:l1])], "female vs male young", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "female vs male young", quote=F, row.names=F, col.names=F)
}


expr_levels = data.frame(geneNames[as.numeric(rownames(res_coef1)[1:l1])],res_coef1$adj.P.Val[1:l1])
write.csv(expr_levels, "expr_levels15_1.csv", row.names = FALSE)

##########################################

res_coef2 = topTable(fit2, coef = 2, n= 10000)



l2= length(which(res_coef2$adj.P.Val < 0.05))
l2

if (l2 > 0) {
  write.table(geneNames[as.numeric(rownames(res_coef2)[1:l2])],  "female vs male middle", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0),  "female vs male middle", quote=F, row.names=F, col.names=F)
}


expr_levels2 = data.frame(geneNames[as.numeric(rownames(res_coef2)[1:l2])],res_coef2$adj.P.Val[1:l2])
write.csv(expr_levels2, "expr_levels15_2.csv", row.names = FALSE)

################################################

res_coef3 = topTable(fit2, coef = 3, n= 10000)



l3= length(which(res_coef3$adj.P.Val < 0.05))
l3

if (l3 > 0) {
  write.table(geneNames[as.numeric(rownames(res_coef3)[1:l3])],  "female vs male old", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0),  "female vs male old", quote=F, row.names=F, col.names=F)
}


expr_levels3 = data.frame(geneNames[as.numeric(rownames(res_coef3)[1:l3])],res_coef3$adj.P.Val[1:l3])
write.csv(expr_levels3, "expr_levels15_3.csv", row.names = FALSE)

###############################################################

res_coef4 = topTable(fit2, coef = 4, n= 10000)


l4= length(which(res_coef4$adj.P.Val < 0.05))
l4

if (l4 > 0) {
  write.table(geneNames[as.numeric(rownames(res_coef4)[1:l4])], "female vs male ancient", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0),  "female vs male ancient", quote=F, row.names=F, col.names=F)
}



#######  with intercept, female young as reference level:

factors$age = relevel(factors$age,"young")

design2 = model.matrix(~factors$gender*factors$age)
design2

colnames(design2)=c("intercept","Males","ancient","middle","Old","Males_ancient","Males_middle","Males_old")
design2

# intercept = female young
# Males =  effect of male in young
# ancient = effect of ancient in females
# middle =  effect of middle females
# old =  effect of old in Females
# Males_ancient = differences in males vs females in ancient
# Males_middle = differences in males vs females in middle
# Males_old = differences in males vs females in old



contrasts2 = makeContrasts(Males, Males+Males_middle, Males+Males_old, Males+Males_ancient , levels = design2)
contrasts2

# coef = 1: Males =  difference between Males and Females in young
# coef = 2: Males + Males_middle = differences in males vs females in middle
# coef = 3: Males + Males_old = differences in males vs females in old
# coef = 4: Males + Males_ancient = differences in males vs females in ancient


fit1 = lmFit(bst,design=design2)
fit2 = contrasts.fit(fit1, contrasts = contrasts2)
fit2 = eBayes(fit2)





res2_1 = topTable(fit2, coef = 1, n= 10000)

l2_1= length(which(res2_1$adj.P.Val < 0.05))
l2_1


if (l2_1 > 0) {
  write.table(geneNames[as.numeric(rownames(res2_1)[1:l2_1])], "with_intercept_female_vs_male_young", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "with_intercept_female_vs_male_young", quote=F, row.names=F, col.names=F)
}

################################################################

res2_2 = topTable(fit2, coef = 2, n= 10000)

l2_2= length(which(res2_2$adj.P.Val < 0.05))
l2_2


if (l2_2 > 0) {
  write.table(geneNames[as.numeric(rownames(res2_2)[1:l2_2])], "with_intercept_female_vs_male_middle", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "with_intercept_female_vs_male_middle", quote=F, row.names=F, col.names=F)
}


#######################################################

res2_3 = topTable(fit2, coef = 3, n= 10000)

l2_3 = length(which(res2_3$adj.P.Val < 0.05))
l2_3


if (l2_3 > 0) {
  write.table(geneNames[as.numeric(rownames(res2_3)[1:l2_3])], "with_intercept_female_vs_male_old", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "with_intercept_female_vs_male_old", quote=F, row.names=F, col.names=F)
}


#############################################

res2_4 = topTable(fit2, coef = 4, n= 10000)

l2_4 = length(which(res2_4$adj.P.Val < 0.05))
l2_4


if (l2_4 > 0) {
  write.table(geneNames[as.numeric(rownames(res2_4)[1:l2_4])], "with_intercept_female_vs_male_ancient", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "with_intercept_female_vs_male_ancient", quote=F, row.names=F, col.names=F)
}


