library (limma)
library (matrixStats)


files = list.files(pattern="*.csv$") # Creates a list of strings containing all the names of the .csv files in our woring directory (factors).
a = read.table(list.files(pattern = "*.clean$"), header=TRUE, sep="\t") # Reads our dataset (should have .clean as extension).

geneNames = a[,2] # create a variable containing only the GeneNames (IDENTIFIER)
bst = as.matrix(a[,-c(1,2)]) # create a matrix containing our dataset minus the first two columns (ID_REF and IDENTIFIER)



#### Normalization

# Plotting the Histograms of b and log2(b) we can decide if our data is normalized or not.
hist(bst)
hist(log2(bst))
boxplot(bst[sample(1:nrow(bst),1000),])

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

# Pavlidis: diavasa kapou oti gia raw data count kai limma kalutera to voom normalization pou ginetai sto telos. ? na to doume?


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







### Notes gia pavlidi:
# to modelo xwris kai to modelo me intercept vgazoun akrivws to idio apotelesma Yay! apla sto sugkekrimeno dataset exw parei P.Value < 0.01 giati me to 
# adjusted P Value < 0.05 den mou evgaze tpt kai ithela na tsekarw an trexoun swsta ta modela.
# To adjusted P value se kapoia datasets vgainei sxedon idia timh sta panw panw toulaxiston.
# voom normalization?


### Without intercept

group = paste(factors$gender,factors$smoking, sep='_')
group

design=model.matrix(~0+group)


############ Coefficients:

# coef = 1 : (female non smoking vs Female smoking) vs (Male non smoking VS Male smoking)
# how much the smoking  effect in males differs from the smoking effect in females
# coef = 2 : (female non smoking VS Male non smoking) vs (Female smoking vs Male smoking)
# how much the gender difference in non smoking differs from the gender difference in smoking  
# mathematicaly equivalent.


contrasts=makeContrasts((groupfemale_nonsmoker-groupfemale_smoker)-(groupmale_nonsmoker-groupmale_smoker),
                        (groupfemale_nonsmoker-groupmale_nonsmoker)-(groupfemale_smoker-groupmale_smoker), levels = colnames(design))
contrasts



fit1 = lmFit(bst,design=design)
fit2 = contrasts.fit(fit1, contrasts = contrasts)
fit2 = eBayes(fit2)



######## (female non smokers vs Female smokers) vs (Male non smokers VS Male smokers))

res_coef1 = topTable(fit2, coef = 1, n= 10000)
res_coef1


l1= length(which(res_coef1$P.Value < 0.01))


if (l1 > 0) {
  write.table(geneNames[as.numeric(rownames(res_coef1)[1:l1])], " (female non smokers vs Female smokers) vs (Male non smokers VS Male smokers)", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), " (female non smokers vs Female smokers) vs (Male non smokers VS Male smokers)", quote=F, row.names=F, col.names=F)
}


######## (female non smokers VS Male non smokers) vs (Female smokers vs Male smokers)

res_coef2 = topTable(fit2, coef = 2, n= 10000)
res_coef2


l2= length(which(res_coef2$P.Value < 0.01))
l2



if (l2 > 0) {
  write.table(geneNames[as.numeric(rownames(res_coef2)[1:l2])], "(female non smokers VS Male non smokers) vs (Female smokers vs Male smokers)", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "(female non smokers VS Male non smokers) vs (Female smokers vs Male smokers)", quote=F, row.names=F, col.names=F)
}



#######  with intercept, female non smokers as reference level:

design2 = model.matrix(~factors$gender*factors$smoking)
design2
colnames(design2)=c("intercept","genderM","smokingS","genderM_smokingS")
design2

# intercept = female non smoking
# genderM =  difference between males and females in non smoking group
# smokingS = effect of smoking in females
# genderM_smokingS = captures how smoking effect differs from males compared to females


contrasts2 = makeContrasts(genderM_smokingS,  levels = design2)
contrasts2


fit1 = lmFit(bst,design=design2)
fit2 = contrasts.fit(fit1, contrasts = contrasts2)
fit2 = eBayes(fit2)



######### captures how smoking effect differs from males compared to females

res2 = topTable(fit2, n= 10000)

l3= length(which(res2$P.Value < 0.01))



if (l3 > 0) {
  write.table(geneNames[as.numeric(rownames(res2)[1:l3])], "with_intercept_female_vs_male_smoking", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "with_intercept_female_vs_male_smoking", quote=F, row.names=F, col.names=F)
}


