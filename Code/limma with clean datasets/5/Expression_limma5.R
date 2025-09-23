library (limma)
library (matrixStats)


files = list.files(pattern="*.csv$") # Creates a list of strings containing all the names of the .csv files in our woring directory (factors).
a = read.table(list.files(pattern = "*.clean$"), header=TRUE, sep="\t") # Reads our dataset (should have .clean as extension).

geneNames = a[,2] # create a variable containing only the GeneNames (IDENTIFIER)
bst = as.matrix(a[,-c(1,2)]) # create a matrix containing our dataset minus the first two columns (ID_REF and IDENTIFIER)



#### Normalization

# Plotting the Histograms or boxplots of bst and log2(bst) we can decide if our data is normalized or not.
hist(bst)
hist(log2(bst))

boxplot(bst[sample(1:nrow(bst),1000),])
boxplot(log2(bst[sample(1:nrow(bst),1000),]))


### !! Only do this step if the data is not normalized beforehand !!!


#bst = log2(bst)

#hist(bst)




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

group = paste(factors$gender,factors$disease_state,factors$tissue, sep='_')
group

design=model.matrix(~0+group)
design
colnames(design) = c("f_a_f_c","f_a_h","f_a_t_c","f_c_f_c","f_c_h","f_c_t_c","m_a_f_c","m_a_h","m_a_t_c","m_c_f_c","m_c_h","m_c_t_c")

design


contrasts=makeContrasts( ( (f_a_f_c + f_a_h + f_a_t_c)/3 - (f_c_f_c + f_c_h + f_c_t_c)/3 ) - ((m_a_f_c + m_a_h + m_a_t_c)/3 - (m_c_f_c + m_c_h + m_c_t_c)/3 ), levels = colnames(design))
contrasts



# (female control vs Female disease) vs (Male control VS Male disease)
# how much the disease  effect in males differs from the disease effect in females



fit1 = lmFit(bst,design=design)
fit2 = contrasts.fit(fit1, contrasts = contrasts)
fit2 = eBayes(fit2)



######## (female control vs Female Disease) vs (Male control VS Male Disease))

res1 = topTable(fit2, n= 10000)
res1


l1= length(which(res1$adj.P.Val < 0.05))
l1

if (l1 > 0) {
  write.table(geneNames[as.numeric(rownames(res1)[1:l1])], " (female control vs Female disease) vs (Male control VS Male disease)", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), " (female control vs Female disease) vs (Male control VS Male disease)", quote=F, row.names=F, col.names=F)
}



#######  with intercept, female control frontal cortex as reference level:

factors$disease_state = relevel(factors$disease,"non_Alzheimers")

design2 = model.matrix(~factors$gender*factors$disease_state*factors$tissue)
design2

colnames(design2)=c("intercept","Males","Alzheimers","Hippocampus","Temporal_Cortex","Males_Alzheimers","Males_Hippocampus","Males_Temporal_cortex",
                    "Alzheimers_hippocampus","Alzheimers_temporal_cortex","Male_Alzheimers_Hippocampus","Male_Alzheimers_Temporal_cortex")
design2

# intercept = female control frontal
# Males =  difference between Males and Females in control frontal
# Alzheimers = effect of Alzheimers in hippocampus in females
# Hippocampus =  effect of hippocampus in control Females
# Temporal_cortex = effect of temporal cortex in control females
# Males_Alzheimers =  captures how alzheimer effect differs from males compared to females in frontal cortex *
# Males_Hippocampus = captures how hippocampus effect differs from male to female in control
# Males_temporal_cortex = captures how temporal_cortex effect differs from male to female in control
# Alzheimers_Hippocampus = captures how alzheimers effect differs from frontal compared to hippocampus in females
# Alzheimers_temporal_cortex = captures how alzheimers effect differs from frontal compared to temporal in females
# Male_Alzheimers_Hippocampus = captures how alzheimers effect differs from males compared to females in hippocampus *
# Male_Alzheimers_Temporal_cortex = captures how alzheimers effect differs from males compared to females in temporal *



contrasts2 = makeContrasts(Males_Alzheimers, Male_Alzheimers_Hippocampus, Male_Alzheimers_Temporal_cortex,  levels = design2)
contrasts2


fit1 = lmFit(bst,design=design2)
fit2 = contrasts.fit(fit1, contrasts = contrasts2)
fit2 = eBayes(fit2)



######### coef = 1 : Males_Alzheimers =  captures how alzheimer effect differs from males compared to females in frontal cortex
######### coef = 2 : Male_Alzheimers_Hippocampus = captures how alzheimers effect differs from males compared to females in hippocampus
######### coef = 3 : Male_Alzheimers_Temporal_cortex = captures how alzheimers effect differs from males compared to females in temporal

res2_1 = topTable(fit2, coef = 1, n= 10000)

res2_1

l2= length(which(res2_1$adj.P.Val < 0.05))



if (l2 > 0) {
  write.table(geneNames[as.numeric(rownames(res2_1)[1:l2])], "with_intercept_female_vs_male_Alzheimers_frontal", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "with_intercept_female_vs_male_Alzheimers_frontal", quote=F, row.names=F, col.names=F)
  
}


res2_2 = topTable(fit2, coef = 2, n= 10000)
res2_2

l3= length(which(res2_2$adj.P.Val < 0.05))



if (l3 > 0) {
  write.table(geneNames[as.numeric(rownames(res2_2)[1:l3])], "with_intercept_female_vs_male_Alzheimers_hippocampus", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "with_intercept_female_vs_male_Alzheimers_hippocampus", quote=F, row.names=F, col.names=F)
}



res2_3 = topTable(fit2, coef = 3, n= 10000)
res2_3

l4= length(which(res2_3$adj.P.Val < 0.05))



if (l4 > 0) {
  write.table(geneNames[as.numeric(rownames(res2_3)[1:l4])], "with_intercept_female_vs_male_Alzheimers_temporal", quote=F, row.names=F, col.names=F)
} else {
  write.table(character(0), "with_intercept_female_vs_male_Alzheimers_temporal", quote=F, row.names=F, col.names=F)
}


