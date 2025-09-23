# Read .soft file in working directory
a = read.table("GDS3884.soft", header=TRUE, nrows=500, sep="\t")

# extract factor name
extract_name = function(x){
  l1 = sub("!subset_description = ","",x)
  l2 = gsub(" ","_",l1)
  return (l2)
}

# extract samples of each factor
extract_sample_id = function(x){
  l = sub("!subset_sample_id = ","",x)
  return (l)
}

# extract level of factor
extract_subset_type = function(x){
  l1 = sub("!subset_type = ","",x)
  l2 = gsub(" ","_",l1)
  return (l2)
}

k=dim(a)[1]

# write factor.level.csv files for each level containing the samples
for (i in 1:k){
  
  if (grepl("!subset_description", a[i,1])==TRUE){
    
    name1 = extract_name(a[i,1])
    name2 = extract_subset_type(a[i+2,1])
    name = paste0(name2,".",name1)
    name_csv = paste(name,".csv",sep="")
    writeLines(extract_sample_id(a[i+1,1]), name_csv)

  }
}
