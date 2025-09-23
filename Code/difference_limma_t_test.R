library(dplyr)

results_limma = read.table(file="diff_expressed_fm",header=FALSE)
results_t_test = read.table(file="results_t_test",header=FALSE)

a= (union(results_limma, results_t_test))
b= inner_join(results_limma, results_t_test)



difference = setdiff(a,b)

write.table(difference, "difference_limma_t_test", quote=F, row.names=F, col.names=F)