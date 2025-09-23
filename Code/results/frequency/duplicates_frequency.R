library(dplyr)

my_list <- readLines("r30.txt")

unique_list <- unique(my_list)

writeLines(unique_list, "r30.txt")
