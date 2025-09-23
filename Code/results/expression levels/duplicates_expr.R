library(dplyr)

df <- read.csv("output_expression.csv")
length(df[,1])
df_unique <- df %>% distinct(df[[1]], .keep_all = TRUE)
df_unique <- df_unique[, 1:2]
length(df_unique[,1])
write.csv(df_unique, "output_expression.csv", row.names = FALSE)
