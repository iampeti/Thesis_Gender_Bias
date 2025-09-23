txt_names = list.files(pattern="*.txt$")
CvT_1=readLines(" control vs gonadotrope")
CvT_2=readLines(" control vs null_cell")

CvT=intersect(CvT_1,CvT_2)

for (i in 1:length(txt_names)){
  name=txt_names[i]
  y=read.table(txt_names[i])
  assign(name,y)
}

all_results = list()

for (i in 1:length(txt_names)){
  y=unlist(get(txt_names[i]))
  all_results = c(all_results, y)
}

all_results
CvT
inter=intersect(all_results,CvT)
inter

write.table(inter, "intersection_9.txt", row.names = FALSE)
