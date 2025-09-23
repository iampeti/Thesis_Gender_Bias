library(ggplot2)
library(gridExtra)
library(grid)
library(formattable)

design_matrix <- data.frame(
  Sample = 1:6,
  B1 = c(1, 1, 1, 0, 0, 0),
  B2 = c(0, 0, 0, 1, 1, 1)
)

contrast_matrix <- data.frame(
  Parameter = c("B1", "B2"),
  C1 = c(1, -1),
  C2 = c(-1, 1)
)

grid.newpage()  # Create a new page for the plot
grid.table(design_matrix)  # Display the table
grid.text("Design Matrix", x = 0.5, y = 0.9, gp = gpar(fontsize = 16, fontface = "bold"))

grid.newpage()  # Create a new page for the plot
grid.table(contrast_matrix)  # Display the table
grid.text("Contrast Matrix", x = 0.5, y = 0.9, gp = gpar(fontsize = 16, fontface = "bold"))

