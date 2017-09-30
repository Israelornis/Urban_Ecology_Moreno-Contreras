#####################################################
# Script to evaluate collinearity between variables #
#####################################################

### Script elaborated by Israel Moreno-Contreras (Ornithology, UNAM)

### Call libraries
library(corrplot) # A graphical display of a correlation matrix

### Read csv file of the clean variables
data <- read.csv("..\\data\\input_data\\clean_variables.csv", header = TRUE, row.names = 1)

### Obtain only predictor variables
predictor_variables <- as.matrix(data[, 10:15])

### Get correlation matrix
cor_variables <- cor(predictor_variables)

### Add p-values in a correlation matrix
# Perform function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# Save the table
p.mat <- cor.mtest(predictor_variables)

# Deal with missing values
diag(cor_variables) = NA

### Save corrplot in png format
tiff(filename="..\\data\\output_data\\corrplot.tif", compression = "lzw", units="mm", width=130, height=110, pointsize=9, res=1000)
corrplot(cor_variables, method = "ellipse", na.label = "NA", number.digits = TRUE, addrect=2, tl.col = "brown", bg = "beige", addgrid.col = "orange")
dev.off()

### Export correlation matrix in csv format
write.csv(cor_variables, file = "..\\data\\input_data\\corr_matrix.csv")

#####################################################################################