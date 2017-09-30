#####################################################
# Script to visualize boxplots of the results       #
#####################################################

### Script elaborated by Israel Moreno-Contreras (Ornithology, UNAM)

### Call libraries
library(ggplot2)
library(Rmisc) # Merge plots

results <- read.csv("..\\data\\input_data\\clean_variables.csv", header = TRUE)

pallete <- c("red",
             "blue",
             "green")

boxplot_ses.PD <- ggplot(results, aes(Habitat, ses.PD)) + geom_boxplot(fill = "white", colour = pallete) + theme(text = element_text(size=6.7), axis.text.x = element_text(colour="grey20",size=5.1, face="bold"), axis.text.y = element_text(colour="grey20",size=5.1, face="bold"), axis.title.y = element_text(colour="grey20",size=7,face="bold"), axis.title.x = element_text(colour="grey20",size=7,face="bold"))
boxplot_NTI <- ggplot(results, aes(Habitat, NTI)) + geom_boxplot(fill = "white", colour = pallete) + theme(text = element_text(size=6.7), axis.text.x = element_text(colour="grey20",size=5.1, face="bold"), axis.text.y = element_text(colour="grey20",size=5.1, face="bold"), axis.title.y = element_text(colour="grey20",size=7,face="bold"), axis.title.x = element_text(colour="grey20",size=7,face="bold"))
boxplot_NRI <- ggplot(results, aes(Habitat, NRI)) + geom_boxplot(fill = "white", colour = pallete) + theme(text = element_text(size=6.7), axis.text.x = element_text(colour="grey20",size=5.1, face="bold"), axis.text.y = element_text(colour="grey20",size=5.1, face="bold"), axis.title.y = element_text(colour="grey20",size=7,face="bold"), axis.title.x = element_text(colour="grey20",size=7,face="bold"))
boxplot_ses.FD <- ggplot(results, aes(Habitat, ses.FD)) + geom_boxplot(fill = "white", colour = pallete) + theme(text = element_text(size=6.7), axis.text.x = element_text(colour="grey20",size=5.1, face="bold"), axis.text.y = element_text(colour="grey20",size=5.1, face="bold"), axis.title.y = element_text(colour="grey20",size=7,face="bold"), axis.title.x = element_text(colour="grey20",size=7,face="bold"))
boxplot_FNTI <- ggplot(results, aes(Habitat, FNTI)) + geom_boxplot(fill = "white", colour = pallete) + theme(text = element_text(size=6.7), axis.text.x = element_text(colour="grey20",size=5.1, face="bold"), axis.text.y = element_text(colour="grey20",size=5.1, face="bold"), axis.title.y = element_text(colour="grey20",size=7,face="bold"), axis.title.x = element_text(colour="grey20",size=7,face="bold"))
boxplot_FNRI <- ggplot(results, aes(Habitat, FNRI)) + geom_boxplot(fill = "white", colour = pallete) + theme(text = element_text(size=6.7), axis.text.x = element_text(colour="grey20",size=5.1, face="bold"), axis.text.y = element_text(colour="grey20",size=5.1, face="bold"), axis.title.y = element_text(colour="grey20",size=7,face="bold"), axis.title.x = element_text(colour="grey20",size=7,face="bold"))


### Save multiplot in png format, ONLY dependent variables
tiff(filename="..\\data\\output_data\\boxplot_results.tif", compression = "lzw", units="mm", 
     width = 174, height = 110, pointsize=2, res=1000)

### Multiplot
multiplot(boxplot_ses.PD, boxplot_ses.FD, boxplot_NRI,  boxplot_FNRI, 
          boxplot_NTI, boxplot_FNTI, cols = 3)

### Provide control over multiple graphics devices
dev.off()






