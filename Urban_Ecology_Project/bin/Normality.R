#####################################################
# Script to evaluate normality of variables         #
#####################################################

### Script elaborated by Israel Moreno-Contreras (Ornithology, UNAM)

### Call libraries
library(maptools) # Allow read shapefiles
library(ggplot2)
library(Rmisc) # Merge plots
library(moments) # To quantify Kurtosis and Skewness
library(GGally) # Nice scatter plots
library(fuzzySim) # Evaluate multicollinearty between independent variables

### Read shapefile containing spatial variables
Points_variables <- readShapeSpatial("..\\data\\raw_data\\Puntos_utm.shp")

### Convert shapfile info into a data frame
Points_variables <- as.data.frame(Points_variables)
Points_variables_2 <- as.data.frame(Points_variables[,6:11])

### Obtain vector for x and y coordinates
x_coor <- Points_variables$coords.x1
y_coor <- Points_variables$coords.x2

### Read csv files of the Phylogenetic and Functional metric results
PD_results <- read.csv("..\\data\\input_data\\PD_metrics.results.csv", header = TRUE)

FD_results <- read.csv("..\\data\\input_data\\FD_metrics.results.csv", header = TRUE)

### Obtain data frame containing metrics (PD and FD) and independent variables
raw_data <- cbind(PD_results, FD_results[,4:7], Points_variables_2, x_coor, y_coor)

### Export data frame containing variables
write.csv(raw_data, file = "..\\data\\input_data\\raw_variables.csv")

### Read csv files of the Phylogenetic and Functional metric results
raw_data <- read.csv("..\\data\\input_data\\raw_variables.csv", header = TRUE, row.names = 1)

### Plot boxplot of the overall set
boxplot_richness <- ggplot(raw_data, aes(Richness, Richness)) + geom_boxplot() + labs(x = "Richness", y = "Number of species (S)") + xlab(expression(bold("Richness"))) + ylab(expression(bold("Number of species (S)"))) + theme(text = element_text(size=6.5), axis.text.x=element_blank())
boxplot_PD <- ggplot(raw_data, aes(PD, PD)) + geom_boxplot() + labs(x = "Phylogenetic Diversity", y = "Millions of years (Ma)") + xlab(expression(bold("Phylogenetic Diversity"))) + ylab(expression(bold("Millions of years (Ma)"))) + theme(text = element_text(size=6.5), axis.text.x=element_blank())
boxplot_ses.PD <- ggplot(raw_data, aes(ses.PD, ses.PD)) + geom_boxplot() + labs(x = "ses Phylogenetic Diversity", y = "") + xlab(expression(bold("ses Phylogenetic Diversity"))) + ylab(expression(bold(""))) + theme(text = element_text(size=6.5), axis.text.x=element_blank())
boxplot_NRI <- ggplot(raw_data, aes(NRI, NRI)) + geom_boxplot() + labs(x = "Net Relatedness Index", y = "") + xlab(expression(bold("Net Relatedness Index"))) + ylab(expression(bold(""))) + theme(text = element_text(size=6.5), axis.text.x=element_blank())
boxplot_NTI <- ggplot(raw_data, aes(NTI, NTI)) + geom_boxplot() + labs(x = "Nearest Taxon Index", y = "") + xlab(expression(bold("Nearest Taxon Index"))) + ylab(expression(bold(""))) + theme(text = element_text(size=6.5), axis.text.x=element_blank())
boxplot_FD <- ggplot(raw_data, aes(FD, FD)) + geom_boxplot() + labs(x = "Functional Diversity", y = "Functional Units") + theme(text = element_text(size=6.5), axis.text.x=element_blank())
boxplot_ses.FD <- ggplot(raw_data, aes(ses.FD, ses.FD)) + geom_boxplot() + labs(x = "ses Functional Diversity", y = "") + xlab(expression(bold("ses Functional Diversity"))) + ylab(expression(bold(""))) + theme(text = element_text(size=6.5), axis.text.x=element_blank())
boxplot_FNRI <- ggplot(raw_data, aes(FNRI, FNRI)) + geom_boxplot() + labs(x = "Functional Net Relatedness Index", y = "") + xlab(expression(bold("Functional Net Relatedness Index"))) + ylab(expression(bold(""))) + theme(text = element_text(size=6.5), axis.text.x=element_blank())
boxplot_FNTI <- ggplot(raw_data, aes(FNTI, FNTI)) + geom_boxplot() + labs(x = "Functional Nearest Taxon Index", y = "") + xlab(expression(bold("Functional Nearest Taxon Index"))) + ylab(expression(bold(""))) + theme(text = element_text(size=6.5), axis.text.x=element_blank())
boxplot_bio_1 <- ggplot(raw_data, aes(bio_1, bio_1)) + geom_boxplot() + labs(x = "Annual Mean Temperature", y = "°C") + xlab(expression(bold("bio 1"))) + ylab(expression(bold("°C"))) + theme(text = element_text(), axis.text.x=element_blank())
boxplot_bio_12 <- ggplot(raw_data, aes(bio_12, bio_12)) + geom_boxplot() + labs(x = "Annual Precipitation", y = "mm") + xlab(expression(bold("bio 12"))) + ylab(expression(bold("mm"))) + theme(text = element_text(), axis.text.x=element_blank())
boxplot_veg_com <- ggplot(raw_data, aes(veg_com_Cl, veg_com_Cl)) + geom_boxplot() + labs(x = "Vegetation Complexity", y = "") + xlab(expression(bold("Vegetation Complexity"))) + ylab(expression(bold(""))) + theme(text = element_text(), axis.text.x=element_blank())
boxplot_h_dem <- ggplot(raw_data, aes(h_dem, h_dem)) + geom_boxplot() + labs(x = "Elevation", y = "m") + xlab(expression(bold("Elevation"))) + ylab(expression(bold("m"))) + theme(text = element_text(), axis.text.x=element_blank())
boxplot_HFP <- ggplot(raw_data, aes(HFP, HFP)) + geom_boxplot() + labs(x = "Human footprint", y = expression(paste("HI per", km^{2}))) + xlab(expression(bold("HFP"))) + ylab(expression(bold(paste("HI per ", km^{2})))) + theme(text = element_text(), axis.text.x=element_blank())    
boxplot_ABS <- ggplot(raw_data, aes(ABS, ABS)) + geom_boxplot() + labs(x = "Alien bird species", y = "Number of species (S)") + xlab(expression(bold("ABS"))) + ylab(expression(bold("Number of species (S)"))) + theme(text = element_text(), axis.text.x=element_blank())      


############# Export dependent variables (boxplots) ################

### Save multiplot in png format, ONLY dependent variables
tiff(filename="..\\data\\output_data\\boxplot_multi_response.tif", compression = "lzw", units="mm", 
     width = 174, height = 110, pointsize=9, res=1000)

### Multiplot
multiplot(boxplot_ses.PD, boxplot_NRI, boxplot_NTI, 
          boxplot_ses.FD, boxplot_FNRI, boxplot_FNTI,boxplot_richness, cols = 4)

### Provide control over multiple graphics devices
dev.off()

####################################################################

############# Export independent variables (boxplots) ##############

### Save multiplot in png format, ONLY independent variables
tiff(filename="..\\data\\output_data\\boxplot_multi_explanatory.tif", compression = "lzw",
     units="mm", width = 174, height = 110, pointsize=9, res=1000)

### Multiplot
multiplot(boxplot_bio_1, boxplot_bio_12, boxplot_veg_com, boxplot_h_dem, boxplot_HFP, 
          boxplot_ABS, cols = 3)

### Provide control over multiple graphics devices
dev.off()

####################################################################

### Shapiro-Wilk Normality Test
shapiro.test(raw_data$pca_Clim_1)
shapiro.test(raw_data$veg_com_Cl)
shapiro.test(raw_data$PET)
shapiro.test(raw_data$cv)
shapiro.test(raw_data$h_dem)
shapiro.test(raw_data$density)
shapiro.test(raw_data$HFP)
shapiro.test(raw_data$ABS)

### Kurtosis
kurtosis(raw_data$bio_1)
kurtosis(raw_data$bio_12)
kurtosis(raw_data$veg_com_Cl)
kurtosis(raw_data$h_dem)
kurtosis(raw_data$HFP)
kurtosis(raw_data$ABS)

### Skewness
skewness(raw_data$bio_1)
skewness(raw_data$bio_12)
skewness(raw_data$veg_com_Cl)
skewness(raw_data$h_dem)
skewness(raw_data$HFP)
skewness(raw_data$ABS)

### Plots showned skewness and expect normal curve ###
# Elevation
skew_ele <- ggplot(raw_data, aes(x = h_dem), binwidth = 2) + 
  geom_histogram(bins = 7, aes(y = ..density..), fill = 'green', alpha = 0.2, color = "black") + 
  geom_density(colour = 'blue') + xlab(expression(bold('Elevation'))) + 
  ylab(expression(bold('Density'))) 

x <- seq(1100, 1350, length.out= 100)
df <- with(raw_data, data.frame(x = x, y = dnorm(x, mean(raw_data$h_dem), sd(raw_data$h_dem))))

skew_ele <- skew_ele + geom_line(data = df, aes(x = x, y = y), color = "red")

skew_ele # View plot

# Bio 1
skew_bio1 <- ggplot(raw_data, aes(x = bio_1), binwidth = 2) + 
  geom_histogram(bins = 7, aes(y = ..density..), fill = 'green', alpha = 0.2, color = "black") + 
  geom_density(colour = 'blue') + xlab(expression(bold('bio 1'))) + 
  ylab(expression(bold('Density'))) 

x <- seq(183, 195, length.out= 100)
df <- with(raw_data, data.frame(x = x, y = dnorm(x, mean(raw_data$bio_1), sd(raw_data$bio_1))))

skew_bio1 <- skew_bio1 + geom_line(data = df, aes(x = x, y = y), color = "red")

skew_bio1 # View plot

# Bio 12
skew_bio12 <- ggplot(raw_data, aes(x = bio_12), binwidth = 2) + 
  geom_histogram(bins = 7, aes(y = ..density..), fill = 'green', alpha = 0.2, color = "black") + 
  geom_density(colour = 'blue') + xlab(expression(bold('bio 12'))) + 
  ylab(expression(bold('Density'))) 

x <- seq(200, 270, length.out= 100)
df <- with(raw_data, data.frame(x = x, y = dnorm(x, mean(raw_data$bio_12), sd(raw_data$bio_12))))

skew_bio12 <- skew_bio12 + geom_line(data = df, aes(x = x, y = y), color = "red")

skew_bio12 # View plot

# Vegetation Complexity
skew_veg <- ggplot(raw_data, aes(x = veg_com_Cl), binwidth = 2) + 
  geom_histogram(bins = 7, aes(y = ..density..), fill = 'green', alpha = 0.2, color = "black") + 
  geom_density(colour = 'blue') + xlab(expression(bold('Vegetation Complexity'))) + 
  ylab(expression(bold('Density'))) 

x <- seq(21, 57, length.out= 100)
df <- with(raw_data, data.frame(x = x, y = dnorm(x, mean(raw_data$veg_com_Cl), sd(raw_data$veg_com_Cl))))

skew_veg <- skew_veg + geom_line(data = df, aes(x = x, y = y), color = "red")

skew_veg # View plot

# HFP
skew_hfp <- ggplot(raw_data, aes(x = HFP), binwidth = 2) + 
  geom_histogram(bins = 7, aes(y = ..density..), fill = 'green', alpha = 0.2, color = "black") + 
  geom_density(colour = 'blue', adjust = 1) + xlab(expression(bold('HFP'))) + 
  ylab(expression(bold('Density'))) 

x <- seq(53, 84, length.out= 100)
df <- with(raw_data, data.frame(x = x, y = dnorm(x, mean(raw_data$HFP), sd(raw_data$HFP))))

skew_hfp <- skew_hfp + geom_line(data = df, aes(x = x, y = y), color = "red")

skew_hfp # View plot

# ABS
skew_abs <- ggplot(raw_data, aes(x = ABS), binwidth = 2) + 
  geom_histogram(bins = 7, aes(y = ..density..), fill = 'green', alpha = 0.2, color = "black") + 
  geom_density(colour = 'blue', adjust = 1) + xlab(expression(bold('ABS'))) + 
  ylab(expression(bold('Density'))) 

x <- seq(1, 5, length.out= 100)
df <- with(raw_data, data.frame(x = x, y = dnorm(x, mean(raw_data$ABS), sd(raw_data$ABS))))

skew_abs <- skew_abs + geom_line(data = df, aes(x = x, y = y), color = "red")

skew_abs # View plot

############# Export independent variables (skewness plots) ##############

### Save multiplot in png format, ONLY independent variables
tiff(filename="..\\data\\output_data\\skewness_multi_explanatory.tif", compression = "lzw",
     units="mm", width=174, height=110, pointsize=9, res= 1000)
multiplot(skew_bio1, skew_bio12, skew_veg, skew_ele, skew_hfp, 
          skew_abs, cols = 3)
dev.off()

##########################################################################

### Plot Q-Q plots

####### Richness
# Find the slope and intercept of the line that passes through the 1st and 3rd
# quartile of the normal q-q plot
y     <- quantile(raw_data$Richness, c(0.25, 0.75)) # Find the 1st and 3rd quartiles
x     <- qnorm( c(0.25, 0.75))         # Find the matching normal values on the x-axis
slope <- diff(y) / diff(x)             # Compute the line slope
int   <- y[1] - slope * x[1]           # Compute the line intercept

qqplot_richness <- ggplot(raw_data, aes(sample = Richness)) + stat_qq() + labs(title="Richness") + geom_abline(intercept=int, slope=slope, col="blue") + theme(text = element_text(size=5, face = "bold"))


####### ses Phylogenetic Diversity
# Find the slope and intercept of the line that passes through the 1st and 3rd
# quartile of the normal q-q plot
y     <- quantile(raw_data$ses.PD, c(0.25, 0.75)) # Find the 1st and 3rd quartiles
x     <- qnorm( c(0.25, 0.75))         # Find the matching normal values on the x-axis
slope <- diff(y) / diff(x)             # Compute the line slope
int   <- y[1] - slope * x[1]           # Compute the line intercept

qqplot_ses.PD <- ggplot(raw_data, aes(sample = ses.PD)) + stat_qq() + labs(title="ses Phylogenetic Diversity") + geom_abline(intercept=int, slope=slope, col="blue") + theme(text = element_text(size=5, face = "bold"))


####### Net Relatedness Index
# Find the slope and intercept of the line that passes through the 1st and 3rd
# quartile of the normal q-q plot
y     <- quantile(raw_data$NRI, c(0.25, 0.75)) # Find the 1st and 3rd quartiles
x     <- qnorm( c(0.25, 0.75))         # Find the matching normal values on the x-axis
slope <- diff(y) / diff(x)             # Compute the line slope
int   <- y[1] - slope * x[1]           # Compute the line intercept

qqplot_NRI <- ggplot(raw_data, aes(sample = NRI)) + stat_qq() + labs(title="Net Relatedness Index") + geom_abline(intercept=int, slope=slope, col="blue") + theme(text = element_text(size=5, face = "bold"))


####### Nearest Taxon Index
# Find the slope and intercept of the line that passes through the 1st and 3rd
# quartile of the normal q-q plot
y     <- quantile(raw_data$NTI, c(0.25, 0.75)) # Find the 1st and 3rd quartiles
x     <- qnorm( c(0.25, 0.75))         # Find the matching normal values on the x-axis
slope <- diff(y) / diff(x)             # Compute the line slope
int   <- y[1] - slope * x[1]           # Compute the line intercept

qqplot_NTI <- ggplot(raw_data, aes(sample = NTI)) + stat_qq() + labs(title="Nearest Taxon Index") + geom_abline(intercept=int, slope=slope, col="blue") + theme(text = element_text(size=5, face = "bold"))


####### ses Functional Diversity
# Find the slope and intercept of the line that passes through the 1st and 3rd
# quartile of the normal q-q plot
y     <- quantile(raw_data$ses.FD, c(0.25, 0.75)) # Find the 1st and 3rd quartiles
x     <- qnorm( c(0.25, 0.75))         # Find the matching normal values on the x-axis
slope <- diff(y) / diff(x)             # Compute the line slope
int   <- y[1] - slope * x[1]           # Compute the line intercept

qqplot_ses.FD <- ggplot(raw_data, aes(sample = ses.FD)) + stat_qq() + labs(title="ses Functional Diversity") + geom_abline(intercept=int, slope=slope, col="blue") + theme(text = element_text(size=5, face = "bold"))


####### Functional Net Relatedness Index
# Find the slope and intercept of the line that passes through the 1st and 3rd
# quartile of the normal q-q plot
y     <- quantile(raw_data$FNRI, c(0.25, 0.75)) # Find the 1st and 3rd quartiles
x     <- qnorm( c(0.25, 0.75))         # Find the matching normal values on the x-axis
slope <- diff(y) / diff(x)             # Compute the line slope
int   <- y[1] - slope * x[1]           # Compute the line intercept

qqplot_FNRI <- ggplot(raw_data, aes(sample = FNRI)) + stat_qq() + labs(title="Functional Net Relatedness Index") + geom_abline(intercept=int, slope=slope, col="blue") + theme(text = element_text(size=5, face = "bold"))


####### Functional Nearest Taxon Index
# Find the slope and intercept of the line that passes through the 1st and 3rd
# quartile of the normal q-q plot
y     <- quantile(raw_data$FNTI, c(0.25, 0.75)) # Find the 1st and 3rd quartiles
x     <- qnorm( c(0.25, 0.75))         # Find the matching normal values on the x-axis
slope <- diff(y) / diff(x)             # Compute the line slope
int   <- y[1] - slope * x[1]           # Compute the line intercept

qqplot_FNTI <- ggplot(raw_data, aes(sample = FNTI)) + stat_qq() + labs(title="Functional Nearest Taxon Index") + geom_abline(intercept=int, slope=slope, col="blue") + theme(text = element_text(size=5, face = "bold"))


####### Annual Mean Temperature
# Find the slope and intercept of the line that passes through the 1st and 3rd
# quartile of the normal q-q plot
y     <- quantile(raw_data$bio_1, c(0.25, 0.75)) # Find the 1st and 3rd quartiles
x     <- qnorm( c(0.25, 0.75))         # Find the matching normal values on the x-axis
slope <- diff(y) / diff(x)             # Compute the line slope
int   <- y[1] - slope * x[1]           # Compute the line intercept

qqplot_bio_1 <- ggplot(raw_data, aes(sample = bio_1)) + stat_qq() + labs(title="bio 1") + geom_abline(intercept=int, slope=slope, col="blue") + theme(text = element_text(size=5.5, face = "bold"))


####### Annual Precipitation
# Find the slope and intercept of the line that passes through the 1st and 3rd
# quartile of the normal q-q plot
y     <- quantile(raw_data$bio_12, c(0.25, 0.75)) # Find the 1st and 3rd quartiles
x     <- qnorm( c(0.25, 0.75))         # Find the matching normal values on the x-axis
slope <- diff(y) / diff(x)             # Compute the line slope
int   <- y[1] - slope * x[1]           # Compute the line intercept

qqplot_bio_12 <- ggplot(raw_data, aes(sample = bio_12)) + stat_qq() + labs(title="bio 12") + geom_abline(intercept=int, slope=slope, col="blue") + theme(text = element_text(size=5.5, face = "bold"))


####### Vegetation Complexity
# Find the slope and intercept of the line that passes through the 1st and 3rd
# quartile of the normal q-q plot
y     <- quantile(raw_data$veg_com_Cl, c(0.25, 0.75)) # Find the 1st and 3rd quartiles
x     <- qnorm( c(0.25, 0.75))         # Find the matching normal values on the x-axis
slope <- diff(y) / diff(x)             # Compute the line slope
int   <- y[1] - slope * x[1]           # Compute the line intercept

qqplot_veg_com <- ggplot(raw_data, aes(sample = veg_com_Cl)) + stat_qq() + labs(title="Vegetation Complexity") + geom_abline(intercept=int, slope=slope, col="blue") + theme(text = element_text(size=5.5, face = "bold"))


####### Elevation
# Find the slope and intercept of the line that passes through the 1st and 3rd
# quartile of the normal q-q plot
y     <- quantile(raw_data$h_dem, c(0.25, 0.75)) # Find the 1st and 3rd quartiles
x     <- qnorm( c(0.25, 0.75))         # Find the matching normal values on the x-axis
slope <- diff(y) / diff(x)             # Compute the line slope
int   <- y[1] - slope * x[1]           # Compute the line intercept

qqplot_h_dem <-ggplot(raw_data, aes(sample = h_dem)) + stat_qq() + labs(title="Elevation") + geom_abline(intercept=int, slope=slope, col="blue") + theme(text = element_text(size=5.5, face = "bold"))


####### Human footprint
# Find the slope and intercept of the line that passes through the 1st and 3rd
# quartile of the normal q-q plot
y     <- quantile(raw_data$HFP, c(0.25, 0.75)) # Find the 1st and 3rd quartiles
x     <- qnorm( c(0.25, 0.75))         # Find the matching normal values on the x-axis
slope <- diff(y) / diff(x)             # Compute the line slope
int   <- y[1] - slope * x[1]           # Compute the line intercept

qqplot_HFP <- ggplot(raw_data, aes(sample = HFP)) + stat_qq() + labs(title="HFP") + geom_abline(intercept=int, slope=slope, col="blue") + theme(text = element_text(size=5.5, face = "bold"))


####### Alien bird species
# Find the slope and intercept of the line that passes through the 1st and 3rd
# quartile of the normal q-q plot
y     <- quantile(raw_data$ABS, c(0.25, 0.75)) # Find the 1st and 3rd quartiles
x     <- qnorm( c(0.25, 0.75))         # Find the matching normal values on the x-axis
slope <- diff(y) / diff(x)             # Compute the line slope
int   <- y[1] - slope * x[1]           # Compute the line intercept

qqplot_ABS <- ggplot(raw_data, aes(sample = ABS)) + stat_qq() + labs(title="ABS") + geom_abline(intercept=int, slope=slope, col="blue") + theme(text = element_text(size=5.5, face = "bold"))


############# Export dependent variables (qqplots) ###############

### Save qqplot (multiplot) in png format, ONLY dependent variables
tiff(filename="..\\data\\output_data\\qqplot_multi_response.tif", compression = "lzw", units="mm", width=174, height=110, pointsize=9, res=1000)

### Multiplot (qqplot), ONLY dependent variables
multiplot(qqplot_ses.PD, qqplot_NRI, qqplot_NTI, 
          qqplot_ses.FD, qqplot_FNRI, qqplot_FNTI, qqplot_richness, cols = 4)

### Provide control over multiple graphics devices
dev.off()

##################################################################

############# Export independent variables (qqplots) #############

### Save qqplot (multiplot) in png format, ONLY independent variables
tiff(filename="..\\data\\output_data\\qqplot_multi_explanatory.tif", compression = "lzw", units="mm", width=174, height=110, pointsize=9, res=1000)

### Multiplot
multiplot(qqplot_bio_1, qqplot_bio_12, qqplot_veg_com, qqplot_h_dem,
          qqplot_HFP, qqplot_ABS, cols = 3)

### Provide control over multiple graphics devices
dev.off()

##################################################################

### Scatterplots ###
sca_ses.PD_tem <- ggplot(raw_data, aes(ses.PD, bio_1, color = bio_1)) + geom_point(shape = 16, size = 5, show.legend = FALSE) + theme_minimal()
sca_ses.PD_pre <- ggplot(raw_data, aes(ses.PD, bio_12, color = bio_12)) + geom_point(shape = 16, size = 5, show.legend = FALSE) + theme_minimal()
sca_ses.PD_veg <- ggplot(raw_data, aes(ses.PD, veg_com_Cl, color = veg_com_Cl)) + geom_point(shape = 16, size = 5, show.legend = FALSE) + theme_minimal()
sca_ses.PD_ele <- ggplot(raw_data, aes(ses.PD, h_dem, color = h_dem)) + geom_point(shape = 16, size = 5, show.legend = FALSE) + theme_minimal()
sca_ses.PD_hfp <- ggplot(raw_data, aes(ses.PD, HFP, color = HFP)) + geom_point(shape = 16, size = 5, show.legend = FALSE) + theme_minimal()
sca_ses.PD_abs <- ggplot(raw_data, aes(ses.PD, ABS, color = ABS)) + geom_point(shape = 16, size = 5, show.legend = FALSE) + theme_minimal()

### Obtain data frames with response variable and explanatory variables
df_ses.PD <- raw_data[, c(5,12,13,14,15,16,17)]
df_NRI <- raw_data[, c(6,12,13,14,15,16,17)]
df_NTI <- raw_data[, c(7,12,13,14,15,16,17)]
df_ses.FD <- raw_data[, c(9,12,13,14,15,16,17)]
df_FNRI <- raw_data[, c(10,12,13,14,15,16,17)]
df_FNTI <- raw_data[, c(11,12,13,14,15,16,17)]
df_response <- raw_data[, c(3,5,6,7,9,10,11)]

### Fuction to edit scatterplots 
lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "blue", size = 1) +
    geom_smooth(method = method, color = "red", ...)
  p
}

### Save plots in tif format

tiff(filename="..\\data\\output_data\\scatter_ses.PD.tif", compression = "lzw", units="mm", width=174, height=110, pointsize=9, res=1000)
ggpairs(df_ses.PD, lower=list(continuous = wrap(lowerFn, method = "lm", size= 1))) + theme(text = element_text(size = 6.5))
dev.off()

tiff(filename="..\\data\\output_data\\scatter_NRI.tif", compression = "lzw", units="mm", width=174, height=110, pointsize=9, res=1000)
ggpairs(df_NRI, lower=list(continuous = wrap(lowerFn, method = "lm", size= 1))) + theme(text = element_text(size = 6.5))
dev.off()

tiff(filename="..\\data\\output_data\\scatter_NTI.tif", compression = "lzw", units="mm", width=174, height=110, pointsize=9, res=1000)
ggpairs(df_NTI, lower=list(continuous = wrap(lowerFn, method = "lm", size= 1))) + theme(text = element_text(size = 6.5))
dev.off()

tiff(filename="..\\data\\output_data\\scatter_ses.FD.tif", compression = "lzw", units="mm", width=174, height=110, pointsize=9, res=1000)
ggpairs(df_ses.FD, lower=list(continuous = wrap(lowerFn, method = "lm", size= 1))) + theme(text = element_text(size = 6.5))
dev.off()

tiff(filename="..\\data\\output_data\\scatter_FNRI.tif", compression = "lzw", units="mm", width=174, height=110, pointsize=9, res=1000)
ggpairs(df_FNRI, lower=list(continuous = wrap(lowerFn, method = "lm", size= 1))) + theme(text = element_text(size = 6.5))
dev.off()

tiff(filename="..\\data\\output_data\\scatter_FNTI.tif", compression = "lzw", units="mm", width=174, height=110, pointsize=9, res=1000)
ggpairs(df_FNTI, lower=list(continuous = wrap(lowerFn, method = "lm", size= 1))) + theme(text = element_text(size = 6.5))
dev.off()

tiff(filename="..\\data\\output_data\\scatter_response_variables.tif", compression = "lzw", units="mm", width=174, height=110, pointsize=9, res=1000)
ggpairs(df_response, lower=list(continuous = wrap(lowerFn, method = "lm", size= 1))) + theme(text = element_text(size = 6.5))
dev.off()