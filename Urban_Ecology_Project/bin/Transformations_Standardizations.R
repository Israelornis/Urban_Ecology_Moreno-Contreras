#####################################################
# Script to transform and standirize variables      #
#####################################################

### Script elaborated by Israel Moreno-Contreras (Ornithology, UNAM)

### Call libraries
library(corrplot) # A graphical display of a correlation matrix
library(car)

### Read csv files of the Phylogenetic and Functional metric results
raw_data <- read.csv("..\\data\\input_data\\raw_variables.csv", header = TRUE, row.names = 1)

############################# Transformations ##################################
### Log transform Elevation 
ele <- log(raw_data$h_dem)

### Square root transform Annual Precipitation due it showed left-skew 
pre <- sqrt(raw_data$bio_12)

# Log transform the Annual Mean Temperature due it showed left-skew
tem <- log(raw_data$bio_1)

# Log transform Vegetation Complexity due it showed right-skew
veg <- log(raw_data$veg_com_Cl)

### Square root transform Human footprint and Alien bird species due 
### they showed right-skew
hfp <- log(raw_data$HFP)
abs <- log(raw_data$ABS)


### Obtain a table containing transformed predictor variables
transform_variables <- cbind(tem, pre, veg, ele, hfp, abs)

### Obtain a table containing transformed and z-standardised predictor variables
standa_transform_variables <- scale(transform_variables)

### Obtain a vector for Localities, richness, and phylogenetic and functional metrics,
### and x and y coordinates
Localities <- c("El Chamizal", "Club Campestre", "Parque Central",
                "Puerta Juarez", "Trepachanga", "Zaragoza", "Loma Blanca")
Habitat <- c("Urban green spaces", "Urban green spaces", "Urban green spaces",
             "Desert scrub", "Desert scrub", "Agricultural fields", "Agricultural fields")
Richness <- raw_data$Richness
ses.PD <- raw_data$ses.PD
NRI <- raw_data$NRI
NTI <- raw_data$NTI
ses.FD <- raw_data$ses.FD
FNRI <- raw_data$FNRI
FNTI <- raw_data$FNTI
x_coor <- raw_data$x_coor
y_coor <- raw_data$y_coor

### Obtain a table containing Localities, Habitat, richness, and phylogenetic and functional metrics
PD_FD_metrics <- cbind(Localities, Habitat, Richness, ses.PD, NRI, NTI, ses.FD, FNRI, FNTI)

### Obtain a table containing transformed and z-standardised predictor variables
### Localities, Habitat, richness, and phylogenetic and functional metrics, and coordinates
trans_stan_data <- cbind(PD_FD_metrics, standa_transform_variables, x_coor, y_coor)

### Export data frame containing "trans_stan_data"
write.csv(trans_stan_data, file = "..\\data\\input_data\\clean_variables.csv")