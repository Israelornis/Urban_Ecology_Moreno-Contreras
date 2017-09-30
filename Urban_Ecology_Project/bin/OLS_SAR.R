#####################################################
# Script to perform OLS and SAR models              #
#####################################################

### Script elaborated by Israel Moreno-Contreras (Ornithology, UNAM)

### Call libraries
library(spdep)
library(ggplot2) # Get Moran and Bubble plots
library(Rmisc) # Merge plots
library(ggmap) # for fortifying shapefiles
library(maptools) # Allow read shapefiles

### Read csv file of the clean variables
clean_data <- read.csv("..\\data\\input_data\\clean_variables.csv", header = TRUE, row.names = 1)

### Get coordinates of the clean data
x_coor <- clean_data$x_coor
y_coor <- clean_data$y_coor

### Merge x and y coordinates
coords <- cbind(x_coor, y_coor)

# Convert coordinates into a distance matrix
coor.dists <- as.matrix(dist(cbind(clean_data$x_coor, clean_data$y_coor))) # Max = 24876.566 Min = 4407.763

#The k = 1 object is also useful in finding the minimum distance
# at which all areas have a distance-based neighbour. Using the nbdists function, we
# can calculate a list of vectors of distances corresponding to the neighbour object, here
# for first nearest neighbours. The greatest value will be the minimum distance needed
# to make sure that all the areas are linked to at least one neighbour
nb.1NN <- knn2nb(knearneigh(coords, k = 1))

nb.1NN

dsts <- unlist(nbdists(nb.1NN, coords))
summary(dsts) # Minimum distance employed in the upper bound is 16230 (17000 m)

### Extract max distance
max_1nn <- max(dsts)
max_1nn

### Obtain a Neighbourhood contiguity by distance
nb_0_17000 <- dnearneigh(coords, 0, 17000)

### Plot Neighbourhood contiguity by distance
plot(nb_0_17000, coords)

### Get spatial weights for neighbours lists
listw <- nb2listw(nb_0_17000, style = "B")
listw

#####################################################################################
#####################################################################################
#####################################################################################

############################# Richness ##############################################

### Model for the Richness (MAM)
ols_Richness <- lm(Richness ~ tem + ele + abs, data = clean_data)
summary(ols_Richness)

### Get residuals
residuals_Richness <- residuals.lm(ols_Richness)

### Moran's I test for spatial autocorrelation in residuals of Richness employing spdep library
moran.mc(residuals_Richness, listw, nsim = 999)
moran.plot(residuals_Richness, listw)

### Spatial lag of a numeric vector
Wsids_Richness = lag.listw(listw, residuals_Richness)

### Get data frame containing residuales and spatial lag
dat_Richness = data.frame(Resid_Richness = residuals_Richness, Wsids = Wsids_Richness)

### Plot Moran test using ggplot2
moran.plot_Richness <- ggplot(dat_Richness, aes(x = Resid_Richness, y = Wsids)) + geom_point() + geom_smooth(method = "lm") + xlab("Richness residuals") + ylab("Spatial lag of Richness residuals") + geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed') + theme(text = element_text(size=5.5))

#####################################################################################

############################# ses.PD ################################################

### Model for the ses.PD (MAM)
ols_ses.PD <- lm(ses.PD ~ ele + abs, data = clean_data) 
summary(ols_ses.PD)

### Get residuals
residuals_ses.PD <- residuals.lm(ols_ses.PD)

### Moran's I test for spatial autocorrelation in residuals of ses.PD employing spdep library
moran.mc(residuals_ses.PD, listw, nsim = 999)
moran.plot(residuals_ses.PD, listw)

### Spatial lag of a numeric vector
Wsids_ses.PD = lag.listw(listw, residuals_ses.PD)

### Get data frame containing residuales and spatial lag
dat_ses.PD = data.frame(Resid_ses.PD = residuals_ses.PD, Wsids = Wsids_ses.PD)

### Plot Moran test using ggplot2
moran.plot_ses.PD <- ggplot(dat_ses.PD, aes(x = Resid_ses.PD, y = Wsids)) + geom_point() + geom_smooth(method = "lm") + xlab("ses Phylogenetic Diversity residuals") + ylab("Spatial lag of ses PD residuals") + geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed') + theme(text = element_text(size=5.5))

#####################################################################################

############################# NRI ###################################################

### Model for the NRI (MAM)
ols_NRI <- lm(NRI ~ tem + ele + hfp + abs, data = clean_data)

### Get residuals
residuals_NRI <- residuals.lm(ols_NRI)

### Moran's I test for spatial autocorrelation in residuals of NRI employing spdep library
moran.test(residuals_NRI, listw)

moran.plot(residuals_NRI, listw)

### Spatial lag of a numeric vector
Wsids_NRI = lag.listw(listw, residuals_NRI)

### Get data frame containing residuales and spatial lag
dat_NRI = data.frame(Resid_NRI = residuals_NRI, Wsids = Wsids_NRI)

### Plot Moran test using ggplot2
moran.plot_NRI <- ggplot(dat_NRI, aes(x = Resid_NRI, y = Wsids)) + geom_point() + geom_smooth(method = "lm") + xlab("NRI residuals") + ylab("Spatial lag of NRI residuals") + geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed') + theme(text = element_text(size=5.5))

#####################################################################################

############################# NTI ###################################################

### Model for the NTI (MAM)
ols_NTI <- lm(NTI ~ veg + ele + hfp + abs, data = clean_data)

### Get residuals
residuals_NTI <- residuals.lm(ols_NTI)

### Moran's I test for spatial autocorrelation in residuals of NTI employing spdep library
moran.test(residuals_NTI, listw)

moran.plot(residuals_NTI, listw)

### Spatial lag of a numeric vector
Wsids_NTI = lag.listw(listw, residuals_NTI)

### Get data frame containing residuales and spatial lag
dat_NTI = data.frame(Resid_NTI = residuals_NTI, Wsids = Wsids_NTI)

### Plot Moran test using ggplot2
moran.plot_NTI <- ggplot(dat_NTI, aes(x = Resid_NTI, y = Wsids)) + geom_point() + geom_smooth(method = "lm") + xlab("NTI residuals") + ylab("Spatial lag of NTI residuals") + geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed') + theme(text = element_text(size=5.5))

#####################################################################################

############################# ses.FD ################################################

### Model for the ses.FD (MAM)
ols_ses.FD <- lm(ses.FD ~ veg + ele + hfp + abs, data = clean_data) 
summary(ols_ses.FD)

### Get residuals
residuals_ses.FD <- residuals.lm(ols_ses.FD)

### Moran's I test for spatial autocorrelation in residuals of ses.PD employing spdep library
moran.test(residuals_ses.FD, listw)

moran.plot(residuals_ses.FD, listw)

### Spatial lag of a numeric vector
Wsids_ses.FD = lag.listw(listw, residuals_ses.FD)

### Get data frame containing residuales and spatial lag
dat_ses.FD = data.frame(Resid_ses.FD = residuals_ses.FD, Wsids = Wsids_ses.FD)

### Plot Moran test using ggplot2
moran.plot_ses.FD <- ggplot(dat_ses.FD, aes(x = Resid_ses.FD, y = Wsids)) + geom_point() + geom_smooth(method = "lm") + xlab("ses Functional Diversity residuals") + ylab("Spatial lag of ses FD residuals") + geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed') + theme(text = element_text(size=5.5))

#####################################################################################

############################# FNRI ##################################################

### Model for the FNRI (MAM)
ols_FNRI <- lm(FNRI ~ ele + abs, data = clean_data) 

### Get residuals
residuals_FNRI <- residuals.lm(ols_FNRI)

### Moran's I test for spatial autocorrelation in residuals of FNRI employing spdep library
moran.test(residuals_FNRI, listw)

moran.plot(residuals_FNRI, listw)

### Spatial lag of a numeric vector
Wsids_FNRI = lag.listw(listw, residuals_FNRI)

### Get data frame containing residuales and spatial lag
dat_FNRI = data.frame(Resid_FNRI = residuals_FNRI, Wsids = Wsids_FNRI)

### Plot Moran test using ggplot2
moran.plot_FNRI <- ggplot(dat_FNRI, aes(x = Resid_FNRI, y = Wsids)) + geom_point() + geom_smooth(method = "lm") + xlab("FNRI residuals") + ylab("Spatial lag of FNRI residuals") + geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed') + theme(text = element_text(size=5.5))

#####################################################################################

############################# FNTI ##################################################

### Model for the FNTI (MAM)
ols_FNTI <- lm(FNTI ~ tem + ele + hfp + abs, data = clean_data) 

### Get residuals
residuals_FNTI <- residuals.lm(ols_FNTI)

### Moran's I test for spatial autocorrelation in residuals of FNTI employing spdep library
moran.test(residuals_FNTI, listw)

moran.plot(residuals_FNTI, listw)

### Spatial lag of a numeric vector
Wsids_FNTI = lag.listw(listw, residuals_FNTI)

### Get data frame containing residuales and spatial lag
dat_FNTI = data.frame(Resid_FNTI = residuals_FNTI, Wsids = Wsids_FNTI)

### Plot Moran test using ggplot2
moran.plot_FNTI <- ggplot(dat_FNTI, aes(x = Resid_FNTI, y = Wsids)) + geom_point() + geom_smooth(method = "lm") + xlab("FNTI residuals") + ylab("Spatial lag of FNTI residuals") + geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed') + theme(text = element_text(size=5.5))

#####################################################################################

### Export Multi Moran plot for OLS models

### Save Multi Moran plot in tif format
tiff(filename="..\\data\\output_data\\OLS_Moran_plot.tif", compression = "lzw", units="mm", width=174, height=110, pointsize=9, res=1000)
multiplot(moran.plot_ses.PD, moran.plot_ses.FD, moran.plot_NRI, moran.plot_FNRI, moran.plot_NTI, moran.plot_FNTI, moran.plot_Richness, cols = 4)
dev.off()

#####################################################################################
#####################################################################################
#####################################################################################


##### SAR models #####

### Spatial simultaneous autoregressive error model estimation for Richness 

spatial.error_Richness <- errorsarlm(Richness ~ tem + ele + abs, data = clean_data, listw)
summary(spatial.error_Richness, Nagelkerke = TRUE)

### Extract SAR Richness residuals
residuals_sar_Richness <- residuals.sarlm(spatial.error_Richness)
residuals_sar_Richness

### Moran's I test for spatial autocorrelation in residuals of Richness employing spdep library
moran.mc(residuals_sar_Richness, listw, nsim = 999)
moran.plot(residuals_sar_Richness, listw)

### Spatial lag of a numeric vector
Wsids_sar_Richness = lag.listw(listw, residuals_sar_Richness)

### Get data frame containing residuales and spatial lag
dat_sar_Richness = data.frame(Resid_sar_Richness = residuals_sar_Richness, Wsids = Wsids_sar_Richness)

### Plot Moran test using ggplot2
moran.plot_sar_Richness <- ggplot(dat_sar_Richness, aes(x = residuals_sar_Richness, y = Wsids)) + geom_point() + geom_smooth(method = "lm") + xlab("Richness residuals") + ylab("Spatial lag of SAR Richness residuals") + geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed') + theme(text = element_text(size=5.5))
moran.plot_sar_Richness

#####################################################################################

### Spatial simultaneous autoregressive error model estimation for ses.PD

spatial.error_ses.PD <- errorsarlm(ses.PD ~ ele + abs, data = clean_data, listw)
summary(spatial.error_ses.PD, Nagelkerke = TRUE)

### Extract ses.PD residuals
residuals_sar_ses.PD <- residuals.sarlm(spatial.error_ses.PD)
residuals_sar_ses.PD

### Moran's I test for spatial autocorrelation in residuals of ses.PD employing spdep library
moran.mc(residuals_sar_ses.PD, listw, nsim = 999)
moran.plot(residuals_sar_ses.PD, listw)

### Spatial lag of a numeric vector
Wsids_sar_ses.PD = lag.listw(listw, residuals_sar_ses.PD)

### Get data frame containing residuales and spatial lag
dat_sar_ses.PD = data.frame(Resid_sar_ses.PD = residuals_sar_ses.PD, Wsids = Wsids_sar_ses.PD)

### Plot Moran test using ggplot2
moran.plot_sar_ses.PD <- ggplot(dat_sar_ses.PD, aes(x = Resid_sar_ses.PD, y = Wsids)) + geom_point() + geom_smooth(method = "lm") + xlab("ses Phylogenetic Diversity residuals") + ylab("Spatial lag of ses PD residuals") + geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed') + theme(text = element_text(size=5.5))
moran.plot_sar_ses.PD

#####################################################################################

### Spatial simultaneous autoregressive error model estimation for NRI
spatial.error_NRI <- errorsarlm(NRI ~ tem + ele + hfp + abs, data = clean_data, listw)
summary(spatial.error_NRI, Nagelkerke = TRUE)

### Extract NRI residuals
residuals_sar_NRI <- residuals.sarlm(spatial.error_NRI)
residuals_sar_NRI

### Moran's I test for spatial autocorrelation in residuals of NRI employing spdep library
moran.mc(residuals_sar_NRI, listw, nsim = 999)
moran.plot(residuals_sar_NRI, listw)

### Spatial lag of a numeric vector
Wsids_sar_NRI = lag.listw(listw, residuals_sar_NRI)

### Get data frame containing residuales and spatial lag
dat_sar_NRI = data.frame(Resid_sar_NRI = residuals_sar_NRI, Wsids = Wsids_sar_NRI)

### Plot Moran test using ggplot2
moran.plot_sar_NRI <- ggplot(dat_sar_NRI, aes(x = Resid_sar_NRI, y = Wsids)) + geom_point() + geom_smooth(method = "lm") + xlab("NRI residuals") + ylab("Spatial lag of NRI residuals") + geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed') + theme(text = element_text(size=5.5))
moran.plot_sar_NRI

#####################################################################################

### Spatial simultaneous autoregressive error model estimation for NTI
spatial.error_NTI <- errorsarlm(NTI ~ veg + ele + hfp + abs, data = clean_data, listw)
summary(spatial.error_NTI, Nagelkerke = TRUE)

### Extract NTI residuals
residuals_sar_NTI <- residuals.sarlm(spatial.error_NTI)
residuals_sar_NTI

### Moran's I test for spatial autocorrelation in residuals of NTI employing spdep library
moran.mc(residuals_sar_NTI, listw, nsim = 999)
moran.plot(residuals_sar_NTI, listw)

### Spatial lag of a numeric vector
Wsids_sar_NTI = lag.listw(listw, residuals_sar_NTI)

### Get data frame containing residuales and spatial lag
dat_sar_NTI = data.frame(Resid_sar_NTI = residuals_sar_NTI, Wsids = Wsids_sar_NTI)

### Plot Moran test using ggplot2
moran.plot_sar_NTI <- ggplot(dat_sar_NTI, aes(x = Resid_sar_NTI, y = Wsids)) + geom_point() + geom_smooth(method = "lm") + xlab("NTI residuals") + ylab("Spatial lag of NTI residuals") + geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed') + theme(text = element_text(size=5.5))
moran.plot_sar_NTI

#####################################################################################

### Spatial simultaneous autoregressive error model estimation for ses.FD

spatial.error_ses.FD <- errorsarlm(ses.FD ~ veg + ele + hfp + abs, data = clean_data, listw)
summary(spatial.error_ses.FD, Nagelkerke = TRUE)

### Extract ses.FD residuals
residuals_sar_ses.FD <- residuals.sarlm(spatial.error_ses.FD)
residuals_sar_ses.FD

### Moran's I test for spatial autocorrelation in residuals of ses.FD employing spdep library
moran.mc(residuals_sar_ses.FD, listw, nsim = 999)
moran.plot(residuals_sar_ses.FD, listw)

### Spatial lag of a numeric vector
Wsids_sar_ses.FD = lag.listw(listw, residuals_sar_ses.FD)

### Get data frame containing residuales and spatial lag
dat_sar_ses.FD = data.frame(Resid_sar_ses.FD = residuals_sar_ses.FD, Wsids = Wsids_sar_ses.FD)

### Plot Moran test using ggplot2
moran.plot_sar_ses.FD <- ggplot(dat_sar_ses.FD, aes(x = Resid_sar_ses.FD, y = Wsids)) + geom_point() + geom_smooth(method = "lm") + xlab("ses Functional Diversity residuals") + ylab("Spatial lag of ses FD residuals") + geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed') + theme(text = element_text(size=5.5))
moran.plot_sar_ses.FD
#####################################################################################

### Spatial simultaneous autoregressive error model estimation for FNRI
spatial.error_FNRI <- errorsarlm(FNRI ~ ele + abs, data = clean_data, listw)
summary(spatial.error_FNRI, Nagelkerke = TRUE)

### Extract FNRI residuals
residuals_sar_FNRI <- residuals.sarlm(spatial.error_FNRI)
residuals_sar_FNRI

### Moran's I test for spatial autocorrelation in residuals of FNRI employing spdep library
moran.mc(residuals_sar_FNRI, listw, nsim = 999)
moran.plot(residuals_sar_FNRI, listw)

### Spatial lag of a numeric vector
Wsids_sar_FNRI = lag.listw(listw, residuals_sar_FNRI)

### Get data frame containing residuales and spatial lag
dat_sar_FNRI = data.frame(Resid_sar_FNRI = residuals_sar_FNRI, Wsids = Wsids_sar_FNRI)

### Plot Moran test using ggplot2
moran.plot_sar_FNRI <- ggplot(dat_sar_FNRI, aes(x = Resid_sar_FNRI, y = Wsids)) + geom_point() + geom_smooth(method = "lm") + xlab("FNRI residuals") + ylab("Spatial lag of FNRI residuals") + geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed') + theme(text = element_text(size=5.5))
moran.plot_sar_FNRI

#####################################################################################

### Spatial simultaneous autoregressive error model estimation for FNTI
spatial.error_FNTI <- errorsarlm(FNTI ~ tem + ele + hfp + abs, data = clean_data, listw)
summary(spatial.error_FNTI, Nagelkerke = TRUE)

### Extract FNTI residuals
residuals_sar_FNTI <- residuals.sarlm(spatial.error_FNTI)
residuals_sar_FNTI

### Moran's I test for spatial autocorrelation in residuals of FNTI employing spdep library
moran.mc(residuals_sar_FNTI, listw, nsim = 999)
moran.plot(residuals_sar_FNTI, listw)

### Spatial lag of a numeric vector
Wsids_sar_FNTI = lag.listw(listw, residuals_sar_FNTI)

### Get data frame containing residuales and spatial lag
dat_sar_FNTI = data.frame(Resid_sar_FNTI = residuals_sar_FNTI, Wsids = Wsids_sar_FNTI)

### Plot Moran test using ggplot2
moran.plot_sar_FNTI <- ggplot(dat_sar_FNTI, aes(x = Resid_sar_FNTI, y = Wsids)) + geom_point() + geom_smooth(method = "lm") + xlab("FNRI residuals") + ylab("Spatial lag of FNTI residuals") + geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed') + theme(text = element_text(size=5.5))
moran.plot_sar_FNTI

#####################################################################################

### Export Multi Moran plot for SAR models

### Save Multi Moran plot in png format
tiff(filename="..\\data\\output_data\\SAR_Moran_plot.tif", compression = "lzw", units="mm", width=174, height=110, pointsize=9, res=1000)

### Multiplot
multiplot(moran.plot_sar_ses.PD, moran.plot_sar_ses.FD, moran.plot_sar_NRI, 
          moran.plot_sar_FNRI, moran.plot_sar_NTI, moran.plot_sar_FNTI, 
          moran.plot_sar_Richness, cols = 4)

### Provide control over multiple graphics devices
dev.off()

#####################################################################################

### Bubble plots ###

### Read shapefile containing spatial variables
Juarez_map <- readShapeSpatial("..\\data\\raw_data\\Juarez_map.shp")

# Next the shapefile has to be converted to a dataframe for use in ggplot2
Juarez_map <- fortify(Juarez_map)

### Richness
OLS_richness_resid <- data.frame(x_coor, y_coor, resids = residuals_Richness)

bub_OLS_richness_resid <- ggplot(OLS_richness_resid, aes(x = x_coor, y = y_coor, size = resids, fill = resids)) + geom_point(shape = 21) + ggtitle("OLS Richness residuals") + 
  labs(x = "X-coordinates", y = "Y-coordinates", size = "Residuals", fill = "Residuals") + scale_fill_continuous(low = "yellow", high = "red") + guides(size = FALSE) + theme(text = element_text(size=5))
ggplot(OLS_richness_resid, aes(x = x_coor, y = y_coor, size = resids, fill = resids)) + geom_polygon(data=Juarez_map, aes(long, lat, group = group, fill = hole), colour = alpha("darkred", 1/2), size = 0.7) + scale_fill_manual(values = c("skyblue", "white")) + theme(legend.position="none")

### ses.PD
OLS_ses.PD_resid <- data.frame(x_coor, y_coor, resids = residuals_ses.PD)

bub_OLS_ses.PD_resid <- ggplot(OLS_ses.PD_resid, aes(x = x_coor, y = y_coor, size = resids, fill = resids)) + geom_point(shape = 21) + ggtitle("OLS ses PD residuals") + 
  labs(x = "X-coordinates", y = "Y-coordinates", size = "Residuals", fill = "Residuals") + scale_fill_continuous(low = "yellow", high = "red") + guides(size = FALSE) + theme(text = element_text(size=5))

### NRI
OLS_NRI_resid <- data.frame(x_coor, y_coor, resids = residuals_NRI)

bub_OLS_NRI_resid <- ggplot(OLS_NRI_resid, aes(x = x_coor, y = y_coor, size = resids, fill = resids)) + geom_point(shape = 21) + ggtitle("OLS NRI residuals") + 
  labs(x = "X-coordinates", y = "Y-coordinates", size = "Residuals", fill = "Residuals") + scale_fill_continuous(low = "yellow", high = "red") + guides(size = FALSE) + theme(text = element_text(size=5))

### NTI
OLS_NTI_resid <- data.frame(x_coor, y_coor, resids = residuals_NTI)

bub_OLS_NTI_resid <- ggplot(OLS_NTI_resid, aes(x = x_coor, y = y_coor, size = resids, fill = resids)) + geom_point(shape = 21) + ggtitle("OLS NTI residuals") + 
  labs(x = "X-coordinates", y = "Y-coordinates", size = "Residuals", fill = "Residuals") + scale_fill_continuous(low = "yellow", high = "red") + guides(size = FALSE) + theme(text = element_text(size=5))

### ses.FD
OLS_ses.FD_resid <- data.frame(x_coor, y_coor, resids = residuals_ses.FD)

bub_OLS_ses.FD_resid <- ggplot(OLS_ses.FD_resid, aes(x = x_coor, y = y_coor, size = resids, fill = resids)) + geom_point(shape = 21) + ggtitle("OLS ses FD residuals") + 
  labs(x = "X-coordinates", y = "Y-coordinates", size = "Residuals", fill = "Residuals") + scale_fill_continuous(low = "yellow", high = "red") + guides(size = FALSE) + theme(text = element_text(size=5))

### FNRI
OLS_FNRI_resid <- data.frame(x_coor, y_coor, resids = residuals_FNRI)

bub_OLS_FNRI_resid <- ggplot(OLS_FNRI_resid, aes(x = x_coor, y = y_coor, size = resids, fill = resids)) + geom_point(shape = 21) + ggtitle("OLS FNRI residuals") + 
  labs(x = "X-coordinates", y = "Y-coordinates", size = "Residuals", fill = "Residuals") + scale_fill_continuous(low = "yellow", high = "red") + guides(size = FALSE) + theme(text = element_text(size=5))

### FNTI
OLS_FNTI_resid <- data.frame(x_coor, y_coor, resids = residuals_FNTI)

bub_OLS_FNTI_resid <- ggplot(OLS_FNTI_resid, aes(x = x_coor, y = y_coor, size = resids, fill = resids)) + geom_point(shape = 21) + ggtitle("OLS FNTI residuals") + 
  labs(x = "X-coordinates", y = "Y-coordinates", size = "Residuals", fill = "Residuals") + scale_fill_continuous(low = "yellow", high = "red") + guides(size = FALSE) + theme(text = element_text(size=5))

### Export Multi Bubble plots for OLS models

### Save Multi Moran plot in tif format
tiff(filename="..\\data\\output_data\\bub_OLS_multiplot.tif", compression = "lzw", units="mm", width=174, height=110, pointsize=9, res=1000)
multiplot(bub_OLS_ses.PD_resid , bub_OLS_ses.FD_resid , bub_OLS_NRI_resid , bub_OLS_FNRI_resid , bub_OLS_NTI_resid, bub_OLS_FNTI_resid, bub_OLS_richness_resid, cols = 4)
dev.off()

tiff(filename="..\\data\\output_data\\bub_OLS_multiplot.tif", compression = "lzw", units="mm", width=174, height=110, pointsize=5, res=1000)
multiplot(bub_OLS_ses.PD_resid , bub_OLS_ses.FD_resid , bub_OLS_NRI_resid , bub_OLS_FNRI_resid, bub_OLS_NTI_resid, bub_OLS_FNTI_resid, cols = 3)
dev.off()