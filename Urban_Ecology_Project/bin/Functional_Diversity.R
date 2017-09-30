#####################################################
# Script to quantify functional diversity metrics   #
#####################################################

### Script elaborated by Israel Moreno-Contreras (Ornithology, UNAM)

### Install libraries
library(ape)
library(picante)
library(FD) ### To obtain Gower distances

### Read traits file without row names
traits <- read.csv("..\\data\\raw_data\\Traits_Juarez_birds.csv", header = TRUE)

### Convert continuous variables into z-standardized variables
Body_mass_g_scaled <- as.vector(scale(traits$Body_mass_g))
Duration_scaled <- as.vector(scale(traits$Duration))
Mean_frequency_scaled <- as.vector(scale(traits$Mean_frequency))
Entropy_scaled <- as.vector(scale(traits$Entropy))

### Rounding continuous variables previously z-standardized
Body_mass_g_scaled <- round(Body_mass_g_scaled, digits = 3)
Duration_scaled <- round(Duration_scaled, digits = 3)
Mean_frequency_scaled <- round(Mean_frequency_scaled, digits = 3)
Entropy_scaled <- round(Entropy_scaled, digits = 3)

### Obtain vectors related with functional traits
Species <- as.vector(traits$Species)
Scavenger <- as.vector(traits$Scavenger)
Vertebrates <- as.vector(traits$Vertebrates)
Insects <- as.vector(traits$Insects)
Worms <- as.vector(traits$Worms)
Crustaceans <- as.vector(traits$Crustaceans)
Molluscs <- as.vector(traits$Molluscs)
Fish <- as.vector(traits$Fish)
Seeds <- as.vector(traits$Seeds)
Fruits <- as.vector(traits$Fruits)
Plants <- as.vector(traits$Plants)
Nectar <- as.vector(traits$Nectar)
Sap <- as.vector(traits$Sap)
Foraging_technique <- as.vector(traits$Foraging_technique)
Foraging_stra_terr_aqua <- as.vector(traits$Foraging_stratum_terr.aquatic)
Foraging_stra_betw_1_2 <- as.vector(traits$Foraging_stratum_betw_1_and_2)
Foraging_stra_above_2m <- as.vector(traits$Foraging_stratum_above_2m)
Foraging_stra_air <- as.vector(traits$Foraging_stratum_air)
Social_Foraging <- as.vector(traits$Social_Foraging)
Diurnal <- as.vector(traits$Diurnal)
Crepuscular <- as.vector(traits$Crepuscular)
Nocturnal <- as.vector(traits$Nocturnal)
Seasonality <- as.vector(traits$Seasonality)
Irruptive <- as.vector(traits$Irruptive)

### Obtain a table containing functional z-standardized
traits_zstandardized <- cbind(Species, Body_mass_g_scaled, Scavenger, Vertebrates, Insects,
                       Worms, Crustaceans, Molluscs, Fish, Seeds, Fruits, Plants, Nectar, Sap,
                       Foraging_technique, Foraging_stra_terr_aqua, Foraging_stra_betw_1_2,
                       Foraging_stra_above_2m, Foraging_stra_air, Social_Foraging, Diurnal,
                       Crepuscular, Nocturnal, Seasonality, Irruptive, Duration_scaled,
                       Mean_frequency_scaled, Entropy_scaled)

### Convert traits_zstandardized into data frame
traits_zstandardized <- as.data.frame(traits_zstandardized)

### Export functional traits z-standardized without row names
write.csv(traits_zstandardized, file = "..\\data\\input_data\\traits_zstandardized.csv", row.names = FALSE)

### Read traits file without row names
traits_zstandardized <- read.csv("..\\data\\input_data\\traits_zstandardized.csv", header = TRUE, row.names = 1)

### Read csv archive of the Juarez bird communities 
comm <- read.csv("..\\data\\raw_data\\Aves_Juarez2.csv", header = TRUE, row.names = 1)

### Create a Gower Dissimilarity matrix from traits data
gower_distance <- gowdis(traits_zstandardized, ord = "podani")

### Obtain functional dendogram using a UPGMA algorithm
func_den <- hclust(gower_distance, method = "average")

### Plot a cluster (dendogram) using a UPGMA algorithm
plot(func_den)

### Convert dendrogram into "phylo" format
functional_tree <- as.phylo(func_den)

### Export Write tree file in parenthetic Format
write.tree(functional_tree, file = "..\\data\\input_data\\functional_tree_scaled.tre")

### Convert "gower_distance" into a matrix
matrix_functional <- as.matrix(gower_distance)

### Read functional tree 
functional_tree <- read.tree("..\\data\\input_data\\functional_tree_scaled.tre")

### Calculate Standardized effect size of functional diversity (Petchey and Gaston 2006) in communities
ses.fd.result <- ses.pd(comm, functional_tree, null.model = "taxa.labels",
       runs = 999, iterations = 1000)

### Visualize Standardized effect size of functional diversity (Petchey and Gaston 2006) in communities
View(ses.fd.result$pd.obs.z)

### Obtain a vector for FD and SES.FD
FD <- ses.fd.result$ntaxa

ses.FD <- ses.fd.result$pd.obs.z

### Obtain a table containing Localities, Richnnes, FD, SES.FD, FNRI, and FNTI
ses.fd.result <- cbind(FD, ses.FD)

### Convert table into data frame
ses.fd.result <- as.data.frame(ses.fd.result)

### Export ses.fd results
write.csv(ses.fd.result, file = "..\\data\\input_data\\ses_fd_result.csv")

### Calculate Standardized effect size of FMPD 
ses.fmpd.result <- ses.mpd(comm, matrix_functional, null.model = "taxa.labels",
                          runs = 999, iterations = 1000)

### Visualize Standardized effect size of phylogenetic diversity (Faith's PD) in communities
View(ses.fmpd.result)

### Export ses.mpd results
write.csv(ses.fmpd.result, file = "..\\data\\input_data\\ses_fmpd_result.csv")

### Calculate Standardized effect size of FMNTD
ses.fmntd.result <- ses.mntd(comm, matrix_functional, null.model = "taxa.labels",
                            runs = 999, iterations = 1000)

### Visualize Standardized effect size oof FMNTD
View(ses.fmntd.result)

### Export ses.fmntd results
write.csv(ses.fmntd.result, file = "..\\data\\input_data\\ses_fmntd_result.csv")

### Calculate FNRI values
FNRI <- ses.fmpd.result$mpd.obs.z * -1

### Calculate FNTI values
FNTI <- ses.fmntd.result$mntd.obs.z * -1

### Obtain a vector for Localities, Richnnes, FD, SES.FD
Localities <- c("El Chamizal", "Club Campestre", "Parque Central",
                "Puerta Juarez", "Trepachanga", "Zaragoza", "Loma Blanca")

Richness <- ses.fmpd.result$ntaxa

FD <- ses.fd.result$FD

ses.FD <- ses.fd.result$ses.FD

### Obtain a table containing Localities, Richnnes, FD, SES.FD, FNRI, and FNTI
FD_metrics.results <- cbind(Localities, Richness, FD, ses.FD, FNRI, FNTI)

### Export data frame containing Localities, Richnnes, FD, SES.FD, FNRI, and FNTI
write.csv(FD_metrics.results, file = "..\\data\\input_data\\FD_metrics.results.csv")