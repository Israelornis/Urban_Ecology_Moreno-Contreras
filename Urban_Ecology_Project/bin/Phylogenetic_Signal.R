#####################################################
# Script to quantify phylogenetic signal            #
#####################################################

### Script elaborated by Israel Moreno-Contreras (Ornithology, UNAM)

############### Calculate phylogenetic signal in continuous traits ####################

### Install libraries
library(ape)
library(picante)
library(phytools) ### K's Blomberg, phylogenetic signal

### Read traits file 
traits_zstandardized <- read.csv("..\\data\\input_data\\traits_zstandardized.csv", header = TRUE, row.names = 1)

### Read csv archive of the Juarez bird communities 
comm <- read.csv("..\\data\\raw_data\\Aves_Juarez2.csv", header = TRUE, row.names = 1)

### Read phylogenetic tree 
phylogenetic_tree <- read.tree("..\\data\\input_data\\prunned_tree.tre")

### First, you need to define which trait you want to test and give names to each value according to species
body_mass <- traits_zstandardized[,1]
names(body_mass)<-rownames(traits_zstandardized)

Duration <- traits_zstandardized[,25]
names(Duration)<-rownames(traits_zstandardized)

Mean_frequency <- traits_zstandardized[,26]
names(Mean_frequency)<-rownames(traits_zstandardized)

Entropy <- traits_zstandardized[,27]
names(Entropy)<-rownames(traits_zstandardized)

### Calculate phylogenetic signal, K's Blomberg for single continuous traits
body_mass_sig <- phylosig(phylogenetic_tree, body_mass, method="K", test=TRUE, nsim=999)

duration <- phylosig(phylogenetic_tree, Duration, method="K", test=TRUE, nsim=999)

mean_frequency_sig <- phylosig(phylogenetic_tree, Mean_frequency, method="K", test=TRUE, nsim=999)

entropy_sig <- phylosig(phylogenetic_tree, Entropy, method="K", test=TRUE, nsim=999)

### Obtain vectors (K and P values) for each single continuous trait
K <- c(body_mass_sig$K, duration$K, mean_frequency_sig$K, entropy_sig$K) 

P <- c(body_mass_sig$P, duration$P, mean_frequency_sig$P, entropy_sig$P) 

### Obtain a vector naming each single continuous trait
Traits <- c("Body mass", "Duration", "Mean frequency",
                "Entropy")

### Obtain a data frame containing traits, K and P values
K_Blomberg.results <- cbind(Traits, K, P)

### Export data frame containing Localities, Richnnes, and NRI
write.csv(K_Blomberg.results, file = "..\\data\\input_data\\K_Blomberg.results.csv", row.names = FALSE)

############### Calculate phylogenetic signal in binary traits ########################

### Call libraries
library(ape)
library(caper)

### Read traits file 
traits_zstandardized <- read.csv("..\\data\\input_data\\traits_zstandardized.csv", header = TRUE)

### Read functional tree 
phylogenetic_tree <- read.tree("..\\data\\input_data\\prunned_tree.tre")

### Make label
phylogenetic_tree <- makeLabel(phylogenetic_tree)

### View node.label
phylogenetic_tree$node.label

### Obtain comparative data
Traits_compa <- comparative.data(phy = phylogenetic_tree, data = traits_zstandardized, names.col = Species)

### Calculate the phylogenetic D statisti
Scavenger_PhyloD <- phylo.d(Traits_compa, binvar = Scavenger)
Vertebrates_PhyloD <- phylo.d(Traits_compa, binvar = Vertebrates)
Insects_PhyloD <- phylo.d(Traits_compa, binvar = Insects)
Worms_PhyloD <- phylo.d(Traits_compa, binvar = Worms)
Crustaceans_PhyloD <- phylo.d(Traits_compa, binvar = Crustaceans)
Molluscs_PhyloD <- phylo.d(Traits_compa, binvar = Molluscs)
Fish_PhyloD <- phylo.d(Traits_compa, binvar = Fish)
Seeds_PhyloD <- phylo.d(Traits_compa, binvar = Seeds)
Fruits_PhyloD <- phylo.d(Traits_compa, binvar = Fruits)
Plants_PhyloD <- phylo.d(Traits_compa, binvar = Plants)
Nectar_PhyloD <- phylo.d(Traits_compa, binvar = Nectar)
Sap_PhyloD <- phylo.d(Traits_compa, binvar = Sap)
Foraging_terr_aqua_PhyloD <- phylo.d(Traits_compa, binvar = Foraging_stra_terr_aqua)
Foraging_1_2_PhyloD <- phylo.d(Traits_compa, binvar = Foraging_stra_betw_1_2)
Foraging_above_2_PhyloD <- phylo.d(Traits_compa, binvar = Foraging_stra_above_2m)
Foraging_air_PhyloD <- phylo.d(Traits_compa, binvar = Foraging_stra_air)
Diurnal_PhyloD <- phylo.d(Traits_compa, binvar = Diurnal)
Crepuscular_PhyloD <- phylo.d(Traits_compa, binvar = Crepuscular)
Nocturnal_PhyloD <- phylo.d(Traits_compa, binvar = Nocturnal)
Irruptive_PhyloD <- phylo.d(Traits_compa, binvar = Irruptive)

### Obtain vectors (D and P values) for each single binary trait
D_Scavenger <- Scavenger_PhyloD$DEstimate
PNR_Scavenger <- Scavenger_PhyloD$Pval1
PBM_Scavenger <- Scavenger_PhyloD$Pval0

D_Vertebrates <- Vertebrates_PhyloD$DEstimate
PNR_Vertebrates <- Vertebrates_PhyloD$Pval1
PBM_Vertebrates <- Vertebrates_PhyloD$Pval0

D_Insects <- Insects_PhyloD$DEstimate
PNR_Insects <- Insects_PhyloD$Pval1
PBM_Insects <- Insects_PhyloD$Pval0

D_Worms <- Worms_PhyloD$DEstimate
PNR_Worms <- Worms_PhyloD$Pval1
PBM_Worms <- Worms_PhyloD$Pval0

D_Crustaceans <- Crustaceans_PhyloD$DEstimate
PNR_Crustaceans <- Crustaceans_PhyloD$Pval1
PBM_Crustaceans <- Crustaceans_PhyloD$Pval0

D_Molluscs <- Molluscs_PhyloD$DEstimate
PNR_Molluscs <- Molluscs_PhyloD$Pval1
PBM_Molluscs <- Molluscs_PhyloD$Pval0

D_Fish <- Fish_PhyloD$DEstimate
PNR_Fish <- Fish_PhyloD$Pval1
PBM_Fish <- Fish_PhyloD$Pval0

D_Seeds <- Seeds_PhyloD$DEstimate 
PNR_Seeds <- Seeds_PhyloD$Pval1
PBM_Seeds <- Seeds_PhyloD$Pval0

D_Fruits <- Fruits_PhyloD$DEstimate
PNR_Fruits <- Fruits_PhyloD$Pval1
PBM_Fruits <- Fruits_PhyloD$Pval0

D_Plants <- Plants_PhyloD$DEstimate 
PNR_Plants <- Plants_PhyloD$Pval1
PBM_Plants <- Plants_PhyloD$Pval0

D_Nectar <- Nectar_PhyloD$DEstimate
PNR_Nectar <- Nectar_PhyloD$Pval1
PBM_Nectar <- Nectar_PhyloD$Pval0

D_Sap <- Sap_PhyloD$DEstimate
PNR_Sap <- Sap_PhyloD$Pval1
PBM_Sap <- Sap_PhyloD$Pval0

D_Foraging_terr_aqua <- Foraging_terr_aqua_PhyloD$DEstimate
PNR_Foraging_terr_aqua <- Foraging_terr_aqua_PhyloD$Pval1
PBM_Foraging_terr_aqua <- Foraging_terr_aqua_PhyloD$Pval0

D_Foraging_1_2 <- Foraging_1_2_PhyloD$DEstimate
PNR_Foraging_1_2 <- Foraging_1_2_PhyloD$Pval1
PBM_Foraging_1_2 <- Foraging_1_2_PhyloD$Pval0

D_Foraging_above_2 <- Foraging_above_2_PhyloD$DEstimate
PNR_Foraging_above_2 <- Foraging_above_2_PhyloD$Pval1
PBM_Foraging_above_2 <- Foraging_above_2_PhyloD$Pval0

D_Foraging_air <- Foraging_air_PhyloD$DEstimate
PNR_Foraging_air <- Foraging_air_PhyloD$Pval1
PBM_Foraging_air <- Foraging_air_PhyloD$Pval0

D_Diurnal <- Diurnal_PhyloD$DEstimate
PNR_Diurnal <- Diurnal_PhyloD$Pval1
PBM_Diurnal <- Diurnal_PhyloD$Pval0

D_Crepuscular <- Crepuscular_PhyloD$DEstimate
PNR_Crepuscular <- Crepuscular_PhyloD$Pval1
PBM_Crepuscular <- Crepuscular_PhyloD$Pval0

D_Nocturnal <- Nocturnal_PhyloD$DEstimate
PNR_Nocturnal <- Nocturnal_PhyloD$Pval1
PBM_Nocturnal <- Nocturnal_PhyloD$Pval0

D_Irruptive <- Irruptive_PhyloD$DEstimate
PNR_Irruptive <- Irruptive_PhyloD$Pval1
PBM_Irruptive <- Irruptive_PhyloD$Pval0

### Obtain a vector naming each single binary trait
Traits <- c("Scavenger", "Vertebrates", "Insects", "Worms", "Crustaceans", "Molluscs",
            "Fish", "Seeds", "Fruits", "Plants", "Nectar", "Sap", 
            "Foraging stratum terrestrial/aquatic", "Foraging stratum between 1 and 2 meters",
            "Foraging stratum above 2 meters", "Foraging stratum in the air", "Diurnal", 
            "Crepuscular", "Nocturnal", "Irruptive")

### Obtain a vector for each D value of a single binary trait
D <- c(D_Scavenger, D_Vertebrates, D_Insects, D_Worms, D_Crustaceans, D_Molluscs,
       D_Fish, D_Seeds, D_Fruits, D_Plants, D_Nectar, D_Sap, D_Foraging_terr_aqua,
       D_Foraging_1_2, D_Foraging_above_2, D_Foraging_air, D_Diurnal, D_Crepuscular,
       D_Nocturnal, D_Irruptive)

### Obtain a vector for each D value modified (Comparable to K's Blomberg) of a single binary trait
D_modified <- -1 * c(D_Scavenger, D_Vertebrates, D_Insects, D_Worms, D_Crustaceans, D_Molluscs,
       D_Fish, D_Seeds, D_Fruits, D_Plants, D_Nectar, D_Sap, D_Foraging_terr_aqua,
       D_Foraging_1_2, D_Foraging_above_2, D_Foraging_air, D_Diurnal, D_Crepuscular,
       D_Nocturnal, D_Irruptive) + 1

### Obtain a vector for each P value (no random) of a single binary trait
P_no_random <- c(PNR_Scavenger, PNR_Vertebrates, PNR_Insects, PNR_Worms, PNR_Crustaceans,
                 PNR_Molluscs, PNR_Fish, PNR_Seeds, PNR_Fruits, PNR_Plants, PNR_Nectar,
                 PNR_Sap, PNR_Foraging_terr_aqua, PNR_Foraging_1_2, PNR_Foraging_above_2,
                 PNR_Foraging_air, PNR_Diurnal, PNR_Crepuscular, PNR_Nocturnal, PNR_Irruptive) 
                 
### Obtain a vector for each P value (Brownian motion) of a single binary trait
P_Brownian <- c(PNR_Scavenger, PBM_Vertebrates, PBM_Insects, PBM_Worms, PBM_Crustaceans,
                PBM_Molluscs, PBM_Fish, PBM_Seeds, PBM_Fruits, PBM_Plants, PBM_Nectar,
                PBM_Sap, PBM_Foraging_terr_aqua, PBM_Foraging_1_2, PBM_Foraging_above_2,
                PBM_Foraging_air, PBM_Diurnal, PBM_Crepuscular, PBM_Nocturnal, PBM_Irruptive)

### Obtain a data frame containing traits, D, P1, P0
D_Estimate.results <- cbind(Traits, D, D_modified, P_no_random, P_Brownian)

### Export data frame containing Localities, Richnnes, and NRI
write.csv(D_Estimate.results, file = "..\\data\\input_data\\D_Estimate.results.csv", row.names = FALSE)

############### Calculate phylogenetic signal in categorical traits ##################
library(ape)
library(geiger)
library(phangorn)
library(phylobase)

### Read traits file 
traits_zstandardized <- read.csv("..\\data\\input_data\\traits_zstandardized.csv", header = TRUE, row.names = 1)

### Read phylogenetic tree 
phylogenetic_tree <- read.tree("..\\data\\input_data\\prunned_tree.tre")

### Function to obtain phylogenetic signal for categorical traits 'phylo.signal.disc' developed by Enrico Rezende
'phylo.signal.disc' <-
  function(trait,phy,rep = 999,cost=NULL)
  {
    lev <- attributes(factor(trait))$levels
    if (length(lev) == length(trait))
      stop("Are you sure this variable is categorical?")
    if(is.null(cost)){
      cost1 <- 1-diag(length(lev))
    }
    else {
      if (length(lev) != dim(cost)[1])
        stop("Dimensions of the character state transition matrix do not agree with the number of levels")
      cost1<- t(cost)
    }
    dimnames(cost1) <- list(lev,lev)
    trait <- as.numeric(trait)
    attributes(trait)$names <- phy$tip
    NULL.MODEL <- matrix(NA,rep,1)
    obs <- t(data.frame(trait))
    obs <- phyDat(t(obs),type="USER",levels=attributes(factor(obs))$levels)
    OBS <- parsimony(phy,obs,method="sankoff",cost=cost1)
    for (i in 1:rep){
      null <- sample(as.numeric(trait))
      attributes(null)$names <- attributes(trait)$names
      null <- t(data.frame(null))
      null <- phyDat(t(null),type="USER",levels=attributes(factor(null))$levels)
      NULL.MODEL[i,]<-parsimony(phy,null,method="sankoff",cost=cost1)
      P.value <- sum(OBS >= NULL.MODEL)/(rep + 1)
    }
    par(mfrow=c(1,2))
    hist(NULL.MODEL,xlab="Transitions.in.Randomizations",xlim=c(min(c(min(NULL.MODEL,OBS-1))),max(NULL.MODEL)+1))
    arrows(OBS,rep/10,OBS,0,angle=20,col="red",lwd=4)
    phy$tip.label <- rep(".",length(trait))
    plot(phy,tip.col=trait+10,cex=250/length(trait),font=1)
    title("Character states")
    par(mfrow=c(1,1))
    
    OUTPUT1 <- t(data.frame(Number.of.Levels = length(attributes(factor(trait))$levels), Evolutionary.Transitions.Observed=OBS,Evolutionary.Transitions.Randomization.Median=median(NULL.MODEL),Evolutionary.Transitions.Randomization.Min=min(NULL.MODEL),Evolutionary.Transitions.Randomization.Max=max(NULL.MODEL),P.value))
    
    if(is.null(cost)){
      list(.Randomization.Results=OUTPUT1,.Levels= lev,.Costs.of.character.state.transition.UNORDERED.PARSIMONY = t(cost1))
    }
    else {
      list(.Randomization.Results=OUTPUT1,.Levels= lev,.Costs.of.character.state.transition.FROM.ROW.TO.COL = t(cost1))        }
  }

### Measure phylogenetic signal for categorical traits
For_tec_phylo_sig <- phylo.signal.disc(traits_zstandardized$Foraging_technique, phylogenetic_tree)

Soc_for_phylo_sig <- phylo.signal.disc(traits_zstandardized$Social_Foraging, phylogenetic_tree)

Sea_phylo_sig <- phylo.signal.disc(traits_zstandardized$Seasonality, phylogenetic_tree)

### Results of the phylogenetic signal for categorical traits
## If evolutionary transitions observed are lower than randomization median, then there are phylogenetic signal
Foraging_technique <- as.vector(For_tec_phylo_sig$.Randomization.Results)

Social_foraging <- as.vector(Soc_for_phylo_sig$.Randomization.Results)

Seasonality <- as.vector(Sea_phylo_sig$.Randomization.Results)

### Obtain a vector for results for each categorical trait
Attributes <- c("Number of Levels", "Evolutionary Transitions Observed",
                "Evolutionary Transitions Randomization Median", 
                "Evolutionary Transitions Randomization Min",
                "Evolutionary Transitions Randomization Max", "P value")

### Obtain a data frame containing traits, D, P1, P0
Phylo_Categorical.results <- cbind(Attributes, Foraging_technique, Social_foraging, 
                               Seasonality)

### Export data frame containing Localities, Richnnes, and NRI
write.csv(Phylo_Categorical.results, file = "..\\data\\input_data\\Phylo_Categorical.results.csv", row.names = FALSE)

### Export high-quality figures corresponding to phylogenetic signal analyses for categorical traits
# Foraging technique
tiff(filename="..\\data\\output_data\\Phylo_sig_For_tec.tif", compression = "lzw", units="mm", width=174, height=110, pointsize=9, res=1000)
phylo.signal.disc(traits_zstandardized$Foraging_technique, phylogenetic_tree)
dev.off()

# Social foraging
tiff(filename="..\\data\\output_data\\Phylo_sig_Soc_for.tif", compression = "lzw", units="mm", width=174, height=110, pointsize=9, res=1000)
phylo.signal.disc(traits_zstandardized$Social_Foraging, phylogenetic_tree)
dev.off()

# Seasonality
tiff(filename="..\\data\\output_data\\Phylo_sig_Sea.tif", compression = "lzw", units="mm", width=174, height=110, pointsize=9, res=1000)
phylo.signal.disc(traits_zstandardized$Seasonality, phylogenetic_tree)
dev.off()