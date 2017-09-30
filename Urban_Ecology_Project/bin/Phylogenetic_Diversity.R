#####################################################
# Script to quantify phylogenetic diversity metrics #
#####################################################

### Script elaborated by Israel Moreno-Contreras (Ornithology, UNAM)

### Install libraries
library(ape)
library(picante)

### Read tree 
tree <- read.tree("..\\data\\raw_data\\RAxML_bipartitions.TOTAL_bipar")

### Read csv archive of the Juarez bird communities 
comm <- read.csv("..\\data\\raw_data\\Aves_Juarez2.csv", header = TRUE, row.names = 1)

### Test if tree is ultrametric 
is.ultrametric(tree)

### Convert tree into ultrametric tree 
dated_tree <- chronos(tree)

### Test if tree is ultrametric ("dated_tree")
is.ultrametric(dated_tree)

### Prune the tree before running pd
prunedTree <- prune.sample(comm, dated_tree)

### Export Write tree file in parenthetic Format
write.tree(prunedTree, file = "..\\data\\input_data\\prunned_tree.tre")

### Read functional tree 
prunedTree <- read.tree("..\\data\\input_data\\prunned_tree.tre")

### computes the pairwise distances between the pairs of 
### tips from a phylogenetic tree using its branch lengths
phydist <- cophenetic.phylo(prunedTree)

### Calculate Faith's Phylogenetic Diversity 
pd.result <- pd(comm, prunedTree, include.root = TRUE)

### Calculate Standardized effect size of phylogenetic diversity (Faith's PD) in communities
ses.pd.result <- ses.pd(comm, prunedTree, null.model = "taxa.labels",
                        runs = 999, iterations = 1000)

### Visualize Standardized effect size of phylogenetic diversity (Faith's PD) in communities
View(ses.pd.result)

### Export ses.pd results
write.csv(ses.pd.result, file = "..\\data\\input_data\\ses_pd_result.csv")

### Calculate Standardized effect size of MPD
ses.mpd.result <- ses.mpd(comm, phydist, null.model = "taxa.labels",
        runs = 999, iterations = 1000)

### Visualize Standardized effect size of phylogenetic diversity (Faith's PD) in communities
View(ses.mpd.result)

### Export ses.mpd results
write.csv(ses.mpd.result, file = "..\\data\\input_data\\ses_mpd_result.csv")

### Calculate Standardized effect size of MNTD
ses.mntd.result <- ses.mntd(comm, phydist, null.model = "taxa.labels",
                          runs = 999, iterations = 1000)

### Visualize Standardized effect size of MNTD
View(ses.mntd.result)

### Export ses.mntd results
write.csv(ses.mntd.result, file = "..\\data\\input_data\\ses_mntd_result.csv")

### Calculate NRI values
NRI <- ses.mpd.result$mpd.obs.z * -1

### Calculate NTI values
NTI <- ses.mntd.result$mntd.obs.z * -1

### Obtain a vector for Localities, Richnnes, PD, SES.PD
Localities <- c("El Chamizal", "Club Campestre", "Parque Central",
                "Puerta Juarez", "Trepachanga", "Zaragoza", "Loma Blanca")

Richness <- ses.mpd.result$ntaxa

PD <- ses.pd.result$pd.obs

ses.PD <- ses.pd.result$pd.obs.z

### Obtain a data frame containing Localities, Richnnes, PD, SES.PD, NRI, and NTI
PD_metrics.results <- cbind(Localities, Richness, PD, ses.PD, NRI, NTI)

### Export data frame containing Localities, Richnnes, and NRI
write.csv(PD_metrics.results, file = "..\\data\\input_data\\PD_metrics.results.csv")
