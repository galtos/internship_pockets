####classification pockets with pharmacophoric parameters####
#Libraries
library(caret)
library(corrplot)
library(car)
library(FactoMineR)
#load data
dt_pharmacophores = read.table("../data/FPCount_save_head_50000.txt", sep = ";", row.names = 1)
dt_72descriptors = read.table("../data/data_72desc.csv", header = T, sep = "\t", row.names = 1, fill=TRUE)
dt_12descriptors = read.table("../data/data_desc.csv", header = T, sep = ",", row.names = 1)

##change row names to PDB_LIG_CHAIN_NBR
#pharmacophores
row_names_pharmacophores = paste(gsub("PDB=| ","",dt_pharmacophores[,1]), gsub(".*/prox5_5/|_prox5_5_res_res-renum.pdb", "", rownames(dt_pharmacophores)), sep = "_")
rownames(dt_pharmacophores) = toupper(row_names_pharmacophores)
#delete coloumn of names
dt_pharmacophores[,1] = NULL

#descripors 72
rownames(dt_72descriptors) 
#descripors 12
rownames(dt_12descriptors)
##check names intersection
#names in both pharmacophores and descriptors
names_all = intersect(rownames(dt_pharmacophores), rownames(dt_12descriptors))




###Random value selection - QUERY DATA SET####
dt = dt_72descriptors
#dt_desc = dt_12descriptors
#dt_desc_ph = merge(dt_desc, dt_pharmacophores, by="row.names")
#row.names(dt_desc_ph) = dt_desc_ph[,1]
#dt_desc_ph[,1] = NULL
#row.names(dt_desc_ph)

#dt = dt_12descriptors
#dt = dt_desc_ph
prct_data = 1/100
index = sample(1:nrow(dt),size = nrow(dt)*prct_data)
dt = dt[index,]
###Valid pockets selection
#delete NA
nrow(dt)
ncol(dt)
dt = na.omit(dt)
nrow(dt)
ncol(dt)
#delete pockets > 60 res and < 14
minimum_size_pocket = 60
maximum_size_pocket = 14
inf_60 = which(dt[,"C_RESIDUES"] <= minimum_size_pocket)
sup_14 = which(dt[,"C_RESIDUES"] >= maximum_size_pocket)
dt = dt[intersect(inf_60,sup_14),]
##Find correlation 12 descriptors
findCorrelation(cor(dt), cutoff = 0.90)
#comparaison données
boxplot(scale(dt))
corrplot(cor(dt))
mat_cor = cor(dt)
which(mat_cor[,] > 0.9)
dt[1,]

####ACP####
dt.acp = PCA(dt, scale.unit = T) #normalisé automatiquement
dt.acp$eig
dt.acp$ind$contrib
dt.acp$var$contrib
plot(dt.acp)

####K-MEANS####
##K-MEANS 10%
dt.kmean = kmeans(scale(dt),5)
library(cluster)
library(fpc)

plot(dt.acp,choix = "ind",col.ind = dt.kmean$cluster, title = "projection kmeans sur ACP")
####K-MEDOIDS####
library(cluster)
library(factoextra)
dt.pam = pam(dt, k = 50)
plot(dt.acp,choix = "ind",col.ind = dt.pam$cluster, title = "projection kmeans sur ACP")
###H clust ###








