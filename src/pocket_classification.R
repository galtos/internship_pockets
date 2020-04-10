####classification pockets with pharmacophoric parameters####
#Libraries
library(caret)
library(corrplot)
library(car)
library(FactoMineR)
####load data####
dt_pharmacophores = read.table("../data/FPCount_save_head_50000.txt", sep = ";", row.names = 1)
dt_72descriptors = read.table("../data/data_72desc.csv", header = T, sep = "\t", row.names = 1, fill=TRUE)
dt_12descriptors = read.table("../data/data_desc.csv", header = T, sep = ",", row.names = 1)
dt_drug = read.table("../data/data_drug.txt", header = T, sep = ";", row.names = 1)
##change row names to PDB_LIG_CHAIN_NBR
#deal with big FP file

processFile = function(filepath, dt, fileConn) {
  con = file(filepath, "r")
  dt_pharmacophores = NULL
  rnames_dt = toupper(rownames(dt))
  i=0
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if (length(line) == 0) {
      close(con)
      break
    }
    line_split = unlist(strsplit(line, split=";"))
    name_pharmacophore = paste(gsub("PDB=| ","",line_split[2]), gsub(".*/prox5_5/|_prox5_5_res_res-renum.pdb", "", line_split[1]), sep = "_")
    name_pharmacophore = toupper(name_pharmacophore)
    #print(name_pharmacophore)
    line_split[1] = name_pharmacophore
    line_split = line_split[-2]
    #print(paste(line_split, collapse =";"))
    if(is.element(name_pharmacophore, rnames_dt)) {
      #dt_pharmacophores = rbind(dt_pharmacophores, line_split[-2])
      print("------------YES------------")
      i=i+1
      print(i)
      cat(paste(line_split, collapse =";"), file = "../data/FPCount_save_all_inter_dt12_bis.txt",append=TRUE)
      cat("\n",file = "../data/FPCount_save_all_inter_dt12_bis.txt",append=TRUE)

    }
  }
  close(con)
}
processFile("../data/FPCount_save_all.txt", dt, fileConn)

dt_pharmacophores = read.table("../data/FPCount_save_all_inter_dt12_bis.txt", sep = ";", row.names = 1, nrows = 20)


rownames(dt_pharmacophores) = toupper(dt_pharmacophores[,1])
dt_pharmacophores = as.data.frame(dt_pharmacophores)
dt_pharmacophores[,1] = NULL
colnames(dt_pharmacophores) = 1:ncol(dt_pharmacophores)
write.table(dt_pharmacophores, "../data/FPCount_save_element_dt12.txt")

#pharmacophores
row_names_pharmacophores = paste(gsub("PDB=| ","",dt_pharmacophores[,1]), gsub(".*/prox5_5/|_prox5_5_res_res-renum.pdb", "", rownames(dt_pharmacophores)), sep = "_")
rownames(dt_pharmacophores) = toupper(row_names_pharmacophores)
#delete coloumn of names
dt_pharmacophores[,1] = NULL

#descripors 72
r_names_dt_72 = toupper(gsub("_prox5_5.desR","",rownames(dt_72descriptors)))
rownames(dt_72descriptors) = toupper(gsub("_prox5_5.desR","",rownames(dt_72descriptors)))

dt_72descriptors = intersect(r_names_dt_72, rownames(dt))
rownames(dt_72descriptors) 
nrow(dt_72descriptors)
#drug like
row.names(dt_drug) = toupper(gsub("_prox5_5.predict","",rownames(dt_drug)))

####Data management Desc72####
#n values
nrow(dt_72descriptors)
#descripors 72
summary(dt_72descriptors)
dt_72descriptors["HEIGHT_MIN_CYLINDER"] = NULL
dt_72descriptors["RADIUS_MIN_CYLINDER"] = NULL
dt_72descriptors["VOLUME_HULL"] = NULL
dt_72descriptors["drugg"] = NULL
summary(dt_72descriptors)
dt_72descriptors = na.omit(dt_72descriptors)
names_all_infos = (which(nchar(rownames(dt_72descriptors))>20))
nrow(dt_72descriptors)
dt_72descriptors = dt_72descriptors[names_all_infos,]

length(intersect(row.names(dt),row.names(dt_drug)))
##check names intersection
#names in both pharmacophores and descriptors
names_all = intersect(rownames(dt_pharmacophores), rownames(dt_12descriptors))
####Data management drug data####
nrow(dt_drug)
length(intersect(row.names(dt_drug),row.names(dt_72descriptors)))
dt_drug = na.omit(dt_drug)
length(which(dt_drug[,"d"] < 0.5))

###Random value selection - QUERY DATA SET####
#dt = dt_72descriptors
#dt_desc = dt_12descriptors
#dt_desc_ph = merge(dt_desc, dt_pharmacophores, by="row.names")
#row.names(dt_desc_ph) = dt_desc_ph[,1]
#dt_desc_ph[,1] = NULL
#row.names(dt_desc_ph)

dt = dt_12descriptors
#dt = dt_desc_ph
prct_data = 1/100
index = sample(1:nrow(dt),size = nrow(dt)*prct_data)
dt = dt[index,]
###Valid pockets selection
delete_clean_data = function(dt){
  print("number pockets NA")
  print(nrow(dt) - nrow(na.omit(dt)))  
  dt = na.omit(dt)
  minimum_size_pocket = 60
  maximum_size_pocket = 14

  inf_60 = which(dt[,"C_RESIDUES"] <= minimum_size_pocket)
  print("number of pockets superior 60:")
  print(nrow(dt)-length(inf_60))
  sup_14 = which(dt[,"C_RESIDUES"] >= maximum_size_pocket)
  print("number of pockets inferior 14:")
  print(nrow(dt)-length(sup_14))
  
  print("Number of pockets >=14 <=60")
  print(length(intersect(inf_60,sup_14)))
  dt = dt[intersect(inf_60,sup_14),]
  
  print("delete DOD, NTN, EDO, SF4:")
  sup_DOD = grep("_DOD_", rownames(dt))
  print("Number DOD:")
  print(length(sup_DOD))
  sup_NTN = grep("_NTN_", rownames(dt))
  print("Number NTN:")
  print(length(sup_NTN))
  sup_EDO = grep("_EDO_", rownames(dt))
  print("Number EDO:")
  print(length(sup_EDO))
  sup_SF4 = grep("_SF4_", rownames(dt))
  print("Number SF4:")
  print(length(sup_SF4))
  
  dt = dt[-c(sup_DOD,sup_NTN,sup_EDO,sup_SF4),]
  
  print("number invalid pockets:")
  load("../results/names_pockets_double")
  print(length(names_pockets_double))
  dt = dt[-which(is.element(rownames(dt), names_pockets_double)),]
  print("nrow jeu final:")
  print(nrow(dt))
  return(dt)
} 
nrow(dt)
ncol(dt)
dt = na.omit(dt)
nrow(dt)
ncol(dt)
#delete pockets >= 60 res and <= 14
minimum_size_pocket = 60
maximum_size_pocket = 14
inf_60 = which(dt[,"C_RESIDUES"] > minimum_size_pocket)
sup_14 = which(dt[,"C_RESIDUES"] < maximum_size_pocket)
length(inf_60)
length(sup_14)
dt = dt[intersect(inf_60,sup_14),]
##Find correlation 12 descriptors
column_cor_remove = findCorrelation(cor(dt), cutoff = 0.60, verbose = T, names = T)
length(column_cor_remove)
#remove correlated columns
dt = dt[,-c(column_cor_remove)]

colnames(dt[,24])
colnames(dt[,53])
#comparaison donnÃ©es
boxplot(scale(dt))
corrplot(cor(dt))
mat_cor = cor(dt)
which(mat_cor[,] > 0.9)
dt[1,]
####add drug columns####
dt_drug[,2] = NULL
dt_drug[which(dt_drug[,1] > 0.5),1]= 1
dt_drug[which(dt_drug[,1] < 0.5),1]= 0
dt_drug = na.omit(dt_drug)
dt = merge(dt, dt_drug, by="row.names")
row.names(dt) =  dt$Row.names
dt$Row.names = NULL
####ACP####
#dt.acp = PCA(dt, scale.unit = T,quali.sup=13) #normalisÃ© automatiquement
dt.acp = PCA(dt_72, scale.unit = T)
dt.acp$eig
dt.acp$ind$contrib
dt.acp$var$contrib
plot(dt.acp)
#plot(dt.acp,habillage = 13, col.hab = c("green","blue"), label="none")
library(rgl)
plot3d(dt.acp$scores[,1:3], col=dt.kmean$cluster)
####K-MEANS####
##K-MEANS nK = 10%; nK = 10
dt = dt_12descriptors
dt = delete_clean_data(dt)
#
nbr_k = nrow(dt)*20/100
nbr_k = as.integer(nbr_k)
nbr_k = 10
start_time <- Sys.time()
dt.kmean = kmeans(scale(dt), nbr_k, nstart = 1)
end_time <- Sys.time()

end_time - start_time
#test
dt.kmean = NULL
dt.kmean = c(dt.kmean, kmeans(scale(dt), nbr_k, nstart = 1))
dt.kmean[1]
##K-MEANS nK = 1-->20
start_time <- Sys.time()
R2 = NULL
inertie.expl = NULL
n=50
for (i in 1:n) {
  k = kmeans(dt, i, nstart = 1, algorithm="MacQueen", iter.max = 200) #change nstart to 10
  R2 = c(R2,k$tot.withinss/k$totss)
  #inertie.expl = c(inertie.expl,k$betweenss/k$totss)
  #Iintra = c(Iintra.expl,k$betweenss)
}
end_time <- Sys.time()

end_time - start_time
#tot.withniss : I intra total
#withinss : I intra
#betweenss : I entre deux
#totss : I totale
#R2 = I.intra/I.total
png(filename=paste0(paste0("../results/kmeans_Iintra_tot_tot_n",n),".png"))
plot(y = R2,x = 1:n, type = 'l')
dev.off()
write(R2, paste0("../results/kmeans_R2_n",n))
png(filename=paste0(paste0("../results/kmeans_Iintra_tot_n",n),".png"))
plot(y = inertie.expl,x = 1:n, type = 'l',xlab="Nb. groups",ylab="% inertia explained")
dev.off()
write(inertie.expl, paste0("../results/kmeans_inertie_expl_n",n))

##Analysis
#number of pockets in each clusters 
#max min clusters
hist(dt.kmean$cluster)

max(table(dt.kmean$cluster))
min(table(dt.kmean$cluster))
plot(table(dt.kmean$cluster))
plot(dt.acp,choix = "ind", col.quali = dt.kmean$cluster, title = "projection kmeans sur ACP", label="none")
#autoplot
library(ggfortify)
autoplot(dt.kmean, data = dt)
## Centroids
dt_kmean_centroids
dt.kmean$centers[1,]
####K MEANS cascade ####
n_iterations = 2
nbr_k_1 = 15
nbr_k_2 = 5
dt.kmean_1 = kmeans(scale(dt), nbr_k_1, nstart = 10)

dt.kmean_2 = NULL
for (i in 1:n_iteration) {
  dt.kmean_2 = c(dt.kmean_2, kmeans(scale(dt[which(dt.kmean_1$cluster == 1),]), nbr_k_2, nstart = 10))

}
####Hierarchical Classification Multiple after K-Means####
dt.hclust_1 = hclust(dist(scale(dt[names(which(dt.kmean$cluster == 1)),])), method = "ward.D2")
dt.hclust_2 = hclust(dist(scale(dt[names(which(dt.kmean$cluster == 2)),])), method = "ward.D2")
dt.hclust_3 = hclust(dist(scale(dt[names(which(dt.kmean$cluster == 3)),])), method = "ward.D2")
dt.hclust_4 = hclust(dist(scale(dt[names(which(dt.kmean$cluster == 4)),])), method = "ward.D2")
dt.hclust_5 = hclust(dist(scale(dt[names(which(dt.kmean$cluster == 5)),])), method = "ward.D2")
dt.hclust_6 = hclust(dist(scale(dt[names(which(dt.kmean$cluster == 6)),])), method = "ward.D2")
dt.hclust_7 = hclust(dist(scale(dt[names(which(dt.kmean$cluster == 7)),])), method = "ward.D2")
dt.hclust_8 = hclust(dist(scale(dt[names(which(dt.kmean$cluster == 8)),])), method = "ward.D2")
dt.hclust_9 = hclust(dist(scale(dt[names(which(dt.kmean$cluster == 9)),])), method = "ward.D2")
dt.hclust_10 = hclust(dist(scale(dt[names(which(dt.kmean$cluster == 10)),])), method = "ward.D2")
dt.hclust_11 = hclust(dist(scale(dt[names(which(dt.kmean$cluster == 11)),])), method = "ward.D2")
dt.hclust_12 = hclust(dist(scale(dt[names(which(dt.kmean$cluster == 12)),])), method = "ward.D2")
dt.hclust_13 = hclust(dist(scale(dt[names(which(dt.kmean$cluster == 13)),])), method = "ward.D2")
dt.hclust_14 = hclust(dist(scale(dt[names(which(dt.kmean$cluster == 14)),])), method = "ward.D2")
dt.hclust_15 = hclust(dist(scale(dt[names(which(dt.kmean$cluster == 15)),])), method = "ward.D2")
dt.hclust_16 = hclust(dist(scale(dt[names(which(dt.kmean$cluster == 16)),])), method = "ward.D2")
dt.hclust_17 = hclust(dist(scale(dt[names(which(dt.kmean$cluster == 17)),])), method = "ward.D2")
dt.hclust_18 = hclust(dist(scale(dt[names(which(dt.kmean$cluster == 18)),])), method = "ward.D2")
dt.hclust_19 = hclust(dist(scale(dt[names(which(dt.kmean$cluster == 19)),])), method = "ward.D2")
dt.hclust_20 = hclust(dist(scale(dt[names(which(dt.kmean$cluster == 20)),])), method = "ward.D2")


pockets_cluster = NULL
dt.visualize = NULL
start_time <- Sys.time()
for (i in 1:nbr_k) {
  print(i)
  pockets_cluster = names(which(dt.kmean$cluster == i))
  dt.hclust = hclust(dist(scale(dt[pockets_cluster,])), method = "ward.D2")
  dt.visualize = rbind(as.data.frame(t(apply(dt[pockets_cluster,],2,mean))),dt.visualize)
  png(filename=paste0(paste0("../results/hclust_",i),".png"))
  plot(dt.hclust,main = "", hang = -0.1)
  dev.off()
}
end_time <- Sys.time()

end_time - start_time

plot(dt.hclust,main = "", hang = -0.1)
####Hierarchical Classification data size dependent ####
dt = delete_clean_data(dt_12descriptors)
prct_data = 40/100
index = sample(1:nrow(dt),size = nrow(dt)*prct_data)
dt = dt[index,]
nrow(dt)
start_time <- Sys.time()
dt.hclust = hclust(dist(scale(dt)), method = "ward.D2")
end_time <- Sys.time()

end_time - start_time

plot(dt.hclust,main = "", hang = -0.1)
####Hierarchical Classification on centroids####
dt_centers.hclust = hclust(dist(scale(dt.kmean$centers)), method = "ward.D2")
png(filename="../results/res_12desc/hclust_centroids_k20.png")
png(filename="../results/res_12desc/hclust_centroids_k10.png")
plot(dt_centers.hclust,main = "", hang = -0.1)
dev.off()
#choix nbr classes
inertie <- sort(dt_centers.hclust$height, decreasing = TRUE)
png(filename="../results/res_12desc/hclust_centroids_k20_inertie.png")
png(filename="../results/res_12desc/hclust_centroids_k10_inertie.png")
plot(inertie[1:30], type = "s", xlab = "Nombre de classes", ylab = "Inertie")
points(c(5, 11, 15), inertie[c(5, 11, 15)], col = c("green3", "red3", 
                                                "blue3"), cex = 2, lwd = 3)
dev.off()
N_classes_selected = 8
##K-MEANS nK = selected
dt_.kmean = kmeans(scale(dt), N_classes_selected)
####Cluster visualization####
library(fmsb)
radarchart(langues.means)
#dt.visualize = as.data.frame(t(apply(scale(dt[1:5,]),2,mean)))
#dt.visualize = NULL
#dt.visualize = rbind(as.data.frame(t(apply(scale(dt[5:10,]),2,mean))),dt.visualize)

#add max and min
dt.visualize = as.data.frame(scale(dt.visualize))
dt.visualize =  rbind(apply(dt.visualize,2,max), apply(dt.visualize,2,min) , dt.visualize)

# Color vector
colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )
radarchart(dt.visualize, axistype=1 , 
           #custom polygon
           pcol=colors_border , plwd=4 , plty=1,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(-1,1,0.5), cglwd=0.8,
           #custom labels
           vlcex=0.8 
)
for (i in 1:nbr_k) {
  png(filename=paste0(paste0("../results/res_12desc/radarchart_cluster_",i),".png"))
  radarchart(dt.visualize[c(1,2,i+2),], axistype=1 , 
             #custom polygon
             plwd=4 , plty=1,
             #custom the grid
             cglcol="grey", caxislabels=seq(-1,1,0.5), cglty=1, axislabcol="grey", cglwd=0.8,
             #custom labels
             vlcex=0.8 
  )
  dev.off()
}
radarchart(dt.visualize[c(1,2,2+2),], axistype=1 , 
           #custom polygon
           plwd=4 , plty=1, pcol = 1,
           #custom the grid
           cglcol="grey", caxislabels=seq(-1,1,0.5), cglty=1, axislabcol="grey", cglwd=0.8,
           #custom labels
           vlcex=0.8 
)
legend(x=1, y=1.2, legend = c("c1","c2","c3"), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)
##convert to newick format
#install.packages("https://cran.r-project.org/src/contrib/Archive/amap/amap_0.8-14.tar.gz", repos=NULL)
#library(amap)

install.packages("https://www.bioconductor.org/packages/release/bioc/src/contrib/ctc_1.60.0.tar.gz", repos=NULL)
#library(ctc)
#write(hc2Newick(dt.hclust_1),file='../results/res_12desc/hclust_1.newick')
#tree
tree.trait<- as.phylo(dt.hclust_1)
write.tree(tree.trait, '../results/res_12desc/hclust_1.newick')


install.packages("phylogram")
library("phylogram")
dendrogram.hclust_1 = read.dendrogram('../results/res_12desc/hclust_1.newick')
library("ape")

pdf("../results/res_12desc/hclust_1.pdf", width=40, height=15)
# Do some plotting
plot(as.phylo(dendrogram.hclust_1))
# Close the PDF file's associated graphics device (necessary to finalize the output)
dev.off()


##TEST for number max visualization
dt.hclust = hclust(dist(scale(dt[1:200,])), method = "ward.D2")
plot(dt.hclust)

####Similarity of pharmacophores####
#based on fuzcav similarity score

mat_similarity_pharmacophores = function(dt){
  mat_similarity = matrix(data = NA, nrow = nrow(dt), ncol = nrow(dt))
  row.names(mat_similarity) = row.names(dt)
  colnames(mat_similarity) = row.names(dt)
  for (i in row.names(dt)) {
    for (j in row.names(dt)) {
      t = table(unlist(dt_pharmacophores[i,]), unlist(dt_pharmacophores[j,]))
      common_non_null = sum(t[2:nrow(t),2:ncol(t)])
      non_null_counts_Fa = sum(t[2:nrow(t),])
      non_null_counts_Fb = sum(t[,2:ncol(t)])
      mat_similarity[i,j] = common_non_null/min(non_null_counts_Fa,non_null_counts_Fb)
    }
  }
  return(mat_similarity)
}
mat_similarity_ph = mat_similarity_pharmacophores(dt_pharmacophores[1:10,])

library(pheatmap)
pheatmap(mat_similarity_ph, cex = -0.9)
##TEST hclust pharmacophores##
dt.kmean = kmeans(scale(dt_pharmacophores), 10, nstart = 10)

###### TRUE WORKFLOW FOR CLASSIFICATION #####
library(data.tree)
library(treemap)
library(gtools)
library(ade4)
#lets suppose we already cleaned the data: dt
nstart = 100
prct_seed = 10/100#10/100
min_size_cluster = 100

pockets_classification_tree = function(dt,
                                       nstart = 100,
                                       path_tree = "alltree",
                                       dt_pharmacophores){ # be careful to scale the data before
  #prct_seed=10/100
  #first k mean
  #nbr_k = nrow(dt)*prct_seed # select number of seed : 10% of the size of the data
  #nbr_k = as.integer(nbr_k)

  #dt.kmean = kmeans(scale(dt), nbr_k, nstart = nstart)
  #hclust on centroids

  #dt_centers.hclust = hclust(dist(scale(dt.kmean$centers)), method = "ward.D2")
  ##
  #TODO:select K optimal
  ##
  nbr_k_optimal = 800
  #seceond kmean
  print("here1")
  #dt.kmean = kmeans(dt, centers = nbr_k_optimal, nstart = nstart, iter.max = 200)
  #dt.kmean = kmeans(dt, nbr_k_optimal, nstart = nstart, algorithm="MacQueen", iter.max = 200)

  print("here2")
  pockets_cluster = list()
  cluster_dend = list()
  #pharmacophores_consensus_mean_cluster = list()
  #pharmacophores_consensus_50_cluster = list()
  for (i in 1:nbr_k_optimal){
    pockets_cluster[[i]] = names(which(dt.kmean$cluster == i))
    cluster_dend[[i]] = as.dendrogram(hclust(dist(dt[which(dt.kmean$cluster == i),]),
                                      method = "ward.D2"))
    #pharmacophores_consensus_mean_cluster[[i]] = apply(dt_pharmacophores[intersect(rownames(dt_pharmacophores),
    #                                  names(which(dt.kmean$cluster == i))),], 2,
    #                                  function(x) v <- as.integer(mean(x)))
    #pharmacophores_consensus_50_cluster[[i]] = apply(dt_pharmacophores[intersect(rownames(dt_pharmacophores),
    #                                                                names(which(dt.kmean$cluster == i))),], 2, function(x) {
    #                                                                if(length(which(x>=1)) >= length(x)/2) {
    #                                                                  return(as.integer(mean(x)))
    #                                                                } else {
    #                                                                  return(0)
    #                                                                }
    #                                                              })
  }
  cluster_dt = data.frame(centers = dt.kmean$centers,
                          withinss = dt.kmean$withinss,
                          size = dt.kmean$size,
                          betweenss = rep(dt.kmean$betweenss, nbr_k_optimal),
                          totss = rep(dt.kmean$totss, nbr_k_optimal),
                          pockets_names = I(pockets_cluster),
                          cluster_dend = I(cluster_dend)
                          #pharmacophores_consensus_mean = I(pharmacophores_consensus_mean_cluster),
                          #pharmacophores_consensus_50 = I(pharmacophores_consensus_50_cluster)
                          )
  
  cluster_dt$pathString = paste(path_tree, 1:nbr_k_optimal,sep = "/")
  
  
  return(cluster_dt)
}
#### with pharmacophores ####
dt_pharmacophores[1:5,1:1000]
length(rownames(dt_pharmacophores))
length(rownames(dt_12descriptors))
names_pharmacophores_descriptors = intersect(rownames(dt_pharmacophores), rownames(dt))
length(names_pharmacophores_descriptors)

dt_pharmacophores_consensus = dt_pharmacophores[1:5,1:1000]
pharmacophores_consensus_mean = apply(dt_pharmacophores[,], 2, function(x) v <- as.integer(mean(x)))
pharmacophores_consensus_max = apply(dt_pharmacophores[,], 2, function(x) v <- as.integer(mean(x)))

pharmacophores_consensus_50 = apply(dt_pharmacophores[,], 2, function(x) {
  if(length(which(x>=1)) >= length(x)/2) {
    return(as.integer(mean(x)))
  } else {
    return(0)
  }
  })
length(which(pharmacophores_consensus_50>= 1))
length(which(dt_pharmacophores[6,]>= 1))

length(which(dt_pharmacophores[,100]>=1))
length(dt_pharmacophores[,100])/2

dt_pharmacophores[1:5,485]
pharmacophores_consensus_mean[485]
pharmacophores_consensus_sum[927]
dt_ph = c(0,2,0,1,0,1)
sum(dt_ph >= 1)

apply(dt_pharmacophores[intersect(rownames(dt_pharmacophores),
                                  names(which(dt.kmean$cluster == 1))),], 2,
                                  function(x) v <- as.integer(mean(x)))
length(intersect(rownames(dt_pharmacophores), names(which(dt.kmean$cluster == 1))))
length(names(which(dt.kmean$cluster == 1)))
length(intersect(rownames(dt_pharmacophores), rownames(dt_12descriptors)))
##
names_physicochemical = c('p_aliphatic_residues','p_Otyr_atom','p_NE2_atom','p_Nlys_atom','p_Ntrp_atom',
                      'p_Ooh_atom','p_ND1_atom')
names_geometrical = c('C_ATOM','RADIUS_CYLINDER','CONVEX.SHAPE_COEFFICIENT','RADIUS_HULL','C_RESIDUES','SURFACE_HULL')
####~ workflow works ~####
#Management dt_72
r_names_dt_72 = toupper(gsub("_prox5_5.desR","",rownames(dt_72descriptors)))
#which(is.element(r_names_dt_72, toupper(rownames(dt))))
dt_72 = dt_72descriptors[which(is.element(r_names_dt_72, toupper(rownames(dt)))),]
rownames(dt_72) = toupper(gsub("_prox5_5.desR","",rownames(dt_72)))
#length(intersect(rownames(dt_72), toupper(rownames(dt))))

#dt = dt_72descriptors[,]
dt = dt_12descriptors[intersect(rownames(dt_pharmacophores), rownames(dt_12descriptors)),]
dt = dt_12descriptors[,]#names_pharmacophores_descriptors,]#c(-8,-2)]
#dt = round(dt,digits=4)
#dt[which(dt[,8] == 0),8] = 0.000000001 #to avoid issue with the k means
#dt[which(dt[,2] == 0),2] = 0.000000001
summary(dt)
dt = delete_clean_data(dt)
#dt = dt[,names_physicochemical]
#dt = dt[,names_geometrical]
#mix data
#dt = dt[sample(1:nrow(dt)),] 
#scale data
dt = scale(dt)
#
pockets_cluster_names = list(row.names(dt))
path_tree = c("alltree")
list_path_tree = NULL
list_path_tree = path_tree
tmp_pockets_cluster_names = NULL
tmp_list_path_tree = NULL
cluster_infos = NULL
nstart = 1#100
n = 1
iter=1
iter_path = 1
for (iter in 1:1) {
  for(i in 1:length(pockets_cluster_names)) {
    if(!is.null(pockets_cluster_names[[i]])) {
      cluster_dt = pockets_classification_tree(
                                      dt = dt[unlist(pockets_cluster_names[[i]]),],
                                      nstart = nstart,
                                      path_tree = list_path_tree[iter_path],
                                      dt_pharmacophores = dt_pharmacophores)
      cluster_infos = rbind(cluster_infos, cluster_dt)
      print("hei")
      print(list_path_tree[iter_path])
      
      for (j in 1:nrow(cluster_dt)) {
        if(cluster_dt[j,"size"] > 20) {
          tmp_pockets_cluster_names[[n]] = unlist(cluster_dt[j,"pockets_names"])
          tmp_list_path_tree = c(tmp_list_path_tree, paste(list_path_tree[iter_path], j, sep ="/"))
        } #TODO:add hierachical clustering
        n=n+1
      }
      iter_path = iter_path + 1
    }
  }
  pockets_cluster_names = tmp_pockets_cluster_names
  list_path_tree = tmp_list_path_tree
  iter_path = 1
  n = 1
  tmp_pockets_cluster_names = NULL
  tmp_list_path_tree = NULL
  print(iter)
}
cluster_alltree = data.frame(
                        withinss = NA,
                        size = NA,
                        betweenss = NA,
                        totss = NA,
                        pockets_names = NA,
                        cluster_dend = I(list(
                                       as.dendrogram(hclust(dist(dt.kmean$centers),
                                                     method = "ward.D2")))),
                        pathString = "alltree"
)
cluster_alltree = cbind(cluster_infos[1,1:12],cluster_alltree)
cluster_alltree[,1:12] = NA
colnames(cluster_alltree)
cluster_infos = rbind(cluster_alltree,cluster_infos)
alltree <- as.Node(cluster_infos)
print(alltree, "size","withinss", "totss", "betweenss")
print(alltree, "size","withinss", "centers.hydrophobic_kyte", "centers.p_Nlys_atom", "centers.p_Ntrp_atom")
alltree$fieldsAll
##save tree
save(alltree, file = "../results/K_Means_comparaison/dt_12clean_tree_classic_seeds800_dend.Rdata", version=2)
#save scale infos
write.table(attr(dt, "scaled:center"), file = "../results/scaled:center_dt12clean.Rdata")
write.table(attr(dt, "scaled:scale"), file = "../results/scaled:scale_dt12clean.Rdata")
#save min max
write.table(rbind(apply(dt,2,max),apply(dt,2,min)), file = "../results/minmax_value_dt12clean.Rdata")
#
load(file = "../results/dt_12dsc_tree_mcqueen_1000.Rdata")
load(file = "../results/res_12desc/dt_12dsc_tree_size400_mcqueen.Rdata")

load(file = "../results/res_12desc/dt_12dsc_tree_size1000_hart.Rdata")

load(file = "../results/res_12desc/dt_12dsc_tree_withinss400_mcqueen.Rdata")

load(file = "../results/res_12desc/dt_12dsc_tree_withinss400_mcqueen_pharm.Rdata")
print(alltree, "size","withinss")

for (i in 1:length(alltree$children)) {
  print(alltree$children[[i]]$size)
  print(alltree$children[[i]]$withinss)
}
S_mq = NULL
W_mq = NULL
for (i in 1:length(alltree$children)) {
  S_mq = c(S_mq,alltree$children[[i]]$size)
  W_mq = c(W_mq,alltree$children[[i]]$withinss)
}
plot(W_mq,S_mq)

#### Compute distance of a new value in the tree ####
new_pocket = data.frame(centers.p_polar_residues = alltree$`1`$`1`$centers.p_polar_residues,
                        centers.p_Nlys_atom = alltree$`1`$`1`$centers.p_Nlys_atom,
                        centers.p_aliphatic_residues = alltree$`1`$`1`$centers.p_aliphatic_residues,
                        centers.VOLUME_HULL = alltree$`1`$`1`$centers.VOLUME_HULL,
                        centers.DIAMETER_HULL = alltree$`1`$`1`$centers.DIAMETER_HULL,
                        centers.p_Ooh_atom = alltree$`1`$`1`$centers.p_Ooh_atom,
                        centers.hydrophobic_kyte = alltree$`1`$`1`$centers.hydrophobic_kyte,
                        centers.p_Ntrp_atom = alltree$`1`$`1`$centers.p_Ntrp_atom,
                        centers.p_Otyr_atom = alltree$`1`$`1`$centers.p_Otyr_atom,
                        centers.C_RESIDUES = alltree$`1`$`1`$centers.C_RESIDUES,
                        centers.p_aromatic_residues = alltree$`1`$`1`$centers.p_aromatic_residues,
                        centers.p_hydrophobic_residues = alltree$`1`$`1`$centers.p_hydrophobic_residues
                        )
#pocket NS1 conf0101_p0
dt_ns1_conf0101_p0 = read.table("../data/pockets_MD_NS1/pocket_PPE_conf0101_p0.txt", header = T, sep = "\t", row.names = 1, fill=TRUE)
dt_ns1_conf0101_p0 = dt_ns1_conf0101_p0["pocket1_atm",]

dt_ns1_conf0101_p0_test = read.table("../data/pockets_MD_NS1/test/pocket0_atm.desR", header = T, sep = "\t")
dt_ns1_conf0101_p0_test_nat = read.table("../data/pockets_MD_NS1/test/pred_user_results_pock_pocket0_atm.desR_NA.csv", header = T, sep = ",", row.names = 1)
dt_ns1_conf0101_p0_test_nat = na.omit(dt_ns1_conf0101_p0_test_nat)
dt_ns1_conf0101_p0_test_nath = dt_ns1_conf0101_p0_test_nat[,2]
names(dt_ns1_conf0101_p0_test_nath) = rownames(dt_ns1_conf0101_p0_test_nat)


new_pocket = data.frame(centers.p_polar_residues = dt_ns1_conf0101_p0[1,"Polar.residues"],
                        #centers.p_Nlys_atom = alltree$`1`$`1`$centers.p_Nlys_atom,
                        centers.p_aliphatic_residues = dt_ns1_conf0101_p0[1,"Aliphatic.residues"],
                        centers.VOLUME_HULL = dt_ns1_conf0101_p0[1,"Volume.hull"],
                        centers.DIAMETER_HULL = dt_ns1_conf0101_p0[1,"Diameter.hull"],
                        centers.p_Ooh_atom = dt_ns1_conf0101_p0[1,"Ooh.atom"],
                        centers.hydrophobic_kyte = dt_ns1_conf0101_p0[1,"Hydrophobic.kyte"],
                        #centers.p_Ntrp_atom = alltree$`1`$`1`$centers.p_Ntrp_atom,
                        centers.p_Otyr_atom = dt_ns1_conf0101_p0[1,"Otyr.atom"],
                        centers.C_RESIDUES = dt_ns1_conf0101_p0[1,"Nb.RES"],
                        centers.p_aromatic_residues = dt_ns1_conf0101_p0[1,"Aromatic.residues"],
                        centers.p_hydrophobic_residues = dt_ns1_conf0101_p0[1,"Hydrophobic.residues"]
)

new_pocket = data.frame(centers.p_polar_residues = dt_ns1_conf0101_p0[1,"p_polar_residues"],
                        centers.p_Nlys_atom = dt_ns1_conf0101_p0[1,"p_Nlys_atom"],
                        centers.p_aliphatic_residues = dt_ns1_conf0101_p0[1,"p_aliphatic_residues"],
                        centers.VOLUME_HULL = dt_ns1_conf0101_p0[1,"VOLUME_HULL"],
                        centers.DIAMETER_HULL = dt_ns1_conf0101_p0[1,"DIAMETER_HULL"],
                        centers.p_Ooh_atom = dt_ns1_conf0101_p0[1,"p_Ooh_atom"],
                        centers.hydrophobic_kyte = dt_ns1_conf0101_p0[1,"hydrophobic_kyte"],
                        centers.p_Ntrp_atom = dt_ns1_conf0101_p0[1,"p_Ntrp_atom"],
                        centers.p_Otyr_atom = dt_ns1_conf0101_p0[1,"p_Otyr_atom"],
                        centers.C_RESIDUES = dt_ns1_conf0101_p0[1,"C_RESIDUES"],
                        centers.p_aromatic_residues = dt_ns1_conf0101_p0[1,"p_aromatic_residues"],
                        centers.p_hydrophobic_residues = dt_ns1_conf0101_p0[1,"p_hydrophobic_residues"]
)
#IMPORTANT scale the data
new_test = dt_12descriptors["4UQX_ACT_A_1",c(-2,-8)]
new_test = scale(new_test, attr(dt, "scaled:center"), attr(dt, "scaled:scale"))
new_pocket = scale(new_pocket, attr(dt, "scaled:center"), attr(dt, "scaled:scale"))

alltree$Do(function(node) node$dist <- dist(rbind(c(
                                                    node$centers.p_polar_residues,
                                                    node$centers.p_Nlys_atom,
                                                    node$centers.p_aliphatic_residues,
                                                    node$centers.VOLUME_HULL,
                                                    node$centers.DIAMETER_HULL,
                                                    node$centers.p_Ooh_atom,
                                                    node$centers.hydrophobic_kyte,
                                                    node$centers.p_Ntrp_atom,
                                                    node$centers.p_Otyr_atom,
                                                    node$centers.C_RESIDUES,
                                                    node$centers.p_aromatic_residues,
                                                    node$centers.p_hydrophobic_residues
                                                    ),
                                                    new_pocket)))
dist(rbind(c(alltree$`1`$`1`$centers.p_polar_residues,
             alltree$`1`$`1`$centers.p_Nlys_atom,
             alltree$`1`$`1`$centers.p_aliphatic_residues,
             alltree$`1`$`1`$centers.VOLUME_HULL,
             alltree$`1`$`1`$centers.DIAMETER_HULL,
             alltree$`1`$`1`$centers.p_Ooh_atom,
             alltree$`1`$`1`$centers.hydrophobic_kyte,
             alltree$`1`$`1`$centers.p_Ntrp_atom,
             alltree$`1`$`1`$centers.p_Otyr_atom,
             alltree$`1`$`1`$centers.C_RESIDUES,
             alltree$`1`$`1`$centers.p_aromatic_residues,
             alltree$`1`$`1`$centers.p_hydrophobic_residues
             ),new_pocket))

dist(rbind(1:10,1:10))
alltree$`1`$`1`$`1`$pockets_names

print(alltree, "size", "withinss", "dist")
for (i in 1:length(alltree$children)) {
  print(alltree$children[[i]], "size", "dist")
}
Sort(alltree, "dist", decreasing = FALSE)
####Found closest pockets ####
name_clust_close_1 = unlist(alltree$`1`$`1`$`4`$pockets_names)
name_clust_close_2 = unlist(alltree$`1`$`1`$`9`$pockets_names)

name_clust_close_test = c(unlist(alltree$`5`$`5`$`10`$pockets_names))#, unlist(alltree$`3`$pockets_names))
name_clust_close_test = c(unlist(alltree$`2`$`8`$pockets_names), unlist(alltree$`2`$`9`$pockets_names))

dist_clust_close_1 = NULL
dist_clust_close_2 = NULL
dist_clust_close_test = NULL
for (name_pocket in rownames(dt[name_clust_close_1,])) {
  dist_clust_close_1 = c(dist_clust_close_1, dist(rbind(dt[name_pocket,],new_pocket)))
}
for (name_pocket in rownames(dt[name_clust_close_2,])) {
  dist_clust_close_2 = c(dist_clust_close_2, dist(rbind(dt[name_pocket,],new_pocket)))
}
for (name_pocket in rownames(dt[name_clust_close_test,])) {
  dist_clust_close_test = c(dist_clust_close_test, dist(rbind(dt[name_pocket,],new_pocket)))
}
names(dist_clust_close_1) = rownames(dt[name_clust_close_1,])
names(dist_clust_close_2) = rownames(dt[name_clust_close_2,])
names(dist_clust_close_test) = rownames(dt[name_clust_close_test,])

min(dist_clust_close_1)
min(dist_clust_close_2)

dist(rbind(dt["1GTV_TMP_B_1",], new_pocket))
intersect("1GTV_TMP_B_1", names(dist_clust_close_test))
dt_12descriptors["1GTV_TMP_B_1",]

dist_clust_close_1[which(dist_clust_close_1 < 3.57)]
dist_clust_close_2[which(dist_clust_close_2 < 3.57)]
names(which(dist_clust_close_1 < 3.57))
names(which(dist_clust_close_2 < 3.57))
#on test
intersect(names(dt_ns1_conf0101_p0_test_nath), names(dist_clust_close_test))
dt_ns1_conf0101_p0_test_nath[names(dist_clust_close_test)]

dist_clust_close_test["1GTV_TMP_B_1"]
dist_clust_close_test["3FVF_1JZ_B_1"]
min(dist_clust_close_test)

length(dist_clust_close_test)
length(names(which(dt_ns1_conf0101_p0_test_nath < 3)))
length(intersect(names(which(dt_ns1_conf0101_p0_test_nath < 3)), names(dist_clust_close_test)))

name_intersect_test = intersect(names(which(dt_ns1_conf0101_p0_test_nath < 2)), names(dist_clust_close_test))
unique(names(dist_clust_close_test))

which(dist_clust_close_test == min(dist_clust_close_test))

length(which(dist_clust_close_test < 2.5))
dt_ns1_test_merge = merge(dt_ns1_conf0101_p0_test_nath, dist_clust_close_test, by = intersect(names(dt_ns1_conf0101_p0_test_nath), names(dist_clust_close_test)))
dt_ns1_test_merge = data.frame(nath = dt_ns1_conf0101_p0_test_nath[names(dist_clust_close_test)], me = dist_clust_close_test, row.names = names(dist_clust_close_test))

####search by pharmacophores ####
v_ph_1 = unlist(alltree$`1`$`1`$pharmacophores_consensus_mean)
v_ph_2 = unlist(alltree$`2`$pharmacophores_consensus_mean)
length(v_ph_1)

alltree$Do(function(node) {
  v_ph_2 = unlist(node$pharmacophores_consensus_mean)
  if(length(v_ph_1) == length(v_ph_2)){
    t = table(v_ph_1, v_ph_2)
    print(node$path)
    print(v_ph_2)
    print(ncol(t))
    common_non_null = sum(t[2:nrow(t),2:ncol(t)])
    non_null_counts_Fa = sum(t[2:nrow(t),])
    non_null_counts_Fb = sum(t[,2:ncol(t)])
    node$val_similarity = common_non_null/min(non_null_counts_Fa,non_null_counts_Fb)
  } else {
    node$val_similarity = 0
  }
})
alltree$`1`$centers.C_RESIDUES
print(alltree$`1`, "size", "dist", "val_similarity", "centers.C_RESIDUES")
print(alltree$`8`, "size", "dist", "val_similarity", "centers.C_RESIDUES")

print(alltree, "size", "dist")#, "val_similarity")
alltree$Get("val_similarity", filterFun = function(node) node$level == 2)
alltree$Get("dist", filterFun = function(node) node$level == 2)
Sort(alltree, "val_similarity", decreasing = TRUE)    

####apply dist ####

start_time <- Sys.time()
d = as.matrix(dt[unlist(alltree$`4`$`8`$`2`$pockets_names),])
res_d = apply(d, 1, function(x) {dist(rbind(d[1,],x))})
end_time <- Sys.time()

end_time - start_time

min(res_d)
which(res_d == min(res_d))
res_d["4X49_G39_A_1"]

#### Get access to informations in the tree ####
alltree$fieldsAll
#ntotal clusters
alltree$totalCount
#sort
Sort(alltree, "size", decreasing = FALSE)
#max
maxSize <- Aggregate(alltree, "size", max)
alltree$Get("name", filterFun = function(x) x$isLeaf && x$size == maxSize)
#plot according Iintra
#plot(as.dendrogram(alltree, heightAttribute = "withinss"))
plot(alltree$`1`)

alltree$`1`$isLeaf
print(alltree$children[[5]]$`2`$`10`,"size","withinss")
print(alltree, "size", "withinss")
length(alltree$children[[1]])

####Save infos withinss and size####

####LOOK  which cluster have a dist < 1####
for (i in 1:length(alltree$children)) {
  if(alltree$children[[i]]$dist < 1) {
    print("level 1")
    print(i)
  }
  for (j in 1:length(alltree$children[[i]]$children)) {
    if(alltree$children[[i]]$children[[j]]$dist < 1) {
      print("level 2")
      print(i)
      print(j)
    }
  }
}
alltree$averageBranchingFactor

length(alltree$children)
alltree$Do()

#test for distance calculation
data(acme)
acme$fieldsAll
print(acme, "cost","p")
new = data.frame(cost = 30000, p = 0.50)
new_bis = data.frame(cost = 40000, p = 0.60)
new = c(51000,0.9)
new_bis = c(40000,0.6)
dist(rbind(new,new_bis))

acme$Do(function(node) node$dist <- dist(rbind(c(node$p,node$cost),new)))
print(acme, "expectedCost", "p", "cost", "dist")
###
alltree$fieldsAll
#to have access to vthe values:
107+87+95+54+96+151+42+137+117+115
print(alltree$`2`$centers.VOLUME_HULL)
alltree$`1`$`1`$`1`$pockets_names
alltree$`1`$`1`$pockets_names

intersect(unlist(alltree$`1`$`1`$`1`$pockets_names), unlist(alltree$`1`$`1`$pockets_names))

#
v=0
for(i in 1:length(pockets_cluster_names)) {
  if(!is.null(pockets_cluster_names[[i]])) {
    v = v+1
    print(length(pockets_cluster_names[[i]]))
  }
}
v


#tester :mettre dans infps le n iteration et le numero de cluster

#
library(pvclust)
library(MASS)
dt.pv <- pvclust(t(scale(dt.kmean$centers)))
plot(dt.pv)
#
data(acme)
print(acme, prob = acme$Get("p", format = function(x) FormatFixedDecimal(x, 4)))
print(acme)
acme$Climb()
Climb(acme, position = c(2))
Climb(acme, name = 'IT')

itClone <- Clone(acme$IT)

dev.off()
dt.hclust = dt_centers.hclust
#alltree <- as.Node(dt.hclust)
dend1 <- as.dendrogram(dt.hclust)
tree <- as.Node(dend1)
tree
tree $fieldsAll
tree $totalCount
tree $leafCount
tree $height
#plot(alltree)
library(DiagrammR)

data(GNI2014)
head(GNI2014)
GNI2014$pathString <- paste("world", 
                            GNI2014$continent, 
                            GNI2014$country, 
                            sep = "/")

population <- as.Node(GNI2014)
print(population, "iso3", "population", "GNI", limit = 20)

#test kmean
dt = dt_12descriptors[,]
dt = delete_clean_data(dt)
dt[,1] = 0.00000000
summary(dt)
dt.kmean = kmeans(scale(dt), 10, nstart = 1)

#grep
dt_12descriptors[grep("2NNL",rownames(dt_12descriptors)),]
if (grep("9999",rownames(dt_12descriptors))) {
  print("r")
}
k=1

is.null(grep("2NNL_ERD", unlist(node$centers.p_polar_residues)))
n_test = function(node){
  print(node)
  print(node$centers.p_polar_residues)
  if (node$centers.p_polar_residues > 0) { 
      print(node$centers.p_polar_residues)
      erturn(1) 
  } else {
      return(0)
  }
}
alltree$Do(function(node) {
  node$test = 0
  print(node$path)
  print(grep("5UPJ", unlist(node$pockets_names)))
  print(grep("6UPJ", unlist(node$pockets_names)))
  print(grep("1MU2", unlist(node$pockets_names)))
  #print(grep("3LM5_QUE", unlist(node$pockets_names)))
  #print(grep("4HKN_LU2", unlist(node$pockets_names)))
  #print(grep("4DGM_AGI", unlist(node$pockets_names)))
},
filterFun = function(node) node$level == 2
)

print(alltree$`8`, "test", "size")
length(alltree$Get("size", filterFun = function(node) node$level == 5))

alltree$Get("path", filterFun = function(node) node$level == 2)
grep("1MU2",rownames(dt_12descriptors))

####check inetie intra, inter, R2 (inter/total), intra/total ####
###intra
##1 er decoupage
#intra
alltree$Get("withinss", filterFun = function(node) node$level == 2)
sum_intra = sum(alltree$Get("withinss", filterFun = function(node) node$level == 2))
#inter
tots = max(alltree$Get("totss", filterFun = function(node) node$level == 2))
inter = tots - sum_intra
tots
inter
#R2
inter/tots
#2eme decoupage
alltree$Get("withinss", filterFun = function(node) node$level == 3)
sum_intra = sum(alltree$Get("withinss", filterFun = function(node) node$level == 3))
#inter
tot = alltree$Get("totss", filterFun = function(node) node$level == 3)
tots = sum(tot[which(names(tot) == "1")])
inter = tots - sum_intra
tots
inter
#R2
inter/tots
#3eme decoupage

#feuilles
w_tree = sqrt(alltree$Get("withinss", filterFun = isLeaf))
s_tree = alltree$Get("size", filterFun = isLeaf)
length(w_tree)
plot(s_tree, w_tree)

sum_intra = sum(alltree$Get("withinss", filterFun = isLeaf))
#inter
tot = alltree$Get("totss", filterFun = isLeaf)
length(tot)
tots = sum(tot[which(names(tot) == "1")])
inter = tots - sum_intra
tots
inter
#R2
inter/tots
alltree$Do(function(node) node$Rws <- node$withinss/node$size)
print(alltree, "R2", "Rws","withinss", "size")
Wrs = alltree$Get("Rws", filterFun = isLeaf)

w_tree = alltree$Get("withinss", filterFun = isLeaf)
bet = tots - sum(w_tree)
bet/tots ######### here total R2

R2 = alltree$Get("R2", filterFun = function(node) node$level == 4)
# plot

w_tree = alltree$Get("withinss", filterFun = isLeaf)
s_tree = alltree$Get("size", filterFun = isLeaf)

png(paste0(paste0("../results/K_Means_comparaison/distance_size_LP_nseeds_", seeds),
           "_nested.png"))

par(mar=c(5, 4, 4, 8), xpd=TRUE)
plot(s_tree, sqrt(w_tree/s_tree), xlim = c(0,300), ylim = c(0,3.5),
     main = paste("distance moyenne NESTED n_seeds = ", seeds))

grp_sup_2 = which(sqrt(w_tree/s_tree) >=2)

nbr_unique_lig = NULL
nbr_unique_prot = NULL
for (n_grp in 1:seeds) {
  print(n_grp)
  names_grp = unlist(alltree$Get("pockets_names", filterFun = isLeaf)[n_grp])
  nbr_unique_lig = c(nbr_unique_lig,
                     length(unique(sapply(strsplit(names_grp, "_"), "[", 2))))
  nbr_unique_prot = c(nbr_unique_prot,
                      length(unique(sapply(strsplit(names_grp, "_"), "[", 1))))
}
nbr_unique_prot_pch = nbr_unique_prot
nbr_unique_prot_pch[which(nbr_unique_prot == 1)] = 13
nbr_unique_prot_pch[which(nbr_unique_prot > 1)] = 10
nbr_unique_prot_pch[which(nbr_unique_prot > 10)] = 1
#points(mean(dt.kmean$size), mean(sqrt(dt.kmean$withinss/dt.kmean$size)), col = "red")

points(s_tree[which(nbr_unique_lig == 1)], 
       sqrt(w_tree/s_tree)[which(nbr_unique_lig == 1)], 
       col = "red", pch = nbr_unique_prot_pch[which(nbr_unique_lig == 1)])
points(s_tree[which(nbr_unique_lig > 1)], 
       sqrt(w_tree/s_tree)[which(nbr_unique_lig > 1)], 
       col = "orange", pch = nbr_unique_prot_pch[which(nbr_unique_lig > 1)])
points(s_tree[which(nbr_unique_lig > 10)], 
       sqrt(w_tree/s_tree)[which(nbr_unique_lig > 10)], 
       col = "green", pch = nbr_unique_prot_pch[which(nbr_unique_lig > 10)])
#text(s_tree, sqrt(w_tree/s_tree), 
#     labels=paste(paste0("L:",nbr_unique_lig),
#                  paste0("P:",nbr_unique_prot), sep = "|"),
#     cex= 0.6, pos=3, col = "blue")

#text(dt.kmean$size[grp_sup_2], sqrt(dt.kmean$withinss/dt.kmean$size)[grp_sup_2], 
#     labels=paste(dt.kmean$size[grp_sup_2],
#                  round(sqrt(dt.kmean$withinss/dt.kmean$size)[grp_sup_2], 2), sep = ","),
#     cex= 0.6, pos=3, col = "blue")
legend("topright", inset=c(-0.35,0), xpd=TRUE, mar(c(7,7,7,7)), cex = 1, bty = "n",
       legend=c(paste0("mean_dist:"),
                paste0(round(mean(sqrt(w_tree/s_tree)),2),
                       paste0("±", 
                              paste0(round(sd(sqrt(w_tree/s_tree)),2)),"sd")),
                paste0("mean_size:",round(mean(mean(s_tree)),2)),
                paste0("n_dist>=2:",length(which(sqrt(w_tree/s_tree) >=2))),
                paste0("n_pock>=2:",sum(s_tree[which(sqrt(w_tree/s_tree) >=2)])),
                paste0("grp_size=1:",length(which(s_tree == 1))),
                paste0("grp_sizeL=1:",length(which(nbr_unique_lig == 1))),
                paste0("grp_sizeP=1:",length(which(nbr_unique_prot == 1))) ))
legend( x="bottomright", inset=c(-0.375,0), xpd=TRUE, mar(c(7,7,7,7)), cex = 1, bty = "n",
        legend=c("grp_sizeL=[1]", 
                 "grp_sizeL=]1,10]",
                 "grp_sizeL=]10,[",
                 "grp_sizeP=[1]",
                 "grp_sizeP=]1,10]",
                 "grp_sizeP=]10,["),
        col=c("red","orange","green", "black", "black", "black"),
        lwd=1, lty=c(0,0), pch=c(1,1,1,13,10,1), x.intersp = 0)
dev.off()

png(paste0(paste0("../results/K_Means_comparaison/size_L_nseeds_", seeds),
           "_nested.png"))
plot(nbr_unique_lig, s_tree, xlim = c(0,300), ylim = c(0,300) , 
     xlab = "Nombre ligands dans groupe", ylab = "Taille du groupe",
     main = paste("taille clusters et nbr ligands associés NESTED n_seeds = ", seeds))
dev.off()
png(paste0(paste0("../results/K_Means_comparaison/size_P_nseeds_", seeds),
           "_nested.png"))
plot(nbr_unique_prot, s_tree, xlim = c(0,300), ylim = c(0,300) , 
     xlab = "Nombre proteines dans groupe", ylab = "Taille du groupe",
     main = paste("taille clusters et nbr proteines associés NESTED n_seeds = ", seeds))
dev.off()

### check bad ligands ###
for (i in 1:alltree$averageBranchingFactor) {
  print(length(grep("DOD", unlist(alltree$children[[i]]$pockets_names))))
}
length(grep("DOD", rownames(dt)))

length(rownames(dt[grep("_DOD_", rownames(dt)),]))
length(rownames(dt[grep("_IOD_", rownames(dt)),]))
length(rownames(dt[grep("_IUM_", rownames(dt)),]))
length(rownames(dt[grep("_D8U_", rownames(dt)),]))
length(rownames(dt[grep("_SO4_", rownames(dt)),]))
length(rownames(dt[grep("_FE2_", rownames(dt)),]))
length(rownames(dt[grep("_GOL_", rownames(dt)),]))
length(rownames(dt[grep("_FES_", rownames(dt)),]))
length(rownames(dt[grep("_PO4_", rownames(dt)),]))
length(rownames(dt[grep("_MSE_", rownames(dt)),]))
length(rownames(dt[grep("_DMS_", rownames(dt)),]))
length(rownames(dt[grep("_URE_", rownames(dt)),]))
length(rownames(dt[grep("_FMT_", rownames(dt)),]))
length(rownames(dt[grep("_TRS_", rownames(dt)),]))
length(rownames(dt[grep("_NCO_", rownames(dt)),]))
length(rownames(dt[grep("_HOH_", rownames(dt)),]))
length(rownames(dt[grep("_H2O_", rownames(dt)),]))
length(rownames(dt[grep("_WAT_", rownames(dt)),]))
length(rownames(dt[grep("_NTN_", rownames(dt)),]))


###check withinss###
dt = dt_12descriptors
dt = delete_clean_data(dt)
dt = scale(dt[1:5,])
dt.kmean = kmeans(dt, 2, nstart = 10)
dt.kmean$centers
dt.kmean$withinss
dt.kmean$cluster

x = dt[dt.kmean$cluster == 2,]

ss <- function(x) sum(apply(x, 2, function(x) x - mean(x))^2)
ss(x)

#kmeans new
dt.kmean = kmeans(dt, 100, nstart = 1, algorithm="MacQueen", iter.max = 200)
#plot withinss/size
plot(dt.kmean$size, dt.kmean$withinss)
boxplot(dt.kmean$size)
boxplot(dt.kmean$withinss)

alltree$Do(function(node) node$sd <- sd(c(
  node$centers.p_polar_residues,
  node$centers.p_Nlys_atom,
  node$centers.p_aliphatic_residues,
  node$centers.VOLUME_HULL,
  node$centers.DIAMETER_HULL,
  node$centers.p_Ooh_atom,
  node$centers.hydrophobic_kyte,
  node$centers.p_Ntrp_atom,
  node$centers.p_Otyr_atom,
  node$centers.C_RESIDUES,
  node$centers.p_aromatic_residues,
  node$centers.p_hydrophobic_residues
  )))

alltree$`1`$`1`$`1`$sd

###check withins dist with means and sd ###
dist(dt[unlist(alltree$`1`$`1`$`1`$pockets_names),])
alltree$Do(function(node) node$node_intra_dist <- dist(dt[unlist(node$pockets_names),]), filterFun = isLeaf)

#mean at each leaf
mean_node_intra_dist = alltree$Get(function(node) mean(node$node_intra_dist), filterFun = isLeaf)
sd_node_intra_dist = alltree$Get(function(node) sd(node$node_intra_dist), filterFun = isLeaf)
size_cluster = alltree$Get("size", filterFun = isLeaf)

mean_node_intra_dist
sd_node_intra_dist
plot(mean_node_intra_dist,sd_node_intra_dist)
length(mean_node_intra_dist)

plot(size_cluster, mean_node_intra_dist)

length(grep("_COC_", rownames(dt)))
rownames(dt[c(99702,91127),])
#within | average within
average_node_within = alltree$Get(function(node) sqrt(node$withinss/node$size), filterFun = function(node) node$level == 4)
node_within = alltree$Get("withinss", filterFun = function(node) node$level == 4)
node_size = alltree$Get("size", filterFun = function(node) node$level == 4)
plot(node_size,average_node_within, main = "NESTED")  

length(node_within)
#in one splitting
dt.kmean_10 = kmeans(dt, 1000, nstart = 1, algorithm="MacQueen", iter.max = 200)

mean_node_intra_dist = NULL
sd_node_intra_dist = NULL
size_cluster = NULL
for (i in 1:1000) {
  mean_node_intra_dist = c(mean_node_intra_dist,mean(dist(dt[which(dt.kmean$cluster == i),])))
  sd_node_intra_dist = c(sd_node_intra_dist, sd(dist(dt[which(dt.kmean$cluster == i),])))
  size_cluster = c(size_cluster, length(which(dt.kmean$cluster == i)))
}
intersect(dt.kmean$size, size_cluster)

plot(dt.kmean_10$size,sqrt(dt.kmean_10$withinss/dt.kmean_10$size),  xlim = c(0,550), ylim = c(0,3), main = "CLASSIC") 
# box plot of both for 2 methods with 1000 clusters

#boxplot average withins
boxplot(cbind(sqrt(dt.kmean$withinss/dt.kmean$size),as.vector( average_node_within)), names = c("k means wtih 1000 seeds","successive k means"), main = "average withinss cluster", ylab = "average withins")
boxplot(sqrt(dt.kmean$withinss/dt.kmean$size))

#boxplot size cluster
boxplot(cbind(dt.kmean$size, node_size), names = c("k means wtih 1000 seeds","successive k means"), main = "size clusters")

#boxplot initial within cluster
boxplot(cbind(dt.kmean$withinss, node_within), names = c("k means wtih 1000 seeds","successive k means"), main = "withinss clusters")

#table 
tree_cluster_leaf = rep(0,length(dt.kmean_10$cluster))
names(tree_cluster_leaf) = names(dt.kmean_10$cluster)
names_cluster_leaf = alltree$Get("pockets_names", filterFun = function(node) node$level == 3)
for (i in 1:length(names_cluster_leaf)) {
  tree_cluster_leaf[unlist(names_cluster_leaf[i])] = i
}
table_cluster = table(dt.kmean_10$cluster,tree_cluster_leaf)

which(table_cluster[1,] > 0)

boxplot(table_cluster[1,])

plot(apply(table_cluster, 2, max))

max(dt.kmean$cluster)


matrix.sort <- function(matrix) {
  
  if (nrow(matrix) != ncol(matrix)) stop("Not diagonal")
  if(is.null(rownames(matrix))) rownames(matrix) <- 1:nrow(matrix)
  
  row.max <- apply(matrix,1,which.max)
  if(all(table(row.max) != 1)) stop("Ties cannot be resolved")
  
  return(matrix[names(sort(row.max)),])
}
table_cluster = matrix.sort(scale(table_cluster))

#library(pheatmap)
pheatmap(table_cluster, cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap(table_cluster[1:100,1:100], cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap(table_cluster[100:200,100:200], cluster_rows = FALSE, cluster_cols = FALSE)

#### ligand protein diversity ####
rownames(dt)[30065]
grep("_CQL_",rownames(dt))

names_prot = sapply(strsplit(rownames(dt), "_"), "[", 1)
names_ligand = sapply(strsplit(rownames(dt), "_"), "[", 2)
length(which(names_ligand == "CLQ"))

length(which(table(names_ligand) > 1))

length(unique(names_ligand))
#names_ligand = unique(names_ligand)

table_names_ligand = table(names_ligand)
which(table_names_ligand >= 100)

plot(sort(table_names_ligand))

hist(sort(table_names_ligand))

table_names_ligand[which(table_names_ligand > 500)]
length(which(table_names_ligand == 1000))

#camember plot:
df <- data.frame(
  group = c("linked to 1 pocket", "linked ]1;10] pockets", "linked to ]10;1000[", "linked to > 1000 pockets"),
  value = c(sum(table_names_ligand[which(table_names_ligand == 1)])/sum(table_names_ligand),
            sum(table_names_ligand[which(table_names_ligand <= 10 & table_names_ligand > 1)])/sum(table_names_ligand),
            sum(table_names_ligand[which(table_names_ligand <= 1000 & table_names_ligand > 10)])/sum(table_names_ligand),
            sum(table_names_ligand[which(table_names_ligand > 1000)])/sum(table_names_ligand)
            )
)
df <- data.frame(
  group = c("linked to 1 pocket:11718 ligands", "linked to ]1;10] pockets:5208 ligands", "linked to ]10;1000]:690 ligands", "linked to > 1000 pockets:6 ligands"),
  value = c(sum(table_names_ligand[which(table_names_ligand == 1)]),
            sum(table_names_ligand[which(table_names_ligand <= 10 & table_names_ligand > 1)]),
            sum(table_names_ligand[which(table_names_ligand <= 1000 & table_names_ligand > 10)]),
            sum(table_names_ligand[which(table_names_ligand > 1000)])
  ),
  size = c(length(table_names_ligand[which(table_names_ligand == 1)]),
            length(table_names_ligand[which(table_names_ligand <= 10 & table_names_ligand > 1)]),
            length(table_names_ligand[which(table_names_ligand <= 1000 & table_names_ligand > 10)]),
            length(table_names_ligand[which(table_names_ligand > 1000)])
  )
)
library(ggplot2)
# Bar plot
bp<- ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")+
  geom_text(aes(y = value/3 + c(0, cumsum(value)[-length(value)]), 
                label = c(11718,6,690,5208)), size=5)

bp + theme(legend.text = element_text(colour="black", size=12, 
                                   face="bold"))


plot(sort(table_names_ligand), breaks = 4)

barplot(df$value,df$size, names.arg = c("[1]","]1;10]","]10;1000]","]1000;["), xlab = "number of ligands linked to either 1;1-10;10-1000;>1000 pockets",ylab="number of pocketss representative of this group of ligands")

names_prot = sapply(strsplit(rownames(dt), "_"), "[", 1)
length(names_prot)
length(unique(names_prot))
names_prot = unique(names_prot)
grep("1AKE", names_prot)
names_prot[145]

which(names_prot == "3AY6")
alltree$Do(function(node) {
  if(length(grep("_093_",unlist(node$pockets_names))) > 0) {
    print(node$path)
   print(length(grep("_093_",unlist(node$pockets_names))))
  }
},filterFun = isLeaf)

### compute distance and fuzcav for ligand families ###
##normal dist
# dist for same ligands
dist_ligs = rep(0,length(which(table(names_ligand) > 1)))#length(unique(names_ligand)))

names_ligand_unique = unique(names_ligand)
flag = 1
for (i in 1:length(names_ligand_unique)){
  #print(i)
  index = grep(paste0(paste0("_",names_ligand_unique[i]),"_"),rownames(dt))
  if(length(index) > 1) {
    print(i)
    if(length(index) > 100) {
      index = sample(index,100)
    }
    dist_pock = NULL
    for (j in 1:(length(index)-1)) {
      for (k in (j+1):length(index)) {
        if(sapply(strsplit(rownames(dt)[index[j]], "_"), "[", 1) == sapply(strsplit(rownames(dt)[index[k]], "_"), "[", 1)) { #CHANGE TO != FOR DIFFERENT POCKETS
          dist_pock = c(dist_pock, dist(dt[c(index[j],index[k]),]))
          #if(dist(dt[c(index[j],index[k]),]) > 3) {
          #  print(paste(index[j],index[k]))
          #}
          #print(dt[c(index[j],index[k]),2])
          #print(dt[c(index[i],index[j]),])
        }
      }
    }
    if(!is.null(dist_pock)) {
      dist_ligs[flag] = mean(dist_pock)
      flag = flag+1      
    }
  }
}
length(which(dist_ligs > 0))
#check lgands with dist > 3
which(table_names_ligand > 1)[which(dist_ligs > 3)]
#dist_ligs_same_prot = dist_ligs
#test 2
dist_ligs_test2 = rep(0,length(which(table(names_ligand) > 1)))#length(unique(names_ligand)))
names(dist_ligs_test2) =  names(which(table(names_ligand) > 1))
for (i in names(which(table(names_ligand) > 1))){
  #print(i)
  index = grep(paste0(paste0("_",i),"_"),rownames(dt))
  print(i)
  dist_ligs_test2[i] = mean(dist(dt[index,]))
}
dist_ligs_test2

boxplot(dist_ligs_test2)

summary(dist_ligs)
length(dist_ligs)
names(dist_ligs) = names(which(table(names_ligand) > 1))


index = grep(paste0(paste0("_","2CW"),"_"),rownames(dt))
dist(dt[index,])
length(names(which(table(names_ligand) == 2)))
length(which(dist_ligs_test2 == 0))
#comparaison avec prot
names_prot = sapply(strsplit(rownames(dt), "_"), "[", 1)
#check number of pockets dist = 0 for same ligands in same prot like the pocket could be evaluated 2 times
count_same_prot_same_lig_dist0 = NULL
for (i in names(which(table(names_ligand) > 1))){
  #print(i)
  index = grep(paste0(paste0("_",i),"_"),rownames(dt))
  names_prot = sapply(strsplit(rownames(dt[index,]), "_"), "[", 1)
  #print(names_prot)
  #print(setequal(names_prot[1],names_prot))
  if(setequal(names_prot[1],names_prot) == TRUE) {
    print(dt[index,])
    #count_same_prot_same_lig_dist0 = c(count_same_prot_same_lig_dist0,i)
  }
}
length(count_same_prot_same_lig_dist0)
length(names(which(table(names_ligand) > 1)))
setdiff(names(which(dist_ligs_test2 == 0)), count_same_prot_same_lig_dist0)
index = grep(paste0(paste0("_","02N"),"_"),rownames(dt))
dt[index,]
# dist for different ligands
unique_names_ligand = unique(names_ligand)
dist_ligs_random = rep(0,5000)#length(unique(names_ligand)))

for (i in 1:5000){#length(unique(names_ligand))){
  print(i)
  lig = sample(unique_names_ligand, 2)
  dist_ligs_random[i] = dist(rbind(dt[sample(grep(paste0(paste0("_",lig[1]),"_"),rownames(dt)),1),],
                                   dt[sample(grep(paste0(paste0("_",lig[2]),"_"),rownames(dt)),1),]))
}
hist(dist_ligs_random)
d_random = density(dist_ligs_random) # returns the density data
d_lig = density(dist_ligs)#dist_ligs[which(dist_ligs > 0)])
plot(d_random)
plot(d_lig)


#library(sm)

sm.density.compare(c(dist_ligs_random,
                     dist_ligs[which(dist_ligs>0)],
                     dist_ligs_protdiff[which(dist_ligs_protdiff>0)]),
                   c(rep(1,length(dist_ligs_random)),
                     rep(2,length(dist_ligs[which(dist_ligs>0)])),
                     rep(3,length(dist_ligs_protdiff[which(dist_ligs_protdiff>0)]))
                     ),
                   model = "none", xlim=c(0,10)
                   , xlab = "Mean distance bewteen pockets"
                   , main = "Density plot of the distance between pockets from different ligands linking the same ligand")

abline(v=mean(dist_ligs[which(dist_ligs>0)]), col = "green")
abline(v=mean(dist_ligs_protdiff[which(dist_ligs_protdiff>0)]), col = "blue")
abline(v=mean(dist_ligs_random), col = "red")
#abline(v=0.5, col = "green")
#abline(v=2.87, col = "blue")

text(1,0.73,"mean = 1.0",col="green")
#text(1,0.15,"pic = 0.5",col="green")
text(3.5,0.73,"mean = 2.2",col="blue")
text(6,0.73,"mean = 4.6",col="red")
#text(3.5,0.15,"intersect = 2.9",col="blue")
mean(dist_ligs_random)
mean(dist_ligs[which(dist_ligs>0)])
mean(dist_ligs_protdiff[which(dist_ligs_protdiff>0)])

lines(d_lig,col="green")
#SAVE DIST TO GAIN TIME #
save(dist_ligs, file = "../results/dist_pockets/dist_ligs_protsame.Rdata")
save(dist_ligs_protdiff, file = "../results/dist_pockets/dist_ligs_protdiff.Rdata")
save(dist_ligs_random, file = "../results/dist_pockets/dist_ligs_random.Rdata")
#LOAD
load("../results/dist_pockets/dist_ligs_protsame.Rdata")
load("../results/dist_pockets/dist_ligs_protdiff.Rdata")
load("../results/dist_pockets/dist_ligs_random.Rdata")
#
library(sm)
sm.density.compare(d_random, d_lig)

mean(dt[grep("HEM",rownames(dt)),])

mat_hem = dist(dt[grep("HEM",rownames(dt)),])
min(mat_hem)
length(mat_hem)

test_vww = dist(dt[grep("VWW",rownames(dt)),])
####fuzcav dist for pharmacophores ####
dt = as.data.frame(dt_pharmacophores)

# dist for same ligands
dist_ligs = rep(0,length(which(table(names_ligand) > 1)))#length(unique(names_ligand)))
rnames_dt = rownames(dt)

names_ligand_unique = unique(names_ligand)
flag = 1
for (i in 1:length(names_ligand_unique)){
  #print(i)
  index = grep(paste0(paste0("_",names_ligand_unique[i]),"_"),rownames(dt))
  if(length(index) > 1) {
    print(i)
    if(length(index) > 10) {
      index = sample(index,10)
    }
    dist_pock = NULL
    for (j in 1:(length(index)-1)) {
      for (k in (j+1):length(index)) {
        if(sapply(strsplit(rnames_dt[index[j]], "_"), "[", 1) != sapply(strsplit(rnames_dt[index[k]], "_"), "[", 1)) { #CHANGE TO != FOR DIFFERENT POCKETS
          #print(index[j])
          #print(index[k])
          #print(dist_fuzcav_ph(as.integer(dt[index[j],]),as.integer(dt[index[k],])))
          dist_pock = c(dist_pock, 
                        dist_fuzcav_ph(as.integer(dt[index[j],]),as.integer(dt[index[k],])))
          
        }
      }
    }
    if(!is.null(dist_pock)) {
      dist_ligs[flag] = mean(dist_pock)
      flag = flag+1      
    }
  }
}
length(which(dist_ligs > 0))
# dist for different ligands

dist_ligs_random = rep(0,5000)#length(unique(names_ligand)))
for (i in 1:5000){#length(unique(names_ligand))){
  print(i)
  lig = sample(names_ligand_unique, 2)
  f_1=grep(paste0(paste0("_",lig[1]),"_"),rownames(dt))
  f_2=grep(paste0(paste0("_",lig[2]),"_"),rownames(dt))
  if(length(f_1) > 1) {
    f_1 = sample(f_1,1)
  }
  if(length(f_2) > 1) {
    f_2 = sample(f_2,1)
  }
  dist_ligs_random[i] = dist_fuzcav_ph(as.integer(dt[f_1,]),
                                       as.integer(dt[f_2,]))
}
hist(dist_ligs_random)
d_random = density(dist_ligs_random)
plot(d_random)
#
load(file = "../results/pharmacophores_results/dist_ligs.Rdata")
load(file = "../results/pharmacophores_results/dist_ligs_random.Rdata")
sm.density.compare(c(dist_ligs_random,dist_ligs[which(dist_ligs>0)]), 
                   c(rep(1,length(dist_ligs_random)),rep(2,length(dist_ligs[which(dist_ligs>0)]))), 
                   model = "none", xlim=c(0,1))

sm.density.compare(c(dist_ligs_random,
                     dist_ligs[which(dist_ligs>0)]),
                   c(rep(1,length(dist_ligs_random)),
                     rep(2,length(dist_ligs[which(dist_ligs>0)]))
                   ),
                   model = "none", xlim=c(0,1)
                   , xlab = "pharmacophore distance bewteen pockets"
                   , main = "Density plot of the pharmacophore similarity between pockets from different ligands linking the same ligand")


abline(v=mean(dist_ligs[1:flag]), col = "blue")
abline(v=mean(dist_ligs_random), col = "red")
text(0.6,8,"mean = 0.43",col="blue")
text(0.3,8,"mean = 0.15",col="red")
mean(dist_ligs_random)
mean(dist_ligs[1:flag])

#check results
v_ph_1 = as.integer(dt[index[1],])
v_ph_2 = as.integer(dt[index[2],])
t = table(v_ph_1, v_ph_2)

common_non_null = sum(t[2:nrow(t),2:ncol(t)])
non_null_counts_Fa = sum(t[2:nrow(t),])
non_null_counts_Fb = sum(t[,2:ncol(t)])
val_similarity = common_non_null/min(non_null_counts_Fa,non_null_counts_Fb)


#test call c function
helloA <- function() {
  system(paste(getwd(),"helloA",sep="/"))
}
helloA()

dyn.load("helloB.so")
helloB <- function() {
  result <- .C("helloB",
               greeting="")
  return(result$greeting)
}
greeting <- helloA()
class(greeting)
greeting <- helloB()
class(greeting)
greeting

dyn.load("helloC.so")
helloC <- function(greeting) {
  if (!is.character(greeting)) {
    stop("Argument 'greeting' must be of type 'character'.")
  }
  result <- .C("helloC",
               greeting=greeting,
               count=as.integer(1))
  return(result$count)
}
helloC("Bonjour tout le monde!")

library("Rcpp")
sourceCpp("C_code/dist_fuzcav.cpp")



cppFunction("bool isOddCpp(int num = 10) {bool result = (num % 2 == 1);return result;}")
cppFunction("
double min_cpp (double a , double b){
  if (a<b) return a;
  return b;
}")


dist_fuzcav_ph(as.integer(dt_pharmacophores[1,]),as.integer(dt_pharmacophores[2,]))


isOddCpp(42L)
#### number pockets in double ####
dist_ligs = rep(0,length(which(table(names_ligand) > 1)))#length(unique(names_ligand)))
names_ligand_unique = unique(names_ligand)

n_pockets_double = NULL
names_pockets_double = NULL
for (name_lig in names(which(table(names_ligand) > 10))){
  print(name_lig)
  index = grep(paste0(paste0("_",name_lig),"_"),rownames(dt))
  d_t = dist(dt[index,])
  n_pockets_double = c(n_pockets_double, length(which(d_t == 0)))
  
  d_m = as.matrix(d_t)
  for(i in 2:nrow(d_m)) {
    if(min(d_m[i,1:(i-1)]) == 0) {
      #print(i)
      #print(rownames(d_m)[i])
      names_pockets_double = c(names_pockets_double,rownames(d_m)[i])
    }
  }
  if(sum(n_pockets_double) != length(names_pockets_double))
    print(name_lig)
}
boxplot(n_pockets_double)
sum(n_pockets_double)
length(names_pockets_double)

save(names_pockets_double, file = "../results/names_pockets_double", version = 2)
load("../results/names_pockets_double")

attributes(d_t)$Labels
#
#### ENZYM FAMILY in DT ####
#
names_prot = sapply(strsplit(rownames(dt), "_"), "[", 1)
#
pdb_names_hydrolases = scan("../data/enzymes_pdb_names/hydrolases.txt", character(), quote = "")
pdb_names_hydrolases = sub(",","",pdb_names_hydrolases)
length(intersect(pdb_names_hydrolases, names_prot))
n_poches_hydrolases = sum(table(names_prot)[intersect(pdb_names_hydrolases, names_prot)])
#
pdb_names_transferases = scan("../data/enzymes_pdb_names/tansferases.txt", character(), quote = "")
pdb_names_transferases = sub(",","",pdb_names_transferases)
length(intersect(pdb_names_transferases, names_prot))
n_poches_transferases = sum(table(names_prot)[intersect(pdb_names_transferases, names_prot)])
#
pdb_names_oxidoreductases = scan("../data/enzymes_pdb_names/oxidoreductases.txt", character(), quote = "")
pdb_names_oxidoreductases = sub(",","",pdb_names_oxidoreductases)
length(intersect(pdb_names_oxidoreductases, names_prot))
n_poches_oxidoreductases = sum(table(names_prot)[intersect(pdb_names_oxidoreductases, names_prot)])
#
pdb_names_lyases = scan("../data/enzymes_pdb_names/lyases.txt", character(), quote = "")
pdb_names_lyases = sub(",","",pdb_names_lyases)
length(intersect(pdb_names_lyases, names_prot))
n_poches_lyases = sum(table(names_prot)[intersect(pdb_names_lyases, names_prot)])
#
pdb_names_isomerase = scan("../data/enzymes_pdb_names/isomerases.txt", character(), quote = "")
pdb_names_isomerase = sub(",","",pdb_names_isomerase)
length(intersect(pdb_names_isomerase, names_prot))
n_poches_isomerase = sum(table(names_prot)[intersect(pdb_names_isomerase, names_prot)])
#
pdb_names_ligases = scan("../data/enzymes_pdb_names/ligases.txt", character(), quote = "")
pdb_names_ligases = sub(",","",pdb_names_ligases)
length(intersect(pdb_names_ligases, names_prot))
n_poches_ligases = sum(table(names_prot)[intersect(pdb_names_ligases, names_prot)])
#
pdb_names_translocases = scan("../data/enzymes_pdb_names/translocases.txt", character(), quote = "")
pdb_names_translocases = sub(",","",pdb_names_translocases)
length(intersect(pdb_names_translocases, names_prot))
n_poches_translocases = sum(table(names_prot)[intersect(pdb_names_translocases, names_prot)])
#
nrow(dt) - sum(n_poches_hydrolases,
    n_poches_transferases,
    n_poches_oxidoreductases,
    n_poches_lyases,
    n_poches_isomerase,
    n_poches_ligases,
    n_poches_translocases)
### drug bank ligands in dt ###
names_ligand = sapply(strsplit(rownames(dt), "_"), "[", 2)
unique_names_ligand = unique(names_ligand)
pdb_names_drugbank = read.csv("../data/drugbank/drugbank_all_target_polypeptide_ids.csv/all.csv")
colnames(pdb_names_drugbank)
class(pdb_names_drugbank$PDB.ID)

pdb_names_drugbank_PDB_ID = as.list(as.character(pdb_names_drugbank$PDB.ID))
pdb_names_drugbank_PDB_ID = strsplit(pdb_names_drugbank$PDB.ID[4], ";")
pdb_names_drugbank_PDB_ID[1]
unlist(pdb_names_drugbank$PDB.ID[2])
summary(pdb_names_drugbank)

n_drugbank_prot = NULL
nbr_drugbank_prot = 0
name_drugbank = NULL
for (i in 1:length(pdb_names_drugbank_PDB_ID)) {
  #print(any(intersect(toupper(names_prot), strsplit(pdb_names_drugbank_PDB_ID[[i]],"; "))))
  name_drugbank = c(name_drugbank,strsplit(pdb_names_drugbank_PDB_ID[[i]],"; ")[[1]])
  if(any(is.element(names_prot, strsplit(pdb_names_drugbank_PDB_ID[[i]],"; "))) ) {
    n_drugbank_prot = c(n_drugbank_prot,intersect(names_prot, strsplit(pdb_names_drugbank_PDB_ID[[i]],"; ")[[1]]))
    nbr_drugbank_prot = nbr_drugbank_prot + length(strsplit(pdb_names_drugbank_PDB_ID[[i]],"; ")[[1]])
  }
}

length(n_drugbank_prot)
nbr_drugbank_prot
length(unique(name_drugbank))

which(is.element(names_prot, name_drugbank) == TRUE)
names_lig_drgbk = names_ligand[which(is.element(names_prot, name_drugbank) == TRUE)]

length(unique(names_lig_drgbk))
which(table(names_lig_drgbk) > 20)

length(names_prot)
length(unique(names_ligand))
setdiff(unique(names_ligand), unique(names_lig_drgbk))
## on all the data
# Load the package required to read XML files.
library("XML")
# Also load the other required package.
library("methods")

# Give the input file name to the function.
full_drugbank <- xmlParse(file = "../data/drugbank/drugbank_all_full_database.xml/full database.xml")
test_db <- xmlParse(file = "../data/drugbank/drugbank.xsd") 
# Print the result.
print(test_db)
#
xmlRoot(test_db)

#n_drugbank = rep(0, length(pdb_names_drugbank$Ligand.ID))
n_drugbank = NULL
n_drugbank_pockets = 0
for(i in 1:length(pdb_names_drugbank$Ligand.ID)) {
  if(any(is.element(unique_names_ligand,  pdb_names_drugbank$Ligand.ID[i]))) {
    #n_drugbank[i] = 1
    n_drugbank = c(n_drugbank, intersect(unique_names_ligand,  pdb_names_drugbank$Ligand.ID[i]))
    n_drugbank_pockets = n_drugbank_pockets + length(which(names_ligand == intersect(unique_names_ligand,  pdb_names_drugbank$Ligand.ID[i])))
  }
}
#pdb_names_drugbank$Ligand.ID[1]
#sum(n_drugbank)
length(n_drugbank)

#### CLASSIC K MEANS ####
dt = dt_12descriptors[,]#names_pharmacophores_descriptors,]#c(-8,-2)]
dt = delete_clean_data(dt)
#scale data
dt = scale(dt)
#
dt.kmean = kmeans(dt, 10, nstart = 1, algorithm="MacQueen", iter.max = 200)
save(dt.kmean, file = "../results/dt.kmean_dt12clean_nseeds100.Rdata", version = 2)
load("../results/dt.kmean_dt12clean_nseeds100.Rdata")
#plot dist / taille
plot(density(sqrt(dt.kmean$withinss/dt.kmean$size)), main = paste("distance moyenne  n_seeds = ", seeds))
plot(dt.kmean$size, sqrt(dt.kmean$withinss/dt.kmean$size),
     main = paste("distance moyenne  n_seeds = ", seeds))

#
nstart = 100
mean_average_within_clusters_dist = NULL
for (seeds in 30:40) {
  print(seeds)
  dt.kmean = kmeans(dt, seeds, nstart = nstart, algorithm="MacQueen", iter.max = 300)
  average_within_clusters_dist = sqrt(dt.kmean$withinss/dt.kmean$size)
  mean_average_within_clusters_dist = c(mean_average_within_clusters_dist,
                                        mean(average_within_clusters_dist))
  
}
png("../results/K_Means_comparaison/distance_moyenne_cluster_30-40.png")
plot(30:40, mean_average_within_clusters_dist,
     xlab = "nombre de graines",
     ylab = "moyenne de distance moyenne entre les poches protéiques et le centroide de leur groupe")
dev.off()
#
library(colorspace)
c_colors = diverging_hsv(11)
nstart = 1
for (seeds in c(800, 1000)) { #800,1000
  print(seeds)
  load(paste0("../results/K_Means_comparaison/dt.kmean/dt.kmean_nseeds_", seeds))
  #dt.kmean = kmeans(dt, seeds, nstart = nstart, algorithm="MacQueen", iter.max = 300)
  #save(dt.kmean, file = paste0("../results/K_Means_comparaison/dt.kmean/dt.kmean_nseeds_", seeds))
  png(paste0(paste0("../results/K_Means_comparaison/distance_size_LP_nseeds_", seeds),
             ".png"))
  
  par(mar=c(5, 4, 4, 8), xpd=TRUE)
  plot(dt.kmean$size, sqrt(dt.kmean$withinss/dt.kmean$size), xlim = c(0,300), ylim = c(0,3.5),
       main = paste("distance moyenne  n_seeds = ", seeds))
  
  grp_sup_2 = which(sqrt(dt.kmean$withinss/dt.kmean$size) >=2)
  
  nbr_unique_lig = NULL
  nbr_unique_prot = NULL
  for (n_grp in 1:seeds) {
    names_grp = names(which(dt.kmean$cluster == n_grp))
    nbr_unique_lig = c(nbr_unique_lig,
                       length(unique(sapply(strsplit(names_grp, "_"), "[", 2))))
    nbr_unique_prot = c(nbr_unique_prot,
                        length(unique(sapply(strsplit(names_grp, "_"), "[", 1))))
  } 
  nbr_unique_prot_pch = nbr_unique_prot
  nbr_unique_prot_pch[which(nbr_unique_prot == 1)] = 13
  nbr_unique_prot_pch[which(nbr_unique_prot > 1)] = 10
  nbr_unique_prot_pch[which(nbr_unique_prot > 10)] = 1
  #points(mean(dt.kmean$size), mean(sqrt(dt.kmean$withinss/dt.kmean$size)), col = "red")
  
  points(dt.kmean$size[which(nbr_unique_lig == 1)], 
         sqrt(dt.kmean$withinss/dt.kmean$size)[which(nbr_unique_lig == 1)], 
         col = "red", pch = nbr_unique_prot_pch[which(nbr_unique_lig == 1)])
  points(dt.kmean$size[which(nbr_unique_lig > 1)], 
         sqrt(dt.kmean$withinss/dt.kmean$size)[which(nbr_unique_lig > 1)], 
         col = "orange", pch = nbr_unique_prot_pch[which(nbr_unique_lig > 1)])
  points(dt.kmean$size[which(nbr_unique_lig > 10)], 
         sqrt(dt.kmean$withinss/dt.kmean$size)[which(nbr_unique_lig > 10)], 
         col = "green", pch = nbr_unique_prot_pch[which(nbr_unique_lig > 10)])

  #text(dt.kmean$size, sqrt(dt.kmean$withinss/dt.kmean$size), 
  #     labels=paste(paste0("L:",nbr_unique_lig),
  #                  paste0("P:",nbr_unique_prot), sep = "|"),
  #     cex= 0.6, pos=3, col = "blue")
  
  #text(dt.kmean$size[grp_sup_2], sqrt(dt.kmean$withinss/dt.kmean$size)[grp_sup_2], 
  #     labels=paste(dt.kmean$size[grp_sup_2],
  #                  round(sqrt(dt.kmean$withinss/dt.kmean$size)[grp_sup_2], 2), sep = ","),
  #     cex= 0.6, pos=3, col = "blue")
  legend("topright", inset=c(-0.35,0), xpd=TRUE, mar(c(7,7,7,7)), cex = 1, bty = "n",
         legend=c(paste0("mean_dist:"),
                  paste0(round(mean(sqrt(dt.kmean$withinss/dt.kmean$size)),2),
                         paste0("±", 
                         paste0(round(sd(sqrt(dt.kmean$withinss/dt.kmean$size)),2)),"sd")),
                  paste0("mean_size:",round(mean(mean(dt.kmean$size)),2)),
                  paste0("n_dist>=2:",length(which(sqrt(dt.kmean$withinss/dt.kmean$size) >=2))),
                  paste0("n_pock>=2:",sum(dt.kmean$size[which(sqrt(dt.kmean$withinss/dt.kmean$size) >=2)])),
                  paste0("grp_size=1:",length(which(dt.kmean$size == 1))),
                  paste0("grp_sizeL=1:",length(which(nbr_unique_lig == 1))),
                  paste0("grp_sizeP=1:",length(which(nbr_unique_prot == 1))) ))
  legend( x="bottomright", inset=c(-0.375,0), xpd=TRUE, mar(c(7,7,7,7)), cex = 1, bty = "n",
          legend=c("grp_sizeL=[1]", 
                   "grp_sizeL=]1,10]",
                   "grp_sizeL=]10,[",
                   "grp_sizeP=[1]",
                   "grp_sizeP=]1,10]",
                   "grp_sizeP=]10,["),
          col=c("red","orange","green", "black", "black", "black"),
          lwd=1, lty=c(0,0), pch=c(1,1,1,13,10,1), x.intersp = 0)
  dev.off()
}
#
which(dt.kmean$size == 1)

hist(average_within_clusters_dist, main = "average_within_clusters_dist CLASSIC")
plot(density(average_within_clusters_dist), main = "average_within_clusters_dist CLASSIC n_seeds = 10")
abline(v=mean(average_within_clusters_dist), col = "blue")
text(2,1.25,"mean = 2.42",col="blue")

#mean between pockets
mean_node_intra_dist = NULL
sd_node_intra_dist = NULL
size_cluster = NULL
for (i in 1:100) {
  print(i)
  mean_node_intra_dist = c(mean_node_intra_dist,mean(dist(dt[which(dt.kmean$cluster == i),])))
  sd_node_intra_dist = c(sd_node_intra_dist, sd(dist(dt[which(dt.kmean$cluster == i),])))
  size_cluster = c(size_cluster, length(which(dt.kmean$cluster == i)))
}
intersect(dt.kmean$size, size_cluster)

plot(density(mean_node_intra_dist))

## analyses kmean ##
seeds = 800
load(paste0("../results/K_Means_comparaison/dt.kmean/dt.kmean_nseeds_", seeds))
grp_sup_2 = which(sqrt(dt.kmean$withinss/dt.kmean$size) < 0.2)
grp_sup_2 = which(dt.kmean$size == 1)
names_cluster_kmean = names(which(dt.kmean$cluster == 10))
lig_cluster_kmean = sapply(strsplit(names_cluster_kmean, "_"), "[", 2)
table(lig_cluster_kmean)
length(unique(lig_cluster_kmean))


table(lig_cluster_kmean)[names(which(table_names_ligand[unique(lig_cluster_kmean)] > 1))]

table_names_ligand[unique(lig_cluster_kmean)]

names(which(table_names_ligand > 1))
table_names_ligand["PE8"]
#plot size sur nbr ligands
png(paste0(paste0("../results/K_Means_comparaison/size_L_nseeds_", seeds),
".png"))
plot(nbr_unique_lig, dt.kmean$size, xlim = c(0,300), ylim = c(0,300) , 
     xlab = "Nombre ligands dans groupe", ylab = "Taille du groupe",
     main = paste("taille clusters et nbr ligands associés  n_seeds = ", seeds))
dev.off()
png(paste0(paste0("../results/K_Means_comparaison/size_P_nseeds_", seeds),
           ".png"))
plot(nbr_unique_prot, dt.kmean$size, xlim = c(0,300), ylim = c(0,300) , 
     xlab = "Nombre proteines dans groupe", ylab = "Taille du groupe",
     main = paste("taille clusters et nbr proteines associés  n_seeds = ", seeds))
dev.off()
#plot hist
png(paste0(paste0("../results/K_Means_comparaison/hist_L_nseeds_", seeds),
           "_inf10.png"))
hist(nbr_unique_lig[which(nbr_unique_lig < 10)],  
     xlab = "Nombre ligands dans groupe",
     main = paste("Hist nombre ligands associés n_seeds = ", seeds))
dev.off()
png(paste0(paste0("../results/K_Means_comparaison/hist_P_nseeds_", seeds),
           "_inf10.png"))
hist(nbr_unique_prot[which(nbr_unique_prot < 10)], xlim = c(0,300), 
     xlab = "Nombre proteines dans groupe",
     main = paste("Hist nombre proteines associés  n_seeds = ", seeds))
dev.off()
png(paste0(paste0("../results/K_Means_comparaison/hist_size_nseeds_", seeds),
           "_inf10.png"))
hist(dt.kmean$size[which(dt.kmean$size < 10)],
     xlab = "Nombre poches dans groupe",
     main = paste("Hist taille groupes n_seeds = ", seeds))
dev.off()
#HEM
hem_nbr = NULL
for (i in 1:length(dt.kmean$size)) {
  names_cluster_kmean = names(which(dt.kmean$cluster == i))
  lig_cluster_kmean = sapply(strsplit(names_cluster_kmean, "_"), "[", 2)
  hem_nbr = c(hem_nbr , length(which(lig_cluster_kmean == "HEM")))
  print(length(which(lig_cluster_kmean == "HEM")))
  #print(table(lig_cluster_kmean))
}
#MTA
mta_nbr = NULL
for (i in 1:length(dt.kmean$size)) {
  names_cluster_kmean = names(which(dt.kmean$cluster == i))
  lig_cluster_kmean = sapply(strsplit(names_cluster_kmean, "_"), "[", 2)
  mta_nbr = c(mta_nbr , length(which(lig_cluster_kmean == "MTA")))
  print(length(which(lig_cluster_kmean == "MTA")))
  #print(table(lig_cluster_kmean))
}
#REA
rea_nbr = NULL
for (i in 1:length(dt.kmean$size)) {
  names_cluster_kmean = names(which(dt.kmean$cluster == i))
  lig_cluster_kmean = sapply(strsplit(names_cluster_kmean, "_"), "[", 2)
  rea_nbr = c(rea_nbr , length(which(lig_cluster_kmean == "REA")))
  print(length(which(lig_cluster_kmean == "REA")))
  #print(table(lig_cluster_kmean))
}
#RBF
rbf_nbr = NULL
for (i in 1:length(dt.kmean$size)) {
  names_cluster_kmean = names(which(dt.kmean$cluster == i))
  lig_cluster_kmean = sapply(strsplit(names_cluster_kmean, "_"), "[", 2)
  rbf_nbr = c(rbf_nbr , length(which(lig_cluster_kmean == "RBF")))
  print(length(which(lig_cluster_kmean == "RBF")))
  #print(table(lig_cluster_kmean))
}

dist(dt[grep("HEM",rownames(dt)),])

#
#test :
test_clstr_1 = hclust(dist(dt[unlist(alltree$`1`$pockets_names),]), method = "ward.D2")
####hclust centroids####
dt_centers.hclust = hclust(dist(dt.kmean$centers), method = "ward.D2")
#dendrogram
dt_centers.dend = as.dendrogram(dt_centers.hclust)
#
library(dendextend)
library(colorspace)

##select lig##
lig_nbr = hem_nbr
#change label name
lab_cluster_dend = as.integer(labels(dt_centers.dend))
labels(dt_centers.dend) <- paste0(lab_cluster_dend, 
                                  paste0(";Nbr_HEM:",lig_nbr[lab_cluster_dend]))

#change label color
color.gradient <- function(x, colors=c("gold","darkgreen"), colsteps=100) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}
x <- c((1:100)^2, (100:1)^2)
plot(lig_nbr,col=color.gradient(lig_nbr), pch=19,cex=2)
c_colors = color.gradient(lig_nbr)

color.df<-data.frame(COLOR_VALUE=lig_nbr, color.name=color.gradient(lig_nbr))
color.df$color.name = as.character(color.df$color.name)
color.df[which(color.df$COLOR_VALUE == 0), "color.name"] = "#FFFFFF"
labels_colors(dt_centers.dend) <- as.character(color.df[as.character(order.dendrogram(dt_centers.dend)), "color.name"])

# Open a PDF for plotting; units are inches by default
pdf("../results/K_Means_comparaison/hclust_centroids_hem_ph.pdf", width=40, height=15)

# Do some plotting
par(cex=0.3, mar=c(5, 8, 6, 1))
plot(dt_centers.dend, type = "rectangle", ylab = "Height", cex.lab = 0.7, cex.main = 8,
     main = paste0("Nombre de ligands HEM dans les groupes. Nombre total:",sum(lig_nbr)))

legend("topright", title = "Occurence du ligands dans les groupes",
       legend = c(paste0("max:", round(max(lig_nbr),2)),
                  rep("",3),
                  rep("",4),
                  paste0("min:", round(min(lig_nbr),2))), pt.cex = 15, cex = 8, bty = "n")
names(lig_nbr) = 1:length(lig_nbr)
legend("topleft", title = "10 groupes les plus représentatifs:",
       legend = paste("Groupe", 
                      paste(names(sort(lig_nbr, decreasing = T)[1:10]),
                            paste("| Nbr", round(sort(lig_nbr,decreasing = T)[1:10],2)))),
       pt.cex = 15, cex = 8, bty = "n")
l_0 = length(which(lig_nbr == 0))
gradientLegend(valRange=c(-14,14), color=c(rep("#FFFFFF",10),
                                           color.gradient(sort(lig_nbr))[-c(1:l_0)]), 
               pos=c(680,40,730,63), coords=TRUE, 
               border.col=alpha('gray'), side=4)#pos.num = 3, length=.2, depth=.02)
# Close the PDF file's associated graphics device (necessary to finalize the output)
dev.off()

### K MEANS with specific dist for pharmacophores###
library(flexclust)
flexclust::distEuclidean
dist_Euc_test = function (x, centers) 
{
  if (ncol(x) != ncol(centers)) 
    stop(sQuote("x"), " and ", sQuote("centers"), " must have the same number of columns")
  z <- matrix(0, nrow = nrow(x), ncol = nrow(centers))
  for (k in 1:nrow(centers)) {
    z[, k] <- sqrt(colSums((t(x) - centers[k, ])^2))
  }
  z
}

dist_fuzcav = function (x, centers) {
  z <- matrix(0, nrow(x), ncol = nrow(centers))
  for (k in 1:nrow(centers)) {
    d_x = apply(x, 1, function(fingerprint) {1-dist_fuzcav_ph(as.integer(fingerprint), as.integer(centers[k,]))})
    #for (i in 1:nrow(x)) {
    #  d_x = c(d_x, 1 - dist_fuzcav_ph(as.integer(x[i,]), as.integer(centers[k,])))
    #}
    z[,k] <- d_x
  }
  z
}
#

start_time <- Sys.time()
res = kcca(dt_pharmacophores[1:20,], 3, family=kccaFamily(dist=dist_fuzcav))
end_time <- Sys.time()
end_time - start_time
res@cluster


res_test = kcca(dt[1:3,], 2, family=kccaFamily(dist=distEuclidean))
res_test@centers

dist_fuzcav_ph(as.integer(dt_pharmacophores["1A00_HEM_A_2",]), as.integer(dt_pharmacophores["1A01_HEM_A_2",]))

### H CLUSTwith specific dist for pharmacophores###
#dt_pharmacophores = read.table("../data/FPCount_save_all_inter_dt12.txt", sep = ";", row.names = 1, nrows = 20)

hclust_ph_dist = matrix(0, nrow(dt_pharmacophores[1:20,]), ncol = nrow(dt_pharmacophores[1:20,]))
rownames(hclust_ph_dist) = rownames(dt_pharmacophores[1:20,])
colnames(hclust_ph_dist) = rownames(dt_pharmacophores[1:20,])

for (i in 1:nrow(dt_pharmacophores[1:20,])) {
  for (j in i:nrow(dt_pharmacophores[1:20,])) {
    hclust_ph_dist[i,j] = 1 - dist_fuzcav_ph(as.integer(dt_pharmacophores[i,]), as.integer(dt_pharmacophores[j,]))
  }
}

hclust_ph_dist = as.dist(t(hclust_ph_dist))

ph.hclust = hclust(hclust_ph_dist, method = "ward.D2")
#plot(ph.hclust, hang = -1, cex = 0.6)

png(filename="../results/kmeans_hclust_test_k3_20pocks.png")
plot(ph.hclust, labels = res@cluster, hang = -1, cex = 0.6)
dev.off()
#### TEST COLLAPSIBLE TREE ####
library(collapsibleTree) 
# input data must be a nested data frame:
head(warpbreaks)
dim(warpbreaks)
p <- collapsibleTree( warpbreaks, c("wool", "tension", "breaks"))
p
#

c_1= NULL
c_2= NULL
c_3= NULL
for (i in strsplit(cluster_infos$pathString, "/")) {#nrow(cluster_infos)
  c_1 = c(c_1,i[2])
  c_2 = c(c_2,i[3])
  c_3 = c(c_3,i[3])
}
cluster_collaps = cbind(cluster_infos, c_1)
cluster_collaps = cbind(cluster_collaps, c_2)
cluster_collaps = cbind(cluster_collaps, c_3)
  
p <- collapsibleTree( cluster_collaps, c("c_1", "c_2", "c_3"), attribute ="size")
p

collapsibleTree(
  Geography,
  hierarchy = c("continent", "sub_region"),
  width = 800
)
#### pharmacophores kmedoids 800 seeds####

dt.kmedoids.ph.800 = read.table("../results/pharmacophores_results/clusters_nseeds800_Deuc.txt", sep = ",")

dt.kmedoids.ph.800[1:nrow(dt.kmedoids.ph.800),]
which(duplicated(dt.kmedoids.ph.800[1:nrow(dt.kmedoids.ph.800),])==TRUE)


dt.kmedoids.ph.800["4AP3_NAP_A_1",]
which(dt.kmedoids.ph.800[,1] == "4AP3_NAP_A_1")  
dt.kmedoids.ph.800[15120,] = NA

dt.kmedoids.ph.800 <- dt.kmedoids.ph.800[!duplicated(dt.kmedoids.ph.800[,1]),]
rownames(dt.kmedoids.ph.800) = dt.kmedoids.ph.800[1:nrow(dt.kmedoids.ph.800),1]
dt.kmedoids.ph.800[,2] = dt.kmedoids.ph.800[,2]+1
nrow(dt.kmedoids.ph.800)
#dt.kmedoids.ph.800[,1] = NULL
dt.kmedoids.ph.800 = as.data.frame(dt.kmedoids.ph.800)
hem_nbr = NULL
for (i in 1:800) {
  names_cluster_kmean = rownames(dt.kmedoids.ph.800[which(dt.kmedoids.ph.800[,2] == i),])
  lig_cluster_kmean = sapply(strsplit(names_cluster_kmean, "_"), "[", 2)
  hem_nbr = c(hem_nbr , length(which(lig_cluster_kmean == "HEM")))
  print(length(which(lig_cluster_kmean == "HEM")))
  #print(table(lig_cluster_kmean))
}
max(dt.kmedoids.ph.800[,2])
#### VALIDATION PROTOCOLE ####
seeds = 800
load(paste0("../results/K_Means_comparaison/dt.kmean/dt.kmean_nseeds_", seeds))
## POUR 800 GRAINES
names_prot = sapply(strsplit(rownames(dt), "_"), "[", 1)
names_ligand = sapply(strsplit(rownames(dt), "_"), "[", 2)
table_names_ligand = table(names_ligand)
#name ligands seen more than 1 time
names_ligand_unique_sup1 = names(which(table_names_ligand > 1))
#list of ligand names in clusters
lig_cluster_kmean = list()
for (i in 1:length(dt.kmean$size)) {
  names_cluster_kmean = names(which(dt.kmean$cluster == i))
  lig_cluster_kmean[[i]] = sapply(strsplit(names_cluster_kmean, "_"), "[", 2)
}
lig_cluster_kmean
#list of protein names in clusters
prot_cluster_kmean = list()
for (i in 1:length(dt.kmean$size)) {
  names_cluster_kmean = names(which(dt.kmean$cluster == i))
  prot_cluster_kmean[[i]] = sapply(strsplit(names_cluster_kmean, "_"), "[", 1)
}
prot_cluster_kmean
###Si ligand se retrouve dans cluster le plus représenté###
vec_number_list = list()
vec_group_list = list()
flag = 0
for (lig_name in names_ligand_unique_sup1) {
  flag = flag+1
  print(flag)
  vec_number = NULL
  for (group in 1:length(lig_cluster_kmean)) {
    lig_cluster_kmean_group = table(lig_cluster_kmean[[group]])
    if(is.na(lig_cluster_kmean_group[lig_name])) {
      vec_number = c(vec_number, 0)
    } else {
      vec_number = c(vec_number, lig_cluster_kmean_group[lig_name])
    }
  }
  vec_number_list[[lig_name]] = vec_number
  vec_group_list[[lig_name]] = which(vec_number == max(vec_number))
}
length(vec_number_list)
length(vec_group_list)

dt_validation = data.frame( row.names = names_ligand_unique_sup1,
                            number = I(vec_number_list),
                            group = I(vec_group_list))
save(dt_validation, file = "../results/prediction_validation/data/dt_validation.Rdata")
load("../results/prediction_validation/data/dt_validation.Rdata")
#count number pocket in cluster most representative / sum
res = sapply(vec_number_list,max)/sapply(vec_number_list,sum)
boxplot(res)
hist(res)
inf_10 = which(sapply(vec_number_list,sum) <10)
inf_100 = which(sapply(vec_number_list,sum) <100 & sapply(vec_number_list,sum) >10)
sup_10 = which(sapply(vec_number_list,sum) >10)
sup_100 = which(sapply(vec_number_list,sum) >100)
#
hist(res[inf_10])
length(res[inf_10])

hist(res[inf_100])
length(res[inf_100])

hist(res[sup_10])
length(res[sup_10])

hist(res[sup_100])
length(res[sup_100])
#
#Selon seuil de distance entre clusters
dist_seuil = 2
##chercher un ligand
dt_validation["0WP",]
###prospect benchamark###
###WITH PDB NAMES###
#dataset : barelier_structures
barelier_structures = read.csv("../data/benchmark/barelier_structures.tar/barelier_structures/barelier_structures/barelier_structures.csv", header = FALSE)
barelier_structures[,1] = toupper(barelier_structures[,1])
barelier_structures[,2] = toupper(barelier_structures[,2])
#dataset : identical_structures
identical_structures = read.csv("../data/benchmark/identical_structures.tar/identical_structures/identical_structures/identical_structures.csv", header = FALSE)
identical_structures[,1] = substr(toupper(identical_structures[,1]),1,4)
identical_structures[,2] = substr(toupper(identical_structures[,2]),1,4)
###
###
dt_benchmark = barelier_structures
dt_benchmark = identical_structures
###
y_true = rep(0, nrow(dt_benchmark))
y_true[which(dt_benchmark$V3 == "active")] = 1
y_true
#
y_predict = rep(0, nrow(dt_benchmark))
#
count = 0
n_line = 0
for (names_2prot in intersect(dt_in_1,dt_in_2)) {
  # if(is.element(dt_benchmark[names_2prot,1], names_prot)) {
  #   count =  count+1
  #   print(count)
  # }
  n_line = n_line+1
  print(n_line)
  prediction = 0
  for (i in 1:length(prot_cluster_kmean)) {
    if(is.element(dt_benchmark[names_2prot,1], prot_cluster_kmean[[i]]) & 
       is.element(dt_benchmark[names_2prot,2], prot_cluster_kmean[[i]])) {
      prediction = 1
      print(names_2prot)
      print("yes")
    }
  } 
  y_predict[names_2prot] = prediction
}
count
length(which(y_true == 1))
table(y_predict,y_true) 


#count number in data
length(which(is.element(dt_benchmark[,1], names_prot) == TRUE))
length(which(is.element(dt_benchmark[,2], names_prot) == TRUE))
dt_in_1 = which(is.element(dt_benchmark[,1], names_prot) == TRUE)
dt_in_2 = which(is.element(dt_benchmark[,2], names_prot) == TRUE)
length(intersect(dt_in_1,dt_in_2))
###WITH POCKET DESCRIPTORS###



