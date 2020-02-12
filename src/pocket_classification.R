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
#pharmacophores
row_names_pharmacophores = paste(gsub("PDB=| ","",dt_pharmacophores[,1]), gsub(".*/prox5_5/|_prox5_5_res_res-renum.pdb", "", rownames(dt_pharmacophores)), sep = "_")
rownames(dt_pharmacophores) = toupper(row_names_pharmacophores)
#delete coloumn of names
dt_pharmacophores[,1] = NULL

#descripors 72
rownames(dt_72descriptors) = toupper(gsub("_prox5_5.desR","",rownames(dt_72descriptors)))

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
  dt = na.omit(dt)
  minimum_size_pocket = 60
  maximum_size_pocket = 14
  inf_60 = which(dt[,"C_RESIDUES"] <= minimum_size_pocket)
  sup_14 = which(dt[,"C_RESIDUES"] >= maximum_size_pocket)
  dt = dt[intersect(inf_60,sup_14),]
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
#comparaison données
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
#dt.acp = PCA(dt, scale.unit = T,quali.sup=13) #normalisé automatiquement
dt.acp = PCA(dt, scale.unit = T)
dt.acp$eig
dt.acp$ind$contrib
dt.acp$var$contrib
plot(dt.acp)
#plot(dt.acp,habillage = 13, col.hab = c("green","blue"), label="none")

####K-MEANS####
##K-MEANS nK = 10%; nK = 10
dt = dt_12descriptors
dt = delete_clean_data(dt)
#
nbr_k = nrow(dt)*20/100
nbr_k = as.integer(nbr_k)
nbr_k = 20
start_time <- Sys.time()
dt.kmean = kmeans(scale(dt), nbr_k, nstart = 10)
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
n=100
for (i in 1:n) {
  k = kmeans(scale(dt),centers =  i, nstart = 10) #change nstart to 10
  R2 = c(R2,k$tot.withinss/k$totss)
  inertie.expl = c(inertie.expl,k$betweenss/k$totss)
  Iintra = c(Iintra.expl,k$betweenss)
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
plot(dt.acp,choix = "ind",col.quali = dt.kmean$cluster, title = "projection kmeans sur ACP", label="none")

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
#lets suppose we already cleaned the data: dt
nstart = 1
prct_seed = 10/100#10/100
min_size_cluster =100

dt.kmean.cluster = NULL#[]
dt.kmean.cluster.centers = NULL#[]
dt.kmean.cluster.Iintra = NULL#[]
dt.kmean.cluster.size = NULL#[]
threshold_recut_cluser = 0
valid_cluster = c(FALSE)




#first k mean
nbr_k = nrow(dt)*prct_seed # select number of seed : 10% of the size of the data
nbr_k = as.integer(nbr_k)
dt.kmean = kmeans(scale(dt), nbr_k, nstart = nstart)
#hclust on centroids
dt_centers.hclust = hclust(dist(scale(dt.kmean$centers)), method = "ward.D2")
plot(dt_centers.hclust,main = "", hang = -0.1)
nbr_k_optimal = 10
#second k mean with optimal number of seeds
dt.kmean = kmeans(scale(dt), nbr_k_optimal, nstart = nstart)
#hclust on optimal number of centroids
#dt_centers.hclust = hclust(dist(scale(dt.kmean$centers)), method = "ward.D2")
#dend1 <- as.dendrogram(dt_centers.hclust)
#alltree <- as.Node(dend1)
#population$height

#save centroids in data.tree
path_tree = "alltree"
pockets_cluster = NULL

for (i in 1:nbr_k_optimal){
  pockets_cluster[[i]] = names(which(dt.kmean$cluster == i))
}
cluster_dt = data.frame(centers = dt.kmean$centers,
                        withinss = dt.kmean$withinss,
                        size = dt.kmean$size,
                        pockets_names = I(pockets_cluster))

cluster_dt$pathString = paste(path_tree, 1, 1:nbr_k_optimal,sep = "/")

tmp = rbind(tmp, cluster_dt)
population <- as.Node(tmp)
print(population, "size", limit = 20)
#TODO : add list of names in infos of the tree
for (i in 1:nbr_k_optimal){
  pockets_cluster[[i]] = names(which(dt.kmean$cluster == i))
}

pockets_classification_tree = function(dt,
                                       nstart = 1,
                                       prct_seed = 1/100,
                                       path_tree = "alltree"){
  #first k mean
  nbr_k = nrow(dt)*prct_seed # select number of seed : 10% of the size of the data
  nbr_k = as.integer(nbr_k)
  dt.kmean = kmeans(scale(dt), nbr_k, nstart = nstart)
  #hclust on centroids
  print(nrow(dt.kmean$centers))
  dt_centers.hclust = hclust(dist(scale(dt.kmean$centers)), method = "ward.D2")
  #plot(dt_centers.hclust,main = "", hang = -0.1)
  ##
  #TODO:select K optimal
  ##
  nbr_k_optimal = 10
  #seceond kmean
  dt.kmean = kmeans(scale(dt), nbr_k_optimal, nstart = nstart)
  pockets_cluster = list()
  for (i in 1:nbr_k_optimal){
    pockets_cluster[[i]] = names(which(dt.kmean$cluster == i))
  }
  cluster_dt = data.frame(centers = dt.kmean$centers,
                          withinss = dt.kmean$withinss,
                          size = dt.kmean$size,
                          pockets_names = I(pockets_cluster))
  
  cluster_dt$pathString = paste(path_tree, 1:nbr_k_optimal,sep = "/")
  
  
  return(cluster_dt)
}
dt = dt_12descriptors[,]
dt = delete_clean_data(dt)
all_valide_clusters = F
path_tree = c("alltree")
list_path_tree = path_tree
list_n_cluster = NULL
list_pockets_cluster_names = NULL
cluster_infos = NULL
nbr_k= 1
pockets_cluster_names = list(row.names(dt))
iter = 1
valid_cluster = c(FALSE)
while(FALSE %in% valid_cluster) {
  for(i in 1:length(valid_cluster)){
    if (valid_cluster[i] == FALSE) {
      print("test")
      cluster_dt = pockets_classification_tree(
                                       dt = dt[unlist(pockets_cluster_names[[i]]),],
                                       path_tree = list_path_tree[i])
      list_n_cluster = c(list_n_cluster, nrow(cluster_dt))
      pockets_cluster_names = c(pockets_cluster_names, cluster_dt[,"pockets_names"])
      cluster_infos = rbind(cluster_infos, cluster_dt)
    }
  }
  
  valid_cluster = NULL
  list_path_tree = NULL
  for (i in list_n_cluster) {
    print(i)
    for (j in 1:i){
      valid_cluster = c(valid_cluster, FALSE)
      list_path_tree = c(list_path_tree, paste(path_tree , j, sep ="/"))
    }
  }
  list_n_cluster = NULL
  if(iter == 3){
    valid_cluster = c(TRUE)
  }
  iter = iter + 1
  if(iter == 2) {
    print("t3")
    print(valid_cluster)
  }
}
alltree <- as.Node(cluster_infos)
print(alltree, "size")

#tester :mettre dans infps le n iteration et le numero de cluster


pockets_classification_recursion = fucntion(dt, names_clust, n_k, i_k, path_tree) {
  cluster_dt = pockets_classification_tree(
                                          dt = dt[names_clust,],
                                          path_tree = paste0(path_tree,i))
  if(length(cluster_dt) < 100 && n_k == i_k){
    print(path_tree)
    return(TRUE)
  }
  if(length(cluster_dt) < 100 && n_k > i_k){
    return(pockets_classification_recursion(dt,
                                            n_k = n_k,
                                            i_k = i_k+1,
                                            path_tree = paste0(path_tree , i, sep ="/")))
  }
  else {
    pockets_classification_recursion
  }
}

#
library(pvclust)
library(MASS)
dt.pv <- pvclust(t(scale(dt.kmean$centers)))
plot(dt.pv)
#
dev.off()

alltree <- as.Node(dt.hclust)
dend1 <- as.dendrogram(dt.hclust)
alltree <- as.Node(dend1)
alltree$fieldsAll
alltree$totalCount
alltree$leafCount
alltree$height
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








