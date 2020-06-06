#### Analyse K means PYTHON ####
#data
dt = read.csv("../data/dt_72clean-50.csv",row.names = 1)
#dt[,c("drugg",hydrophobicity,polarity,HEIGHT_MIN_CYLINDER)]
##poches 50
#pocket_samePsameL_50
load("../results/kmeans_results_reglog/pocket_samePsameL_50.Rdata")
#pocket_diffPsameL_50
load("../results/kmeans_results_reglog/pocket_diffPsameL_50.Rdata")
#pocket_diffPdiffL_50
load("../results/kmeans_results_reglog/pocket_diffPdiffL_50.Rdata")
##FEATURES
features = c(descriptors_hydrophobicity,
             descriptors_aromatic,
             descriptors_polarity,
             descriptors_physicochemical,
             descriptors_geometrical)
features = c(features,"A","C","E","D","G","F","I","H","K","M","L","N","Q","P","S","R","T","W","V","Y")

##
union(union(pocket_samePsameL_50,pocket_diffPsameL_50),pocket_diffPdiffL_50)
write(union(union(pocket_samePsameL_50,pocket_diffPsameL_50),pocket_diffPdiffL_50), file = "../data/liste_pdb_overlpap.txt")
write(rownames(dt), file = "../data/liste_pdb_overlpap_all.txt")
#### ANALYSIS ####
#path = "../results/kmeans_results_ph_euc_50/"
path = "../data/kmeans_results_reglog/"
#
Nseeds = c(5,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000)
kmean_It = 35354.82457577707#scan(file =  paste0(path,"/kmeans_reglog_TSS.txt")
#
seeds = 5
#
kmean_results = read.table(file = paste0(paste0(paste0(path,"kmeans_reglog_seeds"),seeds),"_clusters.txt"), sep = ",", header = F, row.names = 1)
kmean_size = as.vector(table(kmean_results[,1]))
kmean_SSE =  read.table(file = paste0(paste0(paste0(path,"kmeans_reglog_SSE_seeds"),seeds),".txt"), sep = ",",header = T)
kmean_Ia = kmean_SSE[,"SSE"]/kmean_size
kmean_Ie = kmean_It - sum(kmean_SSE[,"SSE"])
R2 = sum(kmean_SSE[,"SSE"])/kmean_It
#
R2 = NULL
kmean_SSE_SSE = NULL
for (seeds in Nseeds) {
  print(seeds)
  kmean_results = read.table(file = paste0(paste0(paste0(path,"kmeans_reglog_seeds"),seeds),"_clusters.txt"), sep = ",", header = F,row.names = 1)
  kmean_size = as.vector(table(kmean_results[,1]))
  kmean_SSE =  read.table(file = paste0(paste0(paste0(path,"kmeans_reglog_SSE_seeds"),seeds),".txt"), sep = ",",header = T)
  kmean_SSE_SSE = c(kmean_SSE_SSE, mean(sqrt(kmean_SSE[,"SSE"]/kmean_size)))
  kmean_Ia = kmean_SSE[,"SSE"]/kmean_size
  #kmean_Ie = kmean_It - sum(kmean_SSE[,"SSE"])
  #R2 = c(R2,sum(kmean_SSE[,"SSE"])/kmean_It)
  
}
R2
plot(Nseeds,R2)
plot(kmean_size, sqrt(kmean_SSE[,"SSE"]/kmean_size),ylim = c(0,1))
plot(Nseeds, kmean_SSE_SSE)


### PLOT #be ware output file path
for (seeds in Nseeds) { #800,1000
  print(seeds)
  kmean_results = read.table(file = paste0(paste0(paste0(path,"kmeans_reglog_seeds"),seeds),"_clusters.txt"), sep = ",", header = F, row.names = 1)
  kmean_size = as.vector(table(kmean_results))
  kmean_SSE =  read.table(file = paste0(paste0(paste0(path,"kmeans_reglog_SSE_seeds"),seeds),".txt"), sep = ",",header = T)
  png(paste0(paste0("../results/kmeans_results_reglog/fig/distance_size_LP_nseeds_", seeds),".png"))
  
  par(mar=c(5, 4, 4, 8), xpd=TRUE)
  plot(kmean_size, sqrt(kmean_SSE[,"SSE"]/kmean_size),
       main = paste("distance moyenne  n_seeds = ", seeds))
  
  grp_sup_2 = which(sqrt(kmean_SSE[,"SSE"]/kmean_size) >=0.4)
  
  nbr_unique_lig = NULL
  nbr_unique_prot = NULL
  for (n_grp in 0:(seeds-1)) {
    names_grp = rownames(kmean_results)[which(kmean_results == n_grp)]
    nbr_unique_lig = c(nbr_unique_lig,
                       length(unique(sapply(strsplit(names_grp, "_"), "[", 2))))
    nbr_unique_prot = c(nbr_unique_prot,
                        length(unique(sapply(strsplit(names_grp, "_"), "[", 1))))
  } 
  nbr_unique_prot_pch = nbr_unique_prot
  nbr_unique_prot_pch[which(nbr_unique_prot == 1)] = 13
  nbr_unique_prot_pch[which(nbr_unique_prot > 1)] = 10
  nbr_unique_prot_pch[which(nbr_unique_prot > 10)] = 1
  
  points(kmean_size[which(nbr_unique_lig == 1)], 
         sqrt(kmean_SSE[,"SSE"]/kmean_size)[which(nbr_unique_lig == 1)], 
         col = "red", pch = nbr_unique_prot_pch[which(nbr_unique_lig == 1)])
  points(kmean_size[which(nbr_unique_lig > 1)], 
         sqrt(kmean_SSE[,"SSE"]/kmean_size)[which(nbr_unique_lig > 1)], 
         col = "orange", pch = nbr_unique_prot_pch[which(nbr_unique_lig > 1)])
  points(kmean_size[which(nbr_unique_lig > 10)], 
         sqrt(kmean_SSE[,"SSE"]/kmean_size)[which(nbr_unique_lig > 10)], 
         col = "green", pch = nbr_unique_prot_pch[which(nbr_unique_lig > 10)])
  legend("topright", inset=c(-0.35,0), xpd=TRUE, mar(c(7,7,7,7)), cex = 1, bty = "n",
         legend=c(paste0("mean_dist:"),
                  paste0(round(mean(sqrt(kmean_SSE[,"SSE"]/kmean_size)),2),
                         paste0("±", 
                                paste0(round(sd(sqrt(kmean_SSE[,"SSE"]/kmean_size)),2)),"sd")),
                  paste0("mean_size:",round(mean(mean(kmean_size)),2)),
                  paste0("n_dist>=2:",length(which(sqrt(kmean_SSE[,"SSE"]/kmean_size) >=2))),
                  paste0("n_pock>=2:",sum(kmean_size[which(sqrt(kmean_SSE[,"SSE"]/kmean_size) >=2)])),
                  paste0("grp_size=1:",length(which(kmean_size == 1))),
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
### 
### PERFORMANCES 50 POCHES ###
data_samePsameL_dt72_pharmacophores = read.csv("../data/data_structure_comparaison/data_samePsameL_dt72_pharmacophores.csv", colClasses = "character")
data_diffPsameL_dt72_pharmacophores = read.csv("../data/data_structure_comparaison/data_diffPsameL_dt72_pharmacophores.csv", colClasses = "character")
data_diffPdiffL_dt72_pharmacophores = read.csv("../data/data_structure_comparaison/data_diffPdiffL_dt72_pharmacophores.csv", colClasses = "character")
#MODEL
load("../results/kmeans_results_reglog/model.glm.step.Rdata")
#save(dt_predict.glm.step,file="../results/kmeans_results_reglog/model.glm.step.Rdata",version=2)
#
load("../results/kmeans_results_reglog/model.glm.step_overlap.Rdata")
#save(dt_predict.glm.step,file="../results/kmeans_results_reglog/model.glm.step_overlap.Rdata",version=2)

#
pocket_samePsameL_50
pocket_diffPsameL_50
pocket_diffPdiffL_50

##
set.seed(83)
pocket_samePsameL_50 = data_samePsameL_dt72_pharmacophores[sample(nrow(data_samePsameL_dt72_pharmacophores),
                                                                  size = 50),]
length(pocket_samePsameL_50)
pocket_diffPsameL_50 = data_diffPsameL_dt72_pharmacophores[sample(nrow(data_diffPsameL_dt72_pharmacophores),
                                                                  size = 50),]
length(pocket_diffPsameL_50)
pocket_diffPdiffL_50 = data_diffPdiffL_dt72_pharmacophores[sample(nrow(data_diffPdiffL_dt72_pharmacophores),
                                                                  size = 50),]
length(pocket_diffPdiffL_50)
##
#dt=as.data.frame(dt)
#colnames_dt = colnames(dt)
##
seeds = 100
seuil = 0.8
kmean_results = read.table(file = paste0(paste0("../results/kmeans_results_reglog_PproxMfpo/kmeans_reglog_seeds",seeds),"_clusters.txt"), sep = ",", header = F, row.names = 1)
kmean_centroids = read.table(file = paste0(paste0("../results/kmeans_results_reglog_PproxMfpo/kmeans_reglog_seeds",seeds),"_means.txt"), sep = ",", header = F)
colnames(kmean_centroids) = colnames_dt
#sameLsameP
y_predict_sameLsameP = NULL
nclust_predict_sameLsameP = NULL
for (i in 1:nrow(pocket_samePsameL_50)) {
  print(i)
  #pock2 found in cluster
  pock1_cluster = kmean_results[pocket_samePsameL_50[i,3],1]
  #pock2 predicted
  pock2_clusters = NULL
  for (j in 1:nrow(kmean_centroids)) {
    pock2_clusters = c(pock2_clusters, predict.glm(dt_predict.glm.step, newdata = sqrt((kmean_centroids[j,features]-dt[pocket_samePsameL_50[i,2],features])**2), type = "response"))
  }
  pock2_cluster = which(pock2_clusters > seuil|pock2_clusters == max(pock2_clusters))-1
  nclust_predict_sameLsameP = c(nclust_predict_sameLsameP,length(pock2_cluster))
  if(is.element(pock1_cluster,pock2_cluster)) {
    y_predict_sameLsameP = c(y_predict_sameLsameP,1)
  } else {
    y_predict_sameLsameP = c(y_predict_sameLsameP,0)
  }
}
#sameLdiffP
y_predict_sameLdiffP = NULL
nclust_predict_sameLdiffP = NULL
for (i in 1:nrow(pocket_diffPsameL_50)) {
  print(i)
  #pock2 found in cluster
  pock1_cluster = kmean_results[pocket_diffPsameL_50[i,3],1]
  #pock2 predicted
  pock2_clusters = NULL
  for (j in 1:nrow(kmean_centroids)) {
    pock2_clusters = c(pock2_clusters, predict.glm(dt_predict.glm.step, newdata = sqrt((kmean_centroids[j,features]-dt[pocket_diffPsameL_50[i,2],features])**2), type = "response"))
  }
  pock2_cluster = which(pock2_clusters > seuil|pock2_clusters == max(pock2_clusters))-1
  nclust_predict_sameLdiffP = c(nclust_predict_sameLdiffP,length(pock2_cluster))
  if(is.element(pock1_cluster,pock2_cluster)) {
    y_predict_sameLdiffP = c(y_predict_sameLdiffP,1)
  } else {
    y_predict_sameLdiffP = c(y_predict_sameLdiffP,0)
  }
}
#diffLdiffP
y_predict_diffLdiffP = NULL
nclust_predict_diffLdiffP = NULL
for (i in 1:nrow(pocket_diffPsameL_50)) {
  print(i)
  #pock2 found in cluster
  pock1_cluster = kmean_results[pocket_diffPdiffL_50[i,3],1]
  #pock2 predicted
  pock2_clusters = NULL
  for (j in 1:nrow(kmean_centroids)) {
    pock2_clusters = c(pock2_clusters, predict.glm(dt_predict.glm.step, newdata = sqrt((kmean_centroids[j,features]-dt[pocket_diffPdiffL_50[i,2],features])**2), type = "response"))
  }
  pock2_cluster = which(pock2_clusters > seuil|pock2_clusters == max(pock2_clusters))-1
  nclust_predict_diffLdiffP = c(nclust_predict_diffLdiffP,length(pock2_cluster))
  if(is.element(pock1_cluster,pock2_cluster)) {
    y_predict_diffLdiffP = c(y_predict_diffLdiffP,1)
  } else {
    y_predict_diffLdiffP = c(y_predict_diffLdiffP,0)
  }
}
##PERFORMANCES
y_true = c(rep(1, length(y_predict_sameLsameP)),
           rep(1, length(y_predict_sameLdiffP)),
           rep(0, length(y_predict_diffLdiffP)))
y_predict = c(y_predict_sameLsameP,
              y_predict_sameLdiffP,
              y_predict_diffLdiffP)
table_clusters <- table(factor(y_predict,levels = 0:1), factor(y_true,levels=0:1))
Se = table_clusters[2,2] / (table_clusters[2,2] + table_clusters[1,2])
Sp = table_clusters[1,1] / (table_clusters[1,1] + table_clusters[2,1])
Se
Sp
#library(mltools)
sLsP = 40
sLdP = 41
dLdP = 20
TP = sLsP+sLdP
TN = 50-dLdP
FP = dLdP
FN = 50-sLsP+50-sLdP
MCC = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
MCC
#
mcc(preds = y_predict, actuals = y_true)
#
length(which(y_predict_sameLsameP == 1))
length(which(y_predict_sameLdiffP == 1))
length(which(y_predict_diffLdiffP == 1))
mean(nclust_predict_sameLsameP)
mean(nclust_predict_sameLdiffP)
mean(nclust_predict_diffLdiffP)


### POCHES FPOCKETS 50###
dt_fpocket = read.table("../data/data_PDB_72desc_overlap_50.txt", header = T, sep = "", row.names = 1, fill = TRUE)
#
scaled_center_dt_t = read.table(file = "../results/scale/scaled_center_dt72clean.Rdata", col.names = F, row.names = 1)
scaled_scale_dt_t = read.table(file = "../results/scale/scaled_scale_dt72clean.Rdata", col.names = F, row.names = 1)
scaled_center_dt = scaled_center_dt_t[,1]
names(scaled_center_dt) = rownames(scaled_center_dt_t)
scaled_scale_dt = scaled_scale_dt_t[,1]
names(scaled_scale_dt) = rownames(scaled_scale_dt_t)
#
dt_fpocket = as.data.frame(scale(dt_fpocket, scaled_center_dt, scaled_scale_dt))
#
#dt = dt_fpocket
nrow(dt)
#
dt_fpocket[pocket_samePsameL_50,]
dt_fpocket[pocket_diffPsameL_50,]
dt_fpocket[pocket_diffPdiffL_50,]
#
##dist ## attention dt a prendre à partir d protocol_classification_feature
dist_lig_sameP_sameL = NULL
dist_lig_diffP_sameL = NULL
dist_lig_diffP_diffL = NULL
for (i in 1:nrow(pocket_samePsameL_50)) {
  print(i)
  dist_lig_sameP_sameL = c(dist_lig_sameP_sameL, predict.glm(dt_predict.glm.step, newdata = sqrt((dt_fpocket[pocket_samePsameL_50[i,2],]-dt[pocket_samePsameL_50[i,2],])**2), type = "response"))
  dist_lig_diffP_sameL = c(dist_lig_diffP_sameL, predict.glm(dt_predict.glm.step, newdata = sqrt((dt_fpocket[pocket_diffPsameL_50[i,2],]-dt[pocket_diffPsameL_50[i,2],])**2), type = "response"))
  dist_lig_diffP_diffL = c(dist_lig_diffP_diffL, predict.glm(dt_predict.glm.step, newdata = sqrt((dt_fpocket[pocket_diffPdiffL_50[i,2],]-dt[pocket_diffPdiffL_50[i,2],])**2), type = "response"))
}
#
y_true = c(rep(1,length(dist_lig_sameP_sameL)+
                 length(dist_lig_diffP_sameL)),
           rep(0,length(dist_lig_diffP_diffL)))
y_predict = c(dist_lig_sameP_sameL,
              dist_lig_diffP_sameL,
              dist_lig_diffP_diffL)
#
y_true = c(rep(1,length(dist_lig_diffP_sameL)),
           rep(0,length(dist_lig_diffP_diffL)))
y_predict = c(dist_lig_sameP_sameL,
              dist_lig_diffP_sameL)
#library(ROCR)
dt.pred = prediction(y_predict, y_true)
dt.perf = performance(dt.pred, "tpr", "fpr")
plot(dt.perf)
dt.auc = performance(dt.pred, "auc")
attr(dt.auc, "y.values")
#
y_predict[which(y_predict > 0.5)] = 1
y_predict[which(y_predict <= 0.5)] = 0
glm.table <- table(y_predict, y_true)
Se = glm.table[2,2] / (glm.table[2,2] + glm.table[1,2])
Sp = glm.table[1,1] / (glm.table[1,1] + glm.table[2,1])
Se
Sp

y_predict[which(y_predict >= 0.5 & y_true == 0)]
#
pocket_diffPdiffL_50[which(pocket_diffPdiffL_50[,2] == "5HHF_62F_A_1"),]
#
pocket_diffPdiffL_50[which(pocket_diffPdiffL_50[,2] == "2HA6_SCK_B_4"),]
#
pocket_diffPdiffL_50[which(pocket_diffPdiffL_50[,2] == "3IND_593_A_1"),]
#### Comparaison avec score overlap ###
SO_50 = read.table("../data/SO_desc_fpocket_50.txt",row.names = 1,sep=",")
names_SO_50 = rownames(SO_50)
SO_50 = as.vector(t(SO_50))
names(SO_50)=names_SO_50
#
#y_predict_fpocket = y_predict
#y_predict = abs(y_predict_fpocket-y_predict)
#
##poches 50
#pocket_samePsameL_50
load("../results/kmeans_results_reglog/pocket_samePsameL_50.Rdata")
#pocket_diffPsameL_50
load("../results/kmeans_results_reglog/pocket_diffPsameL_50.Rdata")
#pocket_diffPdiffL_50
load("../results/kmeans_results_reglog/pocket_diffPdiffL_50.Rdata")
##
plot(SO_50[c(pocket_samePsameL_50)],y_predict[c(pocket_samePsameL_50)],xlim=c(0,1),ylim=c(0,1))
plot(SO_50[c(pocket_diffPsameL_50)],y_predict[c(pocket_diffPsameL_50)],xlim=c(0,1),ylim=c(0,1))
plot(SO_50[c(pocket_diffPdiffL_50)],y_predict[c(pocket_diffPdiffL_50)],xlim=c(0,1),ylim=c(0,1))

cor(SO_50[c(pocket_samePsameL_50)],y_predict[c(pocket_samePsameL_50)])
cor(SO_50[c(pocket_diffPsameL_50)],y_predict[c(pocket_diffPsameL_50)])
cor(SO_50[c(pocket_diffPdiffL_50)],y_predict[c(pocket_diffPdiffL_50)])
#
t.test(SO_50[c(pocket_samePsameL_50)],y_predict[c(pocket_samePsameL_50)])
t.test(SO_50[c(pocket_diffPsameL_50)],y_predict[c(pocket_diffPsameL_50)])
t.test(SO_50[c(pocket_diffPdiffL_50)],y_predict[c(pocket_diffPdiffL_50)])
#
class(SO_50)
class(y_predict)
#
library(ggplot2)
library(ggExtra)
#data(mpg, package="ggplot2")

# mpg <- read.csv("http://goo.gl/uEeRGu")

# Scatterplot

theme_set(theme_bw())  # pre-set the bw theme.
#mpg_select <- mpg[mpg$hwy >= 35 & mpg$cty > 27, ]
sLsP_50 = cbind(SO_50[c(pocket_samePsameL_50)],y_predict[c(pocket_samePsameL_50)])
sLsP_50 = as.data.frame(sLsP_50)
colnames(sLsP_50) = c("SO","ypred")
#
sLdP_50 = cbind(SO_50[c(pocket_diffPsameL_50)],y_predict[c(pocket_diffPsameL_50)])
sLdP_50 = as.data.frame(sLdP_50)
colnames(sLdP_50) = c("SO","ypred")
#
dLdP_50 = cbind(SO_50[c(pocket_diffPdiffL_50)],y_predict[c(pocket_diffPdiffL_50)])
dLdP_50 = as.data.frame(dLdP_50)
colnames(dLdP_50) = c("SO","ypred")
###
sP_50 = cbind(SO_50,y_predict)
sP_50 = as.data.frame(sP_50)
colnames(sP_50) = c("SO","ypred")
cor.test(sP_50$SO,sP_50$ypred)
###
g <- ggplot(sP_50, aes(SO, ypred)) +
  geom_count() +
  geom_smooth(method="lm", se=F) + xlim(0,0.6)+ylim(0,1)


ggMarginal(g, type = "histogram", fill="transparent")

p1 <- ggplot(sP_50, aes(SO,ypred)) +
  geom_smooth(method="lm", se=F) + xlim(0,0.6)+ylim(0,1)
ggMarginal(p1 + geom_point(), type = "histogram", fill="transparent")

#ggMarginal(g, type = "boxplot", fill="transparent")

# ggMarginal(g, type = "density", fill="transparent")
####pharmacophores ####
#data
dt_ph_50 = read.table("../data/FPCount_50.txt", sep = ";", row.names = 1)
row_names_pharmacophores = gsub(".*/pocket_pdb/|_prox5_5.pdb.*", "", rownames(dt_ph_50))
rownames(dt_ph_50) = toupper(row_names_pharmacophores)
#
dt = dt_ph_50
#
library("Rcpp")
sourceCpp("C_code/dist_fuzcav.cpp")
#sourceCpp("C_code/dist_fuzcav_norm.cpp")
#dist_fuzcav_ph(as.integer(dt_pharmacophores[1,]),as.integer(dt_pharmacophores[2,]))
similarity_ph = function(vec1,vec2) {
  return(dist(rbind(vec1,vec2)))
  #return(dist_fuzcav_ph(as.integer(vec1),as.integer(vec2)))
  #return(dist_fuzcav_ph_norm(as.integer(vec1),as.integer(vec2)))
}
#
seeds = 5
seuil = 90
kmean_results = read.table(file = paste0(paste0("../results/kmeans_results_ph_euc_50/kmeans_reglog_seeds",seeds),"_clusters.txt"), sep = ",", header = F, row.names = 1)
kmean_centroids = read.table(file = paste0(paste0("../results/kmeans_results_ph_euc_50/kmeans_reglog_seeds",seeds),"_means.txt"), sep = ",", header = F)
kmean_centroids[,ncol(kmean_centroids)] = NULL
#
colnames(kmean_centroids) = 1:ncol(kmean_centroids)
colnames(dt) = 1:ncol(dt)
#
#sameLsameP
y_predict_sameLsameP = NULL
nclust_predict_sameLsameP = NULL
for (i in 1:nrow(pocket_samePsameL_50)) {
  print(i)
  #pock2 found in cluster
  pock1_cluster = kmean_results[pocket_samePsameL_50[i,3],1]
  #pock2 predicted
  pock2_clusters = NULL
  for (j in 1:nrow(kmean_centroids)) {
    pock2_clusters = c(pock2_clusters, similarity_ph(kmean_centroids[j,],dt[pocket_samePsameL_50[i,2],]))
  }
  pock2_cluster = which(pock2_clusters < seuil|pock2_clusters == min(pock2_clusters))-1
  nclust_predict_sameLsameP = c(nclust_predict_sameLsameP,length(pock2_cluster))
  if(is.element(pock1_cluster,pock2_cluster)) {
    y_predict_sameLsameP = c(y_predict_sameLsameP,1)
  } else {
    y_predict_sameLsameP = c(y_predict_sameLsameP,0)
  }
}
#sameLdiffP
y_predict_sameLdiffP = NULL
nclust_predict_sameLdiffP = NULL
for (i in 1:nrow(pocket_diffPsameL_50)) {
  print(i)
  #pock2 found in cluster
  pock1_cluster = kmean_results[pocket_diffPsameL_50[i,3],1]
  #pock2 predicted
  pock2_clusters = NULL
  for (j in 1:nrow(kmean_centroids)) {
    pock2_clusters = c(pock2_clusters,similarity_ph(kmean_centroids[j,],dt[pocket_diffPsameL_50[i,2],]))
  }
  pock2_cluster = which(pock2_clusters < seuil|pock2_clusters == min(pock2_clusters))-1
  nclust_predict_sameLdiffP = c(nclust_predict_sameLdiffP,length(pock2_cluster))
  if(is.element(pock1_cluster,pock2_cluster)) {
    y_predict_sameLdiffP = c(y_predict_sameLdiffP,1)
  } else {
    y_predict_sameLdiffP = c(y_predict_sameLdiffP,0)
  }
}
#diffLdiffP
y_predict_diffLdiffP = NULL
nclust_predict_diffLdiffP = NULL
for (i in 1:nrow(pocket_diffPsameL_50)) {
  print(i)
  #pock2 found in cluster
  pock1_cluster = kmean_results[pocket_diffPdiffL_50[i,3],1]
  #pock2 predicted
  pock2_clusters = NULL
  for (j in 1:nrow(kmean_centroids)) {
    pock2_clusters = c(pock2_clusters, similarity_ph(kmean_centroids[j,],dt[pocket_diffPdiffL_50[i,2],]))
  }
  pock2_cluster = which(pock2_clusters < seuil|pock2_clusters == min(pock2_clusters))-1
  nclust_predict_diffLdiffP = c(nclust_predict_diffLdiffP,length(pock2_cluster))
  if(is.element(pock1_cluster,pock2_cluster)) {
    y_predict_diffLdiffP = c(y_predict_diffLdiffP,1)
  } else {
    y_predict_diffLdiffP = c(y_predict_diffLdiffP,0)
  }
}
##PERFORMANCES
y_true = c(rep(1, length(y_predict_sameLsameP)),
           rep(1, length(y_predict_sameLdiffP)),
           rep(0, length(y_predict_diffLdiffP)))
y_predict = c(y_predict_sameLsameP,
              y_predict_sameLdiffP,
              y_predict_diffLdiffP)
table_clusters <- table(factor(y_predict,levels = 0:1), factor(y_true,levels=0:1))
Se = table_clusters[2,2] / (table_clusters[2,2] + table_clusters[1,2])
Sp = table_clusters[1,1] / (table_clusters[1,1] + table_clusters[2,1])
Se
Sp
#library(mltools)
sLsP = 24
sLdP = 20
dLdP = 13
TP = sLsP+sLdP
TN = 50-dLdP
FP = dLdP
FN = 50-sLsP+50-sLdP
MCC = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
MCC
#
mcc(preds = y_predict, actuals = y_true)
#
length(which(y_predict_sameLsameP == 1))
length(which(y_predict_sameLdiffP == 1))
length(which(y_predict_diffLdiffP == 1))
mean(nclust_predict_sameLsameP)
mean(nclust_predict_sameLdiffP)
mean(nclust_predict_diffLdiffP)

##dist ## attention dt a prendre à partir d protocol_classification_feature
dist_lig_sameP_sameL = NULL
dist_lig_diffP_sameL = NULL
dist_lig_diffP_diffL = NULL
for (i in 1:nrow(pocket_samePsameL_50)) {
  print(i)
  dist_lig_sameP_sameL = c(dist_lig_sameP_sameL, similarity_ph(dt_ph_50[pocket_samePsameL_50[i,2],],dt_pharmacophores[pocket_samePsameL_50[i,3],]))
  dist_lig_diffP_sameL = c(dist_lig_diffP_sameL, similarity_ph(dt_ph_50[pocket_diffPsameL_50[i,2],],dt_pharmacophores[pocket_diffPsameL_50[i,3],]))
  dist_lig_diffP_diffL = c(dist_lig_diffP_diffL, similarity_ph(dt_ph_50[pocket_diffPdiffL_50[i,2],],dt_pharmacophores[pocket_diffPdiffL_50[i,3],]))
}
#
y_true = c(rep(1,length(dist_lig_sameP_sameL)+
                 length(dist_lig_diffP_sameL)),
           rep(0,length(dist_lig_diffP_diffL)))
y_predict = c(dist_lig_sameP_sameL,
              dist_lig_diffP_sameL,
              dist_lig_diffP_diffL)
#
y_true = c(rep(1,length(dist_lig_sameP_sameL)),
           rep(0,length(dist_lig_diffP_diffL)))
y_predict = c(dist_lig_sameP_sameL,
              dist_lig_diffP_diffL)
#
y_predict[which(is.na(y_predict) == TRUE)] = 0
#library(ROCR)
dt.pred = prediction(y_predict, y_true)
dt.perf = performance(dt.pred, "tpr", "fpr")
#plot(dt.perf)
dt.auc = performance(dt.pred, "auc")
attr(dt.auc, "y.values")
#
m_y_predict = mean(y_predict)
#
y_predict[which(y_predict > m_y_predict)] = 0
y_predict[which(y_predict != 0)] = 1
glm.table <- table(y_predict, y_true)
Se = glm.table[2,2] / (glm.table[2,2] + glm.table[1,2])
Sp = glm.table[1,1] / (glm.table[1,1] + glm.table[2,1])
Se
Sp

### TREE BUILDING ###
library(data.tree)
##
seeds = 30
#
path = "../data/kmeans_results_reglog/"
#path = "../results/kmeans_results_reglog/"
#
colnames_dt = colnames(dt)
#
kmean_results = read.table(file = paste0(paste0(paste0(path,"kmeans_reglog_seeds"),seeds),"_clusters.txt"), sep = ",", header = F, row.names = 1)
kmean_centroids = read.table(file = paste0(paste0("../data/kmeans_results_reglog/kmeans_reglog_seeds",seeds),"_means.txt"), sep = ",", header = F)
colnames(kmean_centroids) = colnames_dt
kmean_size = as.vector(table(kmean_results[,1]))
kmean_SSE =  read.table(file = paste0(paste0(paste0(path,"kmeans_reglog_SSE_seeds"),seeds),".txt"), sep = ",",header = T)
#
pockets_cluster = list()
cluster_dend = list()
for (i in 1:seeds){
  print(i)
  pockets_cluster[[i]] = rownames(kmean_results)[which(kmean_results == i-1)]
  #mat_dist = matrix(0,nrow = length(pockets_cluster[[i]]),ncol=length(pockets_cluster[[i]]))
  #colnames(mat_dist) = pockets_cluster[[i]]
  #rownames(mat_dist) = pockets_cluster[[i]]
  #print(nrow(mat_dist))
  #v = v + nrow(mat_dist)*nrow(mat_dist)
  #for (j in 1:nrow(mat_dist)) {
  #  print(j)
  #  print("---")
  #  for (k in 1:j) {
  #    #print(k)
  #    v= v+1
  #    print(v)
  #    mat_dist[j,k] = 1-predict.glm(dt_predict.glm.step, newdata = as.data.frame(sqrt((dt[pockets_cluster[[i]][j],]-dt[pockets_cluster[[i]][k],])**2)), type = "response")
  #  }
  #}
  #
  mat_dist = read.csv(paste0(paste0("../results/kmeans_results_reglog/mat_dist_30/mat_dist_",i-1),".csv"),row.names = 1)
  cluster_dend[[i]] = as.dendrogram(hclust(as.dist(mat_dist),method = "ward.D2"))
}
# dt_pock = matrix(rep(dt[pockets_cluster[[i]][j],],each=n),nrow=n)
# dt_pock = dt[pockets_cluster[[i]][j],]
# 
# mat_dist = read.csv("../results/kmeans_results_reglog/mat_dist_30/mat_dist_0.csv",row.names = 1)
# nrow(mat_dist)
# ncol(mat_dist)
# rownames(mat_dist)
# colnames(mat_dist)
# test = hclust(as.dist(mat_dist),method = "ward.D2")
# plot(test)

cluster_dt = data.frame(centers = kmean_centroids,
                        withinss = kmean_SSE$SSE,
                        size = kmean_size,
                        #betweenss = rep(dt.kmean$betweenss, nbr_k_optimal),
                        #totss = rep(dt.kmean$totss, nbr_k_optimal),
                        pockets_names = I(pockets_cluster),
                        cluster_dend = I(cluster_dend)
                        
)
path_tree = c("alltree")
cluster_dt$pathString = paste(path_tree, 1:seeds,sep = "/")
cluster_infos = cluster_dt
#library(data.tree)
alltree <- as.Node(cluster_infos)



