#### Analyse K means PYTHON ####
#data
system.time(read.csv("../data/dt_72clean-50.csv",row.names = 1))
# code R
T1<-Sys.time()
dt = read.csv("../data/dt_72clean-50.csv",row.names = 1)
T2<-Sys.time()
T2-T1
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
path = "../results/kmeans_results_reglog_PproxMfpo/"
path = "../results/kmeans_results_reglog_article/"
#
Nseeds = c(5,10,20,30,40,50,60,70,80,90,100,200,300,400,500)
kmean_It = scan(file =  paste0(path,"/kmeans_reglog_TSS.txt"))#35354.82457577707#
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
I=NULL
kmean_SSE_SSE = NULL
for (seeds in Nseeds) {
  print(seeds)
  kmean_results = read.table(file = paste0(paste0(paste0(path,"kmeans_reglog_seeds"),seeds),"_clusters.txt"), sep = ",", header = F,row.names = 1)
  kmean_size = as.vector(table(kmean_results[,1]))
  kmean_SSE =  read.table(file = paste0(paste0(paste0(path,"kmeans_reglog_SSE_seeds"),seeds),".txt"), sep = ",",header = T)
  kmean_SSE_SSE = c(kmean_SSE_SSE, mean(sqrt(kmean_SSE[,"SSE"]/kmean_size)))
  kmean_Ia = kmean_SSE[,"SSE"]/kmean_size
  #kmean_Ie = kmean_It - sum(kmean_SSE[,"SSE"])
  I = c(I,sum(kmean_SSE[,"SSE"]))
  R2 = c(R2,sum(kmean_SSE[,"SSE"])/kmean_It)
  
}
R2
plot(Nseeds,R2,type = "l")
plot(kmean_size, sqrt(kmean_SSE[,"SSE"]/kmean_size),ylim = c(0,1))

mean(table(kmean_results))
sd(table(kmean_results))

length(table(kmean_results))

#kmean_SSE_SSE = I
plot(Nseeds, kmean_SSE_SSE, ylab = "Moyenne de similarité intra classe", xlab = "K graines", type = "l", col = "blue")
points(Nseeds, kmean_SSE_SSE, col = "red",type = "l")
legend("topright", title="Modèle de similarité :",legend=c("Géométrie", "Proximité"),
       col=c("red", "blue"), lty=c(1,1), cex=0.8)
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
seeds = 200
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
sLsP = 25
sLdP = 22
dLdP = 15
TP = sLsP+sLdP
TN = 50-dLdP
FP = dLdP
FN = 50-sLsP+50-sLdP
VPP = TP/(TP+FP)
VPP
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
##dist ## attention dt a prendre a partir d protocol_classification_feature
dist_lig_sameP_sameL = NULL
dist_lig_diffP_sameL = NULL
dist_lig_diffP_diffL = NULL
for (i in 1:nrow(pocket_samePsameL_50)) {
  print(i)
  dist_lig_sameP_sameL = c(dist_lig_sameP_sameL, 1-dist(rbind(dt_fpocket[pocket_samePsameL_50[i,2],features],dt[pocket_samePsameL_50[i,3],features])))
  dist_lig_diffP_sameL = c(dist_lig_diffP_sameL, 1-dist(rbind(dt_fpocket[pocket_diffPsameL_50[i,2],features],dt[pocket_diffPsameL_50[i,3],features])))
  dist_lig_diffP_diffL = c(dist_lig_diffP_diffL, 1-dist(rbind(dt_fpocket[pocket_diffPdiffL_50[i,2],features],dt[pocket_diffPdiffL_50[i,3],features])))
  
  #dist_lig_sameP_sameL = c(dist_lig_sameP_sameL, predict.glm(dt_predict.glm.step, newdata = sqrt((dt_fpocket[pocket_samePsameL_50[i,2],]-dt[pocket_samePsameL_50[i,3],])**2), type = "response"))
  #dist_lig_diffP_sameL = c(dist_lig_diffP_sameL, predict.glm(dt_predict.glm.step, newdata = sqrt((dt_fpocket[pocket_diffPsameL_50[i,2],]-dt[pocket_diffPsameL_50[i,3],])**2), type = "response"))
  #dist_lig_diffP_diffL = c(dist_lig_diffP_diffL, predict.glm(dt_predict.glm.step, newdata = sqrt((dt_fpocket[pocket_diffPdiffL_50[i,2],]-dt[pocket_diffPdiffL_50[i,3],])**2), type = "response"))
}
#
#dt = as.data.frame(dt)
#dt_overlap = as.data.frame(dt_overlap)
name_dt_overlap = rownames(dt_overlap)
#
y_predict = predict.glm(dt_predict.glm.step, newdata = sqrt((dt_overlap[name_dt_overlap,]-dt[name_dt_overlap,])**2), type = "response")
#
y_predict = NULL
for (i in 1:nrow(dt_overlap)) {
  print(i)
  y_predict = c(y_predict, predict.glm(dt_predict.glm.step, newdata = sqrt((dt_overlap[name_dt_overlap[i],]-dt[name_dt_overlap[i],])**2), type = "response"))
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
y_predict = c(dist_lig_diffP_sameL,
              dist_lig_diffP_diffL)
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
#SO_50 = read.table("../data/SO_desc_fpocket_all.txt",row.names = 1,sep=",")
SO_50[,1] = as.numeric(levels(SO_50[,1] ))[SO_50[,1]]
SO_50 = na.omit(SO_50)
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
sP_50 = cbind(SO_50[name_dt_overlap],y_predict)
sP_50 = as.data.frame(sP_50)
colnames(sP_50) = c("SO","ypred")
cor.test(sP_50$SO,sP_50$ypred)
###
g <- ggplot(sP_50, aes(SO, ypred)) +
  geom_count() +
  geom_smooth(method="lm", se=F) + xlim(0,0.6)+ylim(0,1)


ggMarginal(g, type = "histogram", fill="transparent")

p1 <- ggplot(sP_50, aes(SO,ypred)) +
  geom_smooth(method="lm", se=F) + xlim(0,0.8)+ylim(0,1)
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
  #return(dist(rbind(vec1,vec2)))
  #return(dist_fuzcav_ph(as.integer(vec1),as.integer(vec2)))
  return(dist_fuzcav_ph_norm(as.integer(vec1),as.integer(vec2)))
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

##dist ## attention dt a prendre a partir d protocol_classification_feature
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

### DISTRIBUTION Y_PRED DANS KMEANS INTERET ###
sink('../results/kmeans_results_reglog_article/output_50.txt')
for (seeds in c(5,10,20,30,40,50,60,70,80,90,100,200,300,400,500)) {
  for (seuil in c(0.5,0.6,0.7,0.8,0.9,1)) {#0.5,0.6,0.7,0.8,
    kmean_results = read.table(file = paste0(paste0("../results/kmeans_results_reglog_article/kmeans_reglog_seeds",seeds),"_clusters.txt"), sep = ",", header = F, row.names = 1)
    kmean_centroids = read.table(file = paste0(paste0("../results/kmeans_results_reglog_article/kmeans_reglog_seeds",seeds),"_means.txt"), sep = ",", header = F)
    colnames(kmean_centroids) = colnames(dt)
    #
    y_predict_30 = NULL
    y_predict_30_mean = NULL
    N_y_predict_0.8 = NULL
    N_y_predict_0.9 = NULL
    N_y_predict_0.95 = NULL
    N_y_predict = NULL
    #sameLsameP
    y_predict_sameLsameP = NULL
    nclust_predict_sameLsameP = NULL
    all_predicted_pocket = NULL
    query_predicted_pocket = NULL
    for (i in 1:nrow(pocket_samePsameL_50)) {
      #print(i)
      #pock2 found in cluster
      pock1_cluster = kmean_results[pocket_samePsameL_50[i,3],1]
      #pock2 predicted
      pock2_clusters = NULL
      dt_fpocket_expanded <- dt_fpocket[rep(pocket_samePsameL_50[i,2], nrow(kmean_centroids)),features]
      rownames(dt_fpocket_expanded) = rownames(kmean_centroids)
      pock2_clusters = predict.glm(dt_predict.glm.step, newdata = sqrt((kmean_centroids[,features]-dt_fpocket_expanded[,features])**2), type = "response")
      #
      pock2_cluster = which(pock2_clusters > seuil|pock2_clusters == max(pock2_clusters))-1
      nclust_predict_sameLsameP = c(nclust_predict_sameLsameP,length(pock2_cluster))
      
      #sur toutes les poches
      dt_fpocket_expanded <- dt_fpocket[rep(pocket_samePsameL_50[i,2], nrow(dt)),features]
      rownames(dt_fpocket_expanded) = rownames(dt)
      predicted_pocket = predict.glm(dt_predict.glm.step, newdata = sqrt((dt[,features]-dt_fpocket_expanded[,features])**2), type = "response")
      names_interest = rownames(kmean_results)[which(is.element(kmean_results[,1],pock2_cluster))]
      all_predicted_pocket = c(all_predicted_pocket,predicted_pocket[names_interest])
      query_predicted_pocket = c(query_predicted_pocket,predicted_pocket[pocket_samePsameL_50[i,3]])
      #
      N_y_predict = c(N_y_predict,length(predicted_pocket[names_interest]))
      y_predict_30 = c(y_predict_30, sort(predicted_pocket[names_interest], decreasing = TRUE)[1:30])
      y_predict_30_mean = c(y_predict_30_mean,mean(sort(predicted_pocket[names_interest], decreasing = TRUE)[1:30]))
      N_y_predict_0.8 = c(N_y_predict_0.8, length(which(predicted_pocket[names_interest]>0.8)))
      N_y_predict_0.9 = c(N_y_predict_0.9, length(which(predicted_pocket[names_interest]>0.9)))
      N_y_predict_0.95 = c(N_y_predict_0.95, length(which(predicted_pocket[names_interest]>0.95)))
      #
      if(is.element(pock1_cluster,pock2_cluster)) {
        y_predict_sameLsameP = c(y_predict_sameLsameP,1)
    
      } else {
        y_predict_sameLsameP = c(y_predict_sameLsameP,0)
      }
    }
    
    
    #sameLdiffP
    y_predict_sameLdiffP = NULL
    nclust_predict_sameLdiffP = NULL
    all_predicted_pocket = NULL
    query_predicted_pocket = NULL
    for (i in 1:nrow(pocket_diffPsameL_50)) {
      #print(i)
      #pock2 found in cluster
      pock1_cluster = kmean_results[pocket_diffPsameL_50[i,3],1]
      #pock2 predicted
      pock2_clusters = NULL
      dt_fpocket_expanded <- dt_fpocket[rep(pocket_diffPsameL_50[i,2], nrow(kmean_centroids)),features]
      rownames(dt_fpocket_expanded) = rownames(kmean_centroids)
      pock2_clusters = predict.glm(dt_predict.glm.step, newdata = sqrt((kmean_centroids[,features]-dt_fpocket_expanded[,features])**2), type = "response")
      #
      pock2_cluster = which(pock2_clusters > seuil|pock2_clusters == max(pock2_clusters))-1
      nclust_predict_sameLdiffP = c(nclust_predict_sameLdiffP,length(pock2_cluster))
      
      #sur toutes les poches
      dt_fpocket_expanded <- dt_fpocket[rep(pocket_diffPsameL_50[i,2], nrow(dt)),features]
      rownames(dt_fpocket_expanded) = rownames(dt)
      predicted_pocket = predict.glm(dt_predict.glm.step, newdata = sqrt((dt[,features]-dt_fpocket_expanded[,features])**2), type = "response")
      names_interest = rownames(kmean_results)[which(is.element(kmean_results[,1],pock2_cluster))]
      all_predicted_pocket = c(all_predicted_pocket,predicted_pocket[names_interest])
      query_predicted_pocket = c(query_predicted_pocket,predicted_pocket[pocket_diffPsameL_50[i,3]])
      #
      N_y_predict = c(N_y_predict,length(predicted_pocket[names_interest]))
      y_predict_30 = c(y_predict_30, sort(predicted_pocket[names_interest], decreasing = TRUE)[1:30])
      y_predict_30_mean = c(y_predict_30_mean,mean(sort(predicted_pocket[names_interest], decreasing = TRUE)[1:30]))
      N_y_predict_0.8 = c(N_y_predict_0.8, length(which(predicted_pocket[names_interest]>0.8)))
      N_y_predict_0.9 = c(N_y_predict_0.9, length(which(predicted_pocket[names_interest]>0.9)))
      N_y_predict_0.95 = c(N_y_predict_0.95, length(which(predicted_pocket[names_interest]>0.95)))
      #
      if(is.element(pock1_cluster,pock2_cluster)) {
        y_predict_sameLdiffP = c(y_predict_sameLdiffP,1)
        #print("pock apprise:")
        #print(pocket_diffPsameL_50[i,3])
        #print("pock non apprise:")
        #print(pocket_diffPsameL_50[i,2])
        #print("-----------------")
      } else {
        y_predict_sameLdiffP = c(y_predict_sameLdiffP,0)
      }
    }
    
    #diffLdiffP
    y_predict_diffLdiffP = NULL
    nclust_predict_diffLdiffP = NULL
    all_predicted_pocket = NULL
    query_predicted_pocket = NULL
    for (i in 1:nrow(pocket_diffPdiffL_50)) {
      #print(i)
      #pock2 found in cluster
      pock1_cluster = kmean_results[pocket_diffPdiffL_50[i,3],1]
      #pock2 predicted
      pock2_clusters = NULL
      dt_fpocket_expanded <- dt_fpocket[rep(pocket_diffPdiffL_50[i,2], nrow(kmean_centroids)),features]
      rownames(dt_fpocket_expanded) = rownames(kmean_centroids)
      pock2_clusters = predict.glm(dt_predict.glm.step, newdata = sqrt((kmean_centroids[,features]-dt_fpocket_expanded[,features])**2), type = "response")
      #
      pock2_cluster = which(pock2_clusters > seuil|pock2_clusters == max(pock2_clusters))-1
      nclust_predict_diffLdiffP = c(nclust_predict_diffLdiffP,length(pock2_cluster))
      
      #sur toutes les poches
      dt_fpocket_expanded <- dt_fpocket[rep(pocket_diffPdiffL_50[i,2], nrow(dt)),features]
      rownames(dt_fpocket_expanded) = rownames(dt)
      predicted_pocket = predict.glm(dt_predict.glm.step, newdata = sqrt((dt[,features]-dt_fpocket_expanded[,features])**2), type = "response")
      names_interest = rownames(kmean_results)[which(is.element(kmean_results[,1],pock2_cluster))]
      all_predicted_pocket = c(all_predicted_pocket,predicted_pocket[names_interest])
      query_predicted_pocket = c(query_predicted_pocket,predicted_pocket[pocket_diffPdiffL_50[i,3]])
      #
      N_y_predict = c(N_y_predict,length(predicted_pocket[names_interest]))
      y_predict_30 = c(y_predict_30, sort(predicted_pocket[names_interest], decreasing = TRUE)[1:30])
      y_predict_30_mean = c(y_predict_30_mean,mean(sort(predicted_pocket[names_interest], decreasing = TRUE)[1:30]))
      N_y_predict_0.8 = c(N_y_predict_0.8, length(which(predicted_pocket[names_interest]>0.8)))
      N_y_predict_0.9 = c(N_y_predict_0.9, length(which(predicted_pocket[names_interest]>0.9)))
      N_y_predict_0.95 = c(N_y_predict_0.95, length(which(predicted_pocket[names_interest]>0.95)))
      #
      if(is.element(pock1_cluster,pock2_cluster)) {
        y_predict_diffLdiffP = c(y_predict_diffLdiffP,1)
        #print("pock apprise:")
        #print(pocket_diffPdiffL_50[i,3])
        #print("pock non apprise:")
        #print(pocket_diffPdiffL_50[i,2])
        #print("---------c l dif--------")
        
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
    print("----------------------")
    print(paste("seeds", seeds))
    print(paste("seuil", seuil))
    print(paste("sLsP",length(which(y_predict_sameLsameP == 1))))
    print(paste("sLdP",length(which(y_predict_sameLdiffP == 1))))
    print(paste("dLdP",length(which(y_predict_diffLdiffP == 1))))
    print(paste("N1",mean(c(nclust_predict_sameLsameP,
                            nclust_predict_sameLdiffP,
                            nclust_predict_diffLdiffP))))
    #print(paste("N2",mean(nclust_predict_sameLdiffP)))
    #print(paste("N3",mean(nclust_predict_diffLdiffP)))
    
    print(paste("se", Se))
    print(paste("sp", Sp))
    #
    #VPP = TP/(TP+FP)
    #print(paste("VPP", VPP))
    print(paste("mcc", mcc(preds = y_predict, actuals = y_true)))
    #
    print(paste("mean(y_predict_30)",mean(y_predict_30)))
    print(paste("sd(y_predict_30)",sd(y_predict_30)))
    print(paste("mean(N_y_predict_0.8)",mean(N_y_predict_0.8)))
    print(paste("sd(N_y_predict_0.8)",sd(N_y_predict_0.8)))
    print(paste("mean(N_y_predict_0.9)",mean(N_y_predict_0.9)))
    print(paste("sd(N_y_predict_0.9)",sd(N_y_predict_0.9)))
    print(paste("mean(N_y_predict_0.95)",mean(N_y_predict_0.95)))
    print(paste("sd(N_y_predict_0.95)",sd(N_y_predict_0.95)))
    print(paste("mean(N_y_predict)",mean(N_y_predict)))
    print(paste("sd(N_y_predict)",sd(N_y_predict)))
    #
  }
}
sink()

#
y_predict_30_mean
length(y_predict_30_mean)
which(N_y_predict_0.9 == 0)
library(ggplot2)
dataset = data.frame(x=c(1:150,1:150),seuil = c(N_y_predict_0.95,N_y_predict_0.9),
                     Seuil =c(rep("0.95",150),
                             rep("0.90",150)))

b <- ggplot(dataset, aes(x = x, y = seuil,fill = Seuil)) +
  geom_bar(stat = "identity") +
  ylab("Effectif de poches candidates") +
  xlab("150 poches")
b
barplot(N_y_predict_0.95)
#
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
#
mcc(preds = y_predict, actuals = y_true)
#
length(which(y_predict_sameLsameP == 1))
length(which(y_predict_sameLdiffP == 1))
length(which(y_predict_diffLdiffP == 1))
mean(nclust_predict_sameLsameP)
mean(nclust_predict_sameLdiffP)
mean(nclust_predict_diffLdiffP)
#
mean(y_predict_30)
sd(y_predict_30)
mean(N_y_predict_0.8)
mean(N_y_predict_0.9)
mean(N_y_predict_0.95)
mean(N_y_predict)
#library(mltools)
sLsP = 25
sLdP = 22
dLdP = 15
TP = sLsP+sLdP
TN = 50-dLdP
FP = dLdP
FN = 50-sLsP+50-sLdP
VPP = TP/(TP+FP)
VPP
MCC = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
MCC
#
hist(all_predicted_pocket)
hist(query_predicted_pocket)

#library(sm)
sm.density.compare(c(all_predicted_pocket,
                     query_predicted_pocket),
                   c(rep(1,length(all_predicted_pocket)),
                     rep(2,length(query_predicted_pocket))
                   ),
                   model = "none", xlim=c(0,1)
                   , xlab = "prediction des poches dans les clusters"
                   , main = "densite des poches predites dans les cluster (k 200 seuil 0.5) pour paires poches sLsP")
length(all_predicted_pocket)

dt_fpocket_expanded <- dt_fpocket[rep(pocket_samePsameL_50[i,2], nrow(dt)),features]
rownames(dt_fpocket_expanded) = rownames(dt)
predicted_pocket = predict.glm(dt_predict.glm.step, newdata = sqrt((dt[,features]-dt_fpocket_expanded[,features])**2), type = "response")
names_interest = rownames(kmean_results)[which(is.element(kmean_results[,1],pock2_cluster))]
all_predicted_pocket = c(all_predicted_pocket,predicted_pocket[names_interest])
query_predicted_pocket = c(query_predicted_pocket,predicted_pocket[pocket_samePsameL_50[i,2]])


### QUALITE KMEANS ####
seeds = 30
kmean_results_1 = read.table(file = paste0(paste0("../data/kmeans_results_reglog/kmeans_reglog_seeds",seeds),"_clusters.txt"), sep = ",", header = F, row.names = 1)
kmean_results_2 = read.table(file = paste0(paste0("../results/kmeans_results_reglog_PproxMfpo/kmeans_reglog_seeds",seeds),"_clusters.txt"), sep = ",", header = F, row.names = 1)

names_kmean_results_1 = rownames(kmean_results_1)
names_kmean_results_2 = rownames(kmean_results_2)

kmean_results_1 = as.vector(kmean_results_1[,1])
kmean_results_2 = as.vector(kmean_results_2[,1])

names(kmean_results_1) = names_kmean_results_1
names(kmean_results_2) = names_kmean_results_2

#
matrix.sort <- function(matrix) {
  if (nrow(matrix) != ncol(matrix)) stop("Not diagonal")
  if(is.null(rownames(matrix))) rownames(matrix) <- 1:nrow(matrix)
  row.max <- apply(matrix,1,which.max)
  if(all(table(row.max) != 1)) stop("Ties cannot be resolved")
  return(matrix[names(sort(row.max)),])
}

tables_kmean_results = table(factor(kmean_results_1,levels = 0:seeds),factor(kmean_results_2,levels = 0:seeds))
library(pheatmap)
heatmap(matrix.sort(tables_kmean_results),Colv = NA, Rowv = NA)

library(gplots)
heatmap.2(matrix.sort(tables_kmean_results), col = topo.colors(seeds),dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',density.info = "none")

library(ggplot2)
ggplot(data = as.data.frame(matrix.sort(tables_kmean_results)), aes(x=kmean_results_1, y=kmean_results_2, fill=value)) + 
  geom_tile()


### HCLUST CENTROIDES ###
library(data.tree)
##
seeds = 200
#
path = "../data/kmeans_results_reglog/"
path = "../results/kmeans_results_reglog_PproxMfpo/"
#
colnames_dt = colnames(dt)
#
kmean_results = read.table(file = paste0(paste0(paste0(path,"kmeans_reglog_seeds"),seeds),"_clusters.txt"), sep = ",", header = F, row.names = 1)
kmean_centroids = read.table(file = paste0(paste0("../results/kmeans_results_reglog_PproxMfpo/kmeans_reglog_seeds",seeds),"_means.txt"), sep = ",", header = F)
colnames(kmean_centroids) = colnames_dt
kmean_size = as.vector(table(kmean_results[,1]))
kmean_SSE =  read.table(file = paste0(paste0(paste0(path,"kmeans_reglog_SSE_seeds"),seeds),".txt"), sep = ",",header = T)
#
mat_dist = matrix(0,nrow = nrow(kmean_centroids),ncol=nrow(kmean_centroids))
colnames(mat_dist) = 1:nrow(kmean_centroids)
rownames(mat_dist) = 1:nrow(kmean_centroids)
print(nrow(mat_dist))
pock_dist = NULL
for (i in 1:nrow(mat_dist)) {
  print(i)
  kmean_centroids_expanded <- kmean_centroids[rep(i, nrow(mat_dist)),features]
  rownames(kmean_centroids_expanded) = rownames(mat_dist)
  mat_dist[i,] = 1-predict.glm(dt_predict.glm.step, newdata = sqrt((kmean_centroids[,features]-kmean_centroids_expanded[,features])**2), type = "response")
  pock_dist = c(pock_dist, predict.glm(dt_predict.glm.step, newdata = sqrt((kmean_centroids[i,features]-dt[1,features])**2), type = "response"))
}
hclust_kmean_centroids_200_prox = hclust(as.dist(mat_dist),method = "ward.D2")
pdf("../rapport/fig/hclust_geom_200.pdf",width=40, height=15)
plot(as.dendrogram(hclust_kmean_centroids_200_prox),type = "rectangle", ylab = "Height", cex.lab = 0.7, cex.main = 8, cex.axis = 2.5)
dev.off()
###
dt_centers.dend = as.dendrogram(hclust_kmean_centroids_200_prox)
### PLOT ###
library(dendextend)
library(colorspace)
library(plotfunctions)
#label name
#change label name
#
centroid_green = which(labels(dt_centers.dend) == as.integer(names(which(pock_dist == max(pock_dist)))))
lab_cluster_dend = as.integer(labels(dt_centers.dend))
labels(dt_centers.dend) <- paste0(lab_cluster_dend, 
                                  paste0(";Dist:",
                                         round(pock_dist[as.character(lab_cluster_dend)], 2)))
labels(dt_centers.dend)[centroid_green] <- paste0("POCHE->;GRP",
                                                  paste0(lab_cluster_dend[centroid_green], 
                                                         paste0(";Dist:",
                                                                round(max(pock_dist), 2))))
#label color
color.gradient <- function(x, colors=c("black","gold","green"), colsteps=100) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),0.9, length.out=colsteps)) ] )
}
c_colors = color.gradient(pock_dist)
color.df<-data.frame(COLOR_VALUE=pock_dist, color.name=color.gradient(pock_dist))
labels_colors(dt_centers.dend) <- as.character(color.df[as.character(order.dendrogram(dt_centers.dend)), "color.name"])
#
pdf("../results/results_rapport/test.pdf",width=40, height=15)
par(cex=1, mar=c(13, 9, 5, 5))
plot(dt_centers.dend, type = "rectangle", ylab = "Height", cex.lab = 2, cex.main = 3,
     main = "Distance de la poche ")
legend("right", title = "Distance de la poche",
       legend = c(paste0("max:", round(max(pock_dist),2)),
                  rep("",1),
                  paste0("mean:", round(mean(pock_dist),2)),
                  rep("",1),
                  paste0("min:", round(min(pock_dist),2))), pt.cex = 4, cex = 2, bty = "n")
legend("topleft", title = "10 poches les plus proches:",
       legend = paste("P:", 
                      paste(names(sort(pock_dist, decreasing = TRUE)[1:10]),
                            paste("|D", round(sort(pock_dist, decreasing = TRUE)[1:10],2)))),
       pt.cex = 4, cex = 2, bty = "n")

gradientLegend(valRange=c(-14,14), color=color.gradient(sort(pock_dist)), 
               border.col=alpha('gray'), side=4, pos.num = 3, inside = FALSE, length=.2, depth=.02)#pos=c(120,3.5,125,5.5),
dev.off()
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

### EXEMPLE RAPPORT  ###
# KAN
c("1L8T_KAN_A_1",
"6BFH_KAN_A_1",
"3U6T_KAN_A_1",
"4WQL_KAN_A_1")
kmean_results[c("1L8T_KAN_A_1",
                "6BFH_KAN_A_1",
                "3U6T_KAN_A_1",
                "4WQL_KAN_A_1"),]

which(names_prot == "1LB2")
rownames(dt)[which(names_prot == "1V2Q")]


c("5ISX_PNS_A_1","1QJC_PNS_B_1")

c("4KYK_IMN_B_1",
"4COX_IMN_A_1",
"4JQ4_IMN_A_1")

kmean_results[c("4KYK_IMN_B_1",
                "4COX_IMN_A_1",
                "4JQ4_IMN_A_1"),]
#
c("2Y00_Y01_A_1","5OSC_Y01_B_3")
#
4UE3_NFV_L_1
2IW8_4SP_A_2
1FJ4_TLM_C_3
1M1Y_HCA_A_4
pocket_diffPsameL_50[which(pocket_diffPsameL_50[,3] == "1LB2_CMP_A_1"),]
#
rownames(kmean_results)[which(kmean_results == 19)]
###
#modele
load("../results/kmeans_results_reglog/model.glm.step.Rdata")
#
load("../results/kmeans_results_reglog/model.glm.step_overlap.Rdata")
#
predict.glm(dt_predict.glm.step, newdata = sqrt((dt["5ISX_PNS_A_1",features]-dt["1QJC_PNS_B_1",features])**2), type = "response")
#
T1<-Sys.time()

dt_fpocket=dt["1QJC_PNS_B_1",features]
dt_fpocket_expanded <- dt_fpocket[rep("1QJC_PNS_B_1", nrow(dt)),features]
rownames(dt_fpocket_expanded) = rownames(dt)
predicted_pocket = predict.glm(dt_predict.glm.step, newdata = sqrt((dt[,features]-dt_fpocket_expanded[,features])**2), type = "response")
sort(predicted_pocket,decreasing = T)[1:20]
mean(sort(predicted_pocket,decreasing = T)[1:30])
T2<-Sys.time()
T2-T1




###TEMPS recherche###

T1<-Sys.time()
seeds = 200
#
path = "../results/kmeans_results_reglog_PproxMfpo/"
#
colnames_dt = colnames(dt)
#
kmean_results = read.table(file = paste0(paste0(paste0(path,"kmeans_reglog_seeds"),seeds),"_clusters.txt"), sep = ",", header = F, row.names = 1)
kmean_centroids = read.table(file = paste0(paste0("../results/kmeans_results_reglog_PproxMfpo/kmeans_reglog_seeds",seeds),"_means.txt"), sep = ",", header = F)
colnames(kmean_centroids) = colnames(dt)
pock1_cluster = kmean_results["1QJC_PNS_B_1",1]
#pock2 predicted
pock2_clusters = NULL
dt_fpocket=dt["5ISX_PNS_A_1",features]
T2<-Sys.time()
T2-T1


T1<-Sys.time()
###NEW POCKET
#PATH
path_covid_6lu7="../data/pockets_MD_NS1/Res_pocketConf0101-p0_atm/pocketConf101-0_atm.des"

dt_new_pocket_des = read.table(path_covid_6lu7, row.names = 1)
dt_fpocket = data.frame(A = dt_new_pocket_des["pocket_A",1],
                        C = dt_new_pocket_des["pocket_C",1],
                        C_ATOM = dt_new_pocket_des["pocket_C_ATOM",1],
                        C_RESIDUES = dt_new_pocket_des["pocket_C_RESIDUES",1],
                        charge = dt_new_pocket_des["pocket_charge",1],
                        CONVEX.SHAPE_COEFFICIENT = dt_new_pocket_des["pocket_CONVEX-SHAPE_COEFFICIENT",1],
                        D = dt_new_pocket_des["pocket_D",1],
                        DIAMETER_HULL = dt_new_pocket_des["pocket_DIAMETER_HULL",1],
                        E = dt_new_pocket_des["pocket_E",1],
                        F = dt_new_pocket_des["pocket_F",1],
                        FACE = dt_new_pocket_des["pocket_FACE",1],
                        G = dt_new_pocket_des["pocket_G",1],
                        H = dt_new_pocket_des["pocket_H",1],
                        hydrophobic_kyte = dt_new_pocket_des["pocket_hydrophobic_kyte",1],
                        hydrophobicity = dt_new_pocket_des["pocket_hydrophobicity",1],
                        I = dt_new_pocket_des["pocket_I",1],
                        INERTIA_1 = dt_new_pocket_des["pocket_INERTIA_1",1],
                        INERTIA_2 = dt_new_pocket_des["pocket_INERTIA_2",1],
                        INERTIA_3 = dt_new_pocket_des["pocket_INERTIA_3",1],
                        K = dt_new_pocket_des["pocket_K",1],
                        L = dt_new_pocket_des["pocket_L",1],
                        M = dt_new_pocket_des["pocket_M",1],
                        N = dt_new_pocket_des["pocket_N",1],
                        P = dt_new_pocket_des["pocket_P",1],
                        p_aliphatic_residues = dt_new_pocket_des["pocket_p_aliphatic_residues",1],
                        p_aromatic_residues = dt_new_pocket_des["pocket_p_aromatic_residues",1],
                        p_C_atom = dt_new_pocket_des["pocket_p_C_atom",1],
                        p_Car_atom = dt_new_pocket_des["pocket_p_Car_atom",1],
                        p_carbone_atom = dt_new_pocket_des["pocket_p_carbone_atom",1],
                        p_Carg_atom = dt_new_pocket_des["pocket_p_Carg_atom",1],
                        p_Ccoo_atom = dt_new_pocket_des["pocket_p_Ccoo_atom",1],
                        p_Cgln_atom = dt_new_pocket_des["pocket_p_Cgln_atom",1],
                        p_charged_residues = dt_new_pocket_des["pocket_p_charged_residues",1],
                        p_hyd_atom = dt_new_pocket_des["pocket_p_hyd_atom",1],
                        p_hydrophobic_atom = dt_new_pocket_des["pocket_p_hydrophobic_atom",1],
                        p_hydrophobic_residues = dt_new_pocket_des["pocket_p_hydrophobic_residues",1],
                        p_main_chain_atom = dt_new_pocket_des["pocket_p_main_chain_atom",1],
                        p_N_atom = dt_new_pocket_des["pocket_p_N_atom",1],
                        p_ND1_atom = dt_new_pocket_des["pocket_p_ND1_atom",1],
                        p_NE2_atom = dt_new_pocket_des["pocket_p_NE2_atom",1],
                        p_negative_residues = dt_new_pocket_des["pocket_p_negative_residues",1],
                        p_nitrogen_atom = dt_new_pocket_des["pocket_p_nitrogen_atom",1],
                        p_Nlys_atom = dt_new_pocket_des["pocket_p_Nlys_atom",1],
                        p_Ntrp_atom = dt_new_pocket_des["pocket_p_Ntrp_atom",1],
                        p_O_atom = dt_new_pocket_des["pocket_p_O_atom",1],
                        p_Ocoo_atom = dt_new_pocket_des["pocket_p_Ocoo_atom",1],
                        p_Ooh_atom = dt_new_pocket_des["pocket_p_Ooh_atom",1],
                        p_Otyr_atom = dt_new_pocket_des["pocket_p_Otyr_atom",1],
                        p_oxygen_atom = dt_new_pocket_des["pocket_p_oxygen_atom",1],
                        p_polar_residues = dt_new_pocket_des["pocket_p_polar_residues",1],
                        p_positive_residues = dt_new_pocket_des["pocket_p_positive_residues",1],
                        p_S_atom = dt_new_pocket_des["pocket_p_S_atom",1],
                        p_side_chain_atom = dt_new_pocket_des["pocket_p_side_chain_atom",1],
                        p_small_residues = dt_new_pocket_des["pocket_p_small_residues",1],
                        p_sulfur_atom = dt_new_pocket_des["pocket_p_sulfur_atom",1],
                        p_tiny_residues = dt_new_pocket_des["pocket_p_tiny_residues",1],
                        PCI = dt_new_pocket_des["pocket_PCI",1],
                        polarity = dt_new_pocket_des["pocket_polarity",1],
                        PSI = dt_new_pocket_des["pocket_PSI",1],
                        Q = dt_new_pocket_des["pocket_Q",1],
                        R = dt_new_pocket_des["pocket_R",1],
                        RADIUS_CYLINDER = dt_new_pocket_des["pocket_RADIUS_CYLINDER",1],
                        RADIUS_HULL = dt_new_pocket_des["pocket_RADIUS_HULL",1],
                        S = dt_new_pocket_des["pocket_S",1],
                        SMALLEST_SIZE = dt_new_pocket_des["pocket_SMALLEST_SIZE",1],
                        SURFACE_HULL = dt_new_pocket_des["pocket_SURFACE_HULL",1],
                        T = dt_new_pocket_des["pocket_T",1],
                        V = dt_new_pocket_des["pocket_V",1],
                        W = dt_new_pocket_des["pocket_W",1],
                        X._ATOM_CONVEXE = dt_new_pocket_des["pocket_%_ATOM_CONVEXE",1],
                        Y = dt_new_pocket_des["pocket_Y",1],
                        VOLUME_HULL = dt_new_pocket_des["pocket_VOLUME_HULL",1]
)

#
scaled_center_dt_t = read.table(file = "../results/scale/scaled_center_dt72clean.Rdata", col.names = F, row.names = 1)
scaled_scale_dt_t = read.table(file = "../results/scale/scaled_scale_dt72clean.Rdata", col.names = F, row.names = 1)
scaled_center_dt = scaled_center_dt_t[,1]
names(scaled_center_dt) = rownames(scaled_center_dt_t)
scaled_scale_dt = scaled_scale_dt_t[,1]
names(scaled_scale_dt) = rownames(scaled_scale_dt_t)
#
dt_fpocket = as.data.frame(scale(dt_fpocket[,colnames(dt_fpocket)], scaled_center_dt[colnames(dt_fpocket)],
                                                                    scaled_scale_dt[colnames(dt_fpocket)]))

rownames(dt_fpocket)
colnames(dt_fpocket)
setdiff(features,colnames(dt_fpocket))
##
dt_fpocket_expanded <- dt_fpocket[rep("1", nrow(kmean_centroids)),features]
rownames(dt_fpocket_expanded) = rownames(kmean_centroids)
pock2_clusters = predict.glm(dt_predict.glm.step, newdata = sqrt((kmean_centroids[,features]-dt_fpocket_expanded[,features])**2), type = "response")
name_dt = rownames(kmean_results)[which(is.element(kmean_results[,1],which(pock2_clusters == max(pock2_clusters) | pock2_clusters>0.7)-1))]
#sur toutes les poches
dt_fpocket_expanded <- dt_fpocket[rep("1", nrow(dt[name_dt,])),features]
rownames(dt_fpocket_expanded) = rownames(dt[name_dt,])
predicted_pocket = predict.glm(dt_predict.glm.step, newdata = sqrt((dt[name_dt,features]-dt_fpocket_expanded[,features])**2), type = "response")
sort(predicted_pocket[name_dt], decreasing = TRUE)[1:30]
mean(sort(predicted_pocket[name_dt], decreasing = TRUE)[1:30])
sd(sort(predicted_pocket[name_dt], decreasing = TRUE)[1:30])

names(sort(predicted_pocket[name_dt], decreasing = TRUE)[1:30])
names_prot = sapply(strsplit(names(sort(predicted_pocket[name_dt], decreasing = TRUE)[1:30]), "_"), "[", 1)
names_ligand = sapply(strsplit(names(sort(predicted_pocket[name_dt], decreasing = TRUE)[1:30]), "_"), "[", 2)

length(unique(names_prot))
length(unique(names_ligand))

T2<-Sys.time()
T2-T1



names_interest = rownames(kmean_results)[which(is.element(kmean_results[,1],pock2_cluster))]
all_predicted_pocket = c(all_predicted_pocket,predicted_pocket[names_interest])
query_predicted_pocket = c(query_predicted_pocket,predicted_pocket[pocket_samePsameL_50[i,2]])
###enzymes###
pdb_names_hydrolases = scan("../data/enzymes_pdb_names/hydrolases.txt", character(), quote = "")
pdb_names_hydrolases = sub(",","",pdb_names_hydrolases)
length(intersect(names_prot,pdb_names_hydrolases))
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
