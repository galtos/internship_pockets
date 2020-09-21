#### Analyse K means PYTHON ####
#pocket_diffPsameL_50
load("../data/data_structure_comparaison/pocket_pos_100.Rdata")
#pocket_neg_100
load("../data/data_structure_comparaison/pocket_neg_100.Rdata")

length(unique(as.character(pocket_pos_100)))
length(unique(as.character(pocket_neg_100)))
## SCORE OVERLAP 100 POCKETS ###
value_so_pos = SO_fpocket[which(is.element(SO_fpocket$V1,as.character(pocket_pos_100))),2]
value_so_neg = SO_fpocket[which(is.element(SO_fpocket$V1,as.character(pocket_neg_100))),2]

boxplot(c(value_so_pos,value_so_neg))
hist(c(value_so_pos,value_so_neg))
length(which(c(value_so_pos,value_so_neg) < 0.2286))

summary(SO_fpocket)
summary(SO_fpocket[c(as.character(pocket_pos_100),as.character(pocket_neg_100)),])
t.test(na.omit(SO_fpocket$V2))
#
rownames(SO_fpocket) = as.character(SO_fpocket$V1)
pocket_pos_100 = pocket_pos_100[which(SO_fpocket[as.character(pocket_pos_100),2] > 0.2286)]
pocket_neg_100 = pocket_neg_100[which(SO_fpocket[as.character(pocket_neg_100),2] > 0.2286)]
## POS and NEG
data_samePsameL_positive = read.csv("../data/data_structure_comparaison/data_samePsameL_dt72_pharmacophores_positive.csv", colClasses = "character")
data_diffPdiffL_negative = read.csv("../data/data_structure_comparaison/data_diffPdiffL_dt72_pharmacophores_negative.csv", colClasses = "character")
#
pocket_pos_100 = data_samePsameL_positive[which(is.element(data_samePsameL_positive$name_pock1,as.character(pocket_pos_100))),2:3]
pocket_neg_100 = data_diffPdiffL_negative[which(is.element(data_diffPdiffL_negative$name_pock1,as.character(pocket_neg_100))),2:3]

nrow(pocket_pos_100)
nrow(pocket_neg_100)

length(unique(pocket_pos_100$name_pock1))
length(unique(pocket_neg_100$name_pock1))

#### ANALYSIS ####
path = "../results/kmeans_results_reglog_article/"

#
Nseeds = c(5,10,20,30,40,50,60,70,80,90,100,200,300,400,500)
kmean_It = scan(file =  paste0(path,"/kmeans_reglog_TSS.txt"))#35354.82457577707#


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


### DISTRIBUTION Y_PRED DANS KMEANS INTERET ###
dt_fpocket = read.table("../data/data_PDB_72desc_overlap_all.txt", header = T, sep = "", row.names = 1, fill = TRUE)
#
scaled_center_dt_t = read.table(file = "../results/scale/scaled_center_dt72clean.Rdata", col.names = F, row.names = 1)
scaled_scale_dt_t = read.table(file = "../results/scale/scaled_scale_dt72clean.Rdata", col.names = F, row.names = 1)
scaled_center_dt = scaled_center_dt_t[,1]
names(scaled_center_dt) = rownames(scaled_center_dt_t)
scaled_scale_dt = scaled_scale_dt_t[,1]
names(scaled_scale_dt) = rownames(scaled_scale_dt_t)
#
dt_fpocket = as.data.frame(scale(dt_fpocket, scaled_center_dt, scaled_scale_dt))
# MODEL
load(file = "../results/kmeans_results_reglog/model_article_geom.glm.step.Rdata")
#library
library(mltools)
#Nseeds
sink('../results/kmeans_results_reglog_article/output_geom.txt')
for (seeds in c(100,200,300,400,500,50)) {
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
    for (i in 1:nrow(data_samePsameL_positive)) {
      #print(i)
      #pock2 found in cluster
      pock1_cluster = kmean_results[data_samePsameL_positive[i,3],1]
      #pock2 predicted
      pock2_clusters = NULL
      dt_fpocket_expanded <- dt_fpocket[rep(data_samePsameL_positive[i,2], nrow(kmean_centroids)),features]
      rownames(dt_fpocket_expanded) = rownames(kmean_centroids)
      pock2_clusters = predict.glm(dt_predict.glm.step, newdata = sqrt((kmean_centroids[,features]-dt_fpocket_expanded[,features])**2), type = "response")
      #
      pock2_cluster = which(pock2_clusters > seuil|pock2_clusters == max(pock2_clusters))-1
      nclust_predict_sameLsameP = c(nclust_predict_sameLsameP,length(pock2_cluster))
      
      #sur toutes les poches
      dt_fpocket_expanded <- dt_fpocket[rep(data_samePsameL_positive[i,2], nrow(dt)),features]
      rownames(dt_fpocket_expanded) = rownames(dt)
      predicted_pocket = predict.glm(dt_predict.glm.step, newdata = sqrt((dt[,features]-dt_fpocket_expanded[,features])**2), type = "response")
      names_interest = rownames(kmean_results)[which(is.element(kmean_results[,1],pock2_cluster))]
      all_predicted_pocket = c(all_predicted_pocket,predicted_pocket[names_interest])
      query_predicted_pocket = c(query_predicted_pocket,predicted_pocket[data_samePsameL_positive[i,3]])
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
    
    
    #diffLdiffP
    y_predict_diffLdiffP = NULL
    nclust_predict_diffLdiffP = NULL
    all_predicted_pocket = NULL
    query_predicted_pocket = NULL
    for (i in 1:nrow(data_diffPdiffL_negative)) {
      #print(i)
      #pock2 found in cluster
      pock1_cluster = kmean_results[data_diffPdiffL_negative[i,3],1]
      #pock2 predicted
      pock2_clusters = NULL
      dt_fpocket_expanded <- dt_fpocket[rep(data_diffPdiffL_negative[i,2], nrow(kmean_centroids)),features]
      rownames(dt_fpocket_expanded) = rownames(kmean_centroids)
      pock2_clusters = predict.glm(dt_predict.glm.step, newdata = sqrt((kmean_centroids[,features]-dt_fpocket_expanded[,features])**2), type = "response")
      #
      pock2_cluster = which(pock2_clusters > seuil|pock2_clusters == max(pock2_clusters))-1
      nclust_predict_diffLdiffP = c(nclust_predict_diffLdiffP,length(pock2_cluster))
      
      #sur toutes les poches
      dt_fpocket_expanded <- dt_fpocket[rep(data_diffPdiffL_negative[i,2], nrow(dt)),features]
      rownames(dt_fpocket_expanded) = rownames(dt)
      predicted_pocket = predict.glm(dt_predict.glm.step, newdata = sqrt((dt[,features]-dt_fpocket_expanded[,features])**2), type = "response")
      names_interest = rownames(kmean_results)[which(is.element(kmean_results[,1],pock2_cluster))]
      all_predicted_pocket = c(all_predicted_pocket,predicted_pocket[names_interest])
      query_predicted_pocket = c(query_predicted_pocket,predicted_pocket[data_diffPdiffL_negative[i,3]])
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
               rep(0, length(y_predict_diffLdiffP)))
    y_predict = c(y_predict_sameLsameP,
                  y_predict_diffLdiffP)
    table_clusters <- table(factor(y_predict,levels = 0:1), factor(y_true,levels=0:1))
    Se = table_clusters[2,2] / (table_clusters[2,2] + table_clusters[1,2])
    Sp = table_clusters[1,1] / (table_clusters[1,1] + table_clusters[2,1])
    print("----------------------")
    print(paste("seeds", seeds))
    print(paste("seuil", seuil))
    print(paste("sLsP",length(which(y_predict_sameLsameP == 1))))
    print(paste("dLdP",length(which(y_predict_diffLdiffP == 1))))
    print(paste("N1",mean(c(nclust_predict_sameLsameP,
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

dt_fpocket = as.data.frame(dt)
#sink('../results/kmeans_results_reglog_article/output_prox.txt')
for (seeds in c(100,200,300,400,500,50)) {
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
    for (i in 1:nrow(data_samePsameL_positive)) {
      #print(i)
      #pock2 found in cluster
      pock1_cluster = kmean_results[data_samePsameL_positive[i,3],1]
      #pock2 predicted
      pock2_clusters = NULL
      dt_fpocket_expanded <- dt_fpocket[rep(data_samePsameL_positive[i,2], nrow(kmean_centroids)),features]
      rownames(dt_fpocket_expanded) = rownames(kmean_centroids)
      pock2_clusters = predict.glm(dt_predict.glm.step, newdata = sqrt((kmean_centroids[,features]-dt_fpocket_expanded[,features])**2), type = "response")
      #
      pock2_cluster = which(pock2_clusters > seuil|pock2_clusters == max(pock2_clusters))-1
      nclust_predict_sameLsameP = c(nclust_predict_sameLsameP,length(pock2_cluster))
      
      #sur toutes les poches
      dt_fpocket_expanded <- dt_fpocket[rep(data_samePsameL_positive[i,2], nrow(dt)),features]
      rownames(dt_fpocket_expanded) = rownames(dt)
      predicted_pocket = predict.glm(dt_predict.glm.step, newdata = sqrt((dt[,features]-dt_fpocket_expanded[,features])**2), type = "response")
      names_interest = rownames(kmean_results)[which(is.element(kmean_results[,1],pock2_cluster))]
      all_predicted_pocket = c(all_predicted_pocket,predicted_pocket[names_interest])
      query_predicted_pocket = c(query_predicted_pocket,predicted_pocket[data_samePsameL_positive[i,3]])
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
    
    
    #diffLdiffP
    y_predict_diffLdiffP = NULL
    nclust_predict_diffLdiffP = NULL
    all_predicted_pocket = NULL
    query_predicted_pocket = NULL
    for (i in 1:nrow(data_diffPdiffL_negative)) {
      #print(i)
      #pock2 found in cluster
      pock1_cluster = kmean_results[data_diffPdiffL_negative[i,3],1]
      #pock2 predicted
      pock2_clusters = NULL
      dt_fpocket_expanded <- dt_fpocket[rep(data_diffPdiffL_negative[i,2], nrow(kmean_centroids)),features]
      rownames(dt_fpocket_expanded) = rownames(kmean_centroids)
      pock2_clusters = predict.glm(dt_predict.glm.step, newdata = sqrt((kmean_centroids[,features]-dt_fpocket_expanded[,features])**2), type = "response")
      #
      pock2_cluster = which(pock2_clusters > seuil|pock2_clusters == max(pock2_clusters))-1
      nclust_predict_diffLdiffP = c(nclust_predict_diffLdiffP,length(pock2_cluster))
      
      #sur toutes les poches
      dt_fpocket_expanded <- dt_fpocket[rep(data_diffPdiffL_negative[i,2], nrow(dt)),features]
      rownames(dt_fpocket_expanded) = rownames(dt)
      predicted_pocket = predict.glm(dt_predict.glm.step, newdata = sqrt((dt[,features]-dt_fpocket_expanded[,features])**2), type = "response")
      names_interest = rownames(kmean_results)[which(is.element(kmean_results[,1],pock2_cluster))]
      all_predicted_pocket = c(all_predicted_pocket,predicted_pocket[names_interest])
      query_predicted_pocket = c(query_predicted_pocket,predicted_pocket[data_diffPdiffL_negative[i,3]])
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
               rep(0, length(y_predict_diffLdiffP)))
    y_predict = c(y_predict_sameLsameP,
                  y_predict_diffLdiffP)
    table_clusters <- table(factor(y_predict,levels = 0:1), factor(y_true,levels=0:1))
    Se = table_clusters[2,2] / (table_clusters[2,2] + table_clusters[1,2])
    Sp = table_clusters[1,1] / (table_clusters[1,1] + table_clusters[2,1])
    print("----------------------")
    print(paste("seeds", seeds))
    print(paste("seuil", seuil))
    print(paste("sLsP",length(which(y_predict_sameLsameP == 1))))
    print(paste("dLdP",length(which(y_predict_diffLdiffP == 1))))
    print(paste("N1",mean(c(nclust_predict_sameLsameP,
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
#sink()

### DENSITY PLOT ###
##dist
dist_lig_diffP_diffL = NULL 
dist_lig_sameP_sameL = NULL
for (i in 1:nrow(pocket_neg_100)) {
  #print(i)
  dist_lig_diffP_diffL= c(dist_lig_diffP_diffL, dist(rbind(dt_dt[pocket_neg_100[i,"name_pock1"],features],
                                                           dt_dt[pocket_neg_100[i,"name_pock2"],features])))
}
for (i in 1:nrow(pocket_pos_100)) {
  #print(i)
  dist_lig_sameP_sameL= c(dist_lig_sameP_sameL, dist(rbind(dt_dt[pocket_pos_100[i,"name_pock1"],features],
                                                           dt_dt[pocket_pos_100[i,"name_pock2"],features])))
}
##reglog
dt_dt = as.data.frame(dt)
dt_overlap = as.data.frame(dt)
#
dt_pock_pos = sqrt((dt_overlap[as.character(pocket_pos_100[,"name_pock1"]),features] - dt_dt[as.character(pocket_pos_100[,"name_pock2"]),features])**2)
names(dt_pock_pos) = sapply(strsplit(as.character(pocket_pos_100$name_pock1), "_"), "[", 2)
dt_pock_neg = sqrt((dt_overlap[as.character(pocket_neg_100[,"name_pock1"]),features] - dt_dt[as.character(pocket_neg_100[,"name_pock2"]),features])**2)
#
dist_lig_sameP_sameL = predict.glm(dt_predict.glm.step, newdata=dt_pock_pos[,features],type = "response" )
dist_lig_diffP_diffL = predict.glm(dt_predict.glm.step, newdata=dt_pock_neg[,features],type = "response" )
#
library(sm)
sm.density.compare(c(dist_lig_diffP_diffL,
                     dist_lig_sameP_sameL),
                   c(rep(1,length(dist_lig_diffP_diffL)),
                     rep(2,length(dist_lig_sameP_sameL))),
                   model = "none",xlim = c(0,1)
                   , xlab = " distance bewteen pockets")
#
y_predict = NULL
y_predict = c(dist_lig_sameP_sameL,dist_lig_diffP_diffL)
y_true_app = NULL
y_true_app = c(rep(1,length(dist_lig_sameP_sameL)),rep(0,length(dist_lig_diffP_diffL)))
#
perf_auc(y_predict, y_true_app)


#### BAR PLOT ####
seeds = 200
seuil = 0.6
#
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
for (pocket in c(as.character(pocket_pos_100),as.character(pocket_neg_100))) {
  print(pocket)
  #
  pock2_clusters = NULL
  dt_fpocket_expanded <- dt_fpocket[rep(pocket, nrow(kmean_centroids)),features]
  rownames(dt_fpocket_expanded) = rownames(kmean_centroids)
  pock2_clusters = predict.glm(dt_predict.glm.step, newdata = sqrt((kmean_centroids[,features]-dt_fpocket_expanded[,features])**2), type = "response")
  #
  pock2_cluster = which(pock2_clusters > seuil|pock2_clusters == max(pock2_clusters))-1
  #
  dt_fpocket_expanded <- dt_fpocket[rep(pocket, nrow(dt)),features]
  rownames(dt_fpocket_expanded) = rownames(dt)
  predicted_pocket = predict.glm(dt_predict.glm.step, newdata = sqrt((dt[,features]-dt_fpocket_expanded[,features])**2), type = "response")
  names_interest = rownames(kmean_results)[which(is.element(kmean_results[,1],pock2_cluster))]
  #
  N_y_predict = c(N_y_predict,length(predicted_pocket[names_interest]))
  y_predict_30 = c(y_predict_30, sort(predicted_pocket[names_interest], decreasing = TRUE)[1:30])
  y_predict_30_mean = c(y_predict_30_mean,mean(sort(predicted_pocket[names_interest], decreasing = TRUE)[1:30]))
  N_y_predict_0.8 = c(N_y_predict_0.8, length(which(predicted_pocket[names_interest]>0.8)))
  N_y_predict_0.9 = c(N_y_predict_0.9, length(which(predicted_pocket[names_interest]>0.9)))
  N_y_predict_0.95 = c(N_y_predict_0.95, length(which(predicted_pocket[names_interest]>0.95)))
}
#save(N_y_predict_0.95, file= "../results/kmeans_results_reglog_article/N_y_predict_0.95_geom.Rdata")
#save(N_y_predict_0.9, file= "../results/kmeans_results_reglog_article/N_y_predict_0.90_geom.Rdata")
#save(N_y_predict_0.8, file= "../results/kmeans_results_reglog_article/N_y_predict_0.80_geom.Rdata")
#save(N_y_predict_0.95, file= "../results/kmeans_results_reglog_article/N_y_predict_0.95_prox.Rdata")
#save(N_y_predict_0.9, file= "../results/kmeans_results_reglog_article/N_y_predict_0.90_prox.Rdata")
#save(N_y_predict_0.8, file= "../results/kmeans_results_reglog_article/N_y_predict_0.80_prox.Rdata")
#
load("../results/kmeans_results_reglog_article/N_y_predict_0.95_geom.Rdata")
load("../results/kmeans_results_reglog_article/N_y_predict_0.90_geom.Rdata")
load("../results/kmeans_results_reglog_article/N_y_predict_0.80_geom.Rdata")
#
load("../results/kmeans_results_reglog_article/N_y_predict_0.95_prox.Rdata")
load("../results/kmeans_results_reglog_article/N_y_predict_0.90_prox.Rdata")
load("../results/kmeans_results_reglog_article/N_y_predict_0.80_prox.Rdata")
#
y_predict_30_mean
length(y_predict_30_mean)
which(N_y_predict_0.9 == 0)
library(ggplot2)
dataset = data.frame(x=c(1:200,1:200),seuil = c(N_y_predict_0.95,N_y_predict_0.9),
                     Seuil =c(rep("0.95",200),
                              rep("0.90",200)))

b <- ggplot(dataset, aes(x = x, y = seuil,fill = Seuil)) +
  geom_bar(stat = "identity") +
  ylab("Effectif de poches candidates") +
  xlab("200 poches")
b

