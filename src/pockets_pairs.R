##LOAD TOUGH-M1 DATASET
TOUGH_M1_positive = read.table("../data/TOUGH-M1/TOUGH-M1_positive.list")
TOUGH_M1_negative = read.table("../data/TOUGH-M1/TOUGH-M1_negative.list")
#
##
length(which(TOUGH_M1_positive$V5 > 0.85))
length(which(TOUGH_M1_negative$V5 > 0.85))
##Correspondance TOUGH-M1 et POCKETS names
TOUGH_M1_target = read.table("../data/TOUGH-M1/TOUGH-M1_target.list")
#
#dt_72descriptors = read.table("../data/data_PDB_72desc.txt", header = T, sep = "", row.names = 1, fill = TRUE)
#dt = dt_72descriptors
#
names_prot = sapply(strsplit(rownames(dt), "_"), "[", 1)
names_ligand = sapply(strsplit(rownames(dt), "_"), "[", 2)
names_chain = sapply(strsplit(rownames(dt), "_"), "[", 3)
names_dt_PDB_chain = paste(names_prot,names_chain,sep = "_")
#
names(names_ligand) = paste(names_prot,names_chain,sep="_")
##TARGET
apply(TOUGH_M1_target,1,paste(TOUGH_M1_target))
TOUGH_M1_PDB = paste(toupper(substr(unlist(TOUGH_M1_target),1,4)))
TOUGH_M1_PDB_chain = paste(toupper(substr(unlist(TOUGH_M1_target),1,4)),
                           toupper(substr(unlist(TOUGH_M1_target),5,5)), sep="_")
##POSITIVE
apply(TOUGH_M1_positive$V1,1,paste(TOUGH_M1_positive$V1))
TOUGH_M1_PDB = paste(toupper(substr(unlist(TOUGH_M1_positive$V1),1,4)))
TOUGH_M1_PDB_chain = paste(toupper(substr(unlist(TOUGH_M1_positive$V1),1,4)),
                           toupper(substr(unlist(TOUGH_M1_positive$V1),5,5)), sep="_")
##PDB
length(intersect(names_prot,TOUGH_M1_PDB))
setdiff(TOUGH_M1_PDB,names_prot)
##PDB+chain
length(TOUGH_M1_PDB_chain)
length(intersect(names_dt_PDB_chain,TOUGH_M1_PDB_chain))
setdiff(TOUGH_M1_PDB_chain,names_dt_PDB_chain)
##
names_ligand[setequal(names_prot,TOUGH_M1_PDB)]
table_names_ligand = table(names_ligand[TOUGH_M1_PDB_chain])
head(sort(table_names_ligand, decreasing=T))
sum(table_names_ligand)

#
length(table_names_ligand)
barplot(table_names_ligand)
#
length(which(is.element(names_prot,TOUGH_M1_PDB)==TRUE))
#
############ Get exemple lig for TANIMOTO ##############
## data après nettoyage
rnames_dt = rownames(dt)

length(rnames_dt)

length(unique(names_ligand))

names_pocket_ligand_unique = NULL
for (lig in unique(names_ligand)) {
  print(lig)
  pocket_random = sample(rnames_dt[which(names_ligand == lig)],1)
  names_pocket_ligand_unique = c(names_pocket_ligand_unique,pocket_random)
}
names_prot = sapply(strsplit(names_pocket_ligand_unique, "_"), "[", 1)
names_ligand = sapply(strsplit(names_pocket_ligand_unique, "_"), "[", 2)
names_pocket_ligand_unique = paste(names_prot,names_ligand,sep=" ")
#writeLines(names_pocket_ligand_unique, "../data/names_pocket_ligand_unique.txt")
output.file <- file("../data/list_names_extract_sdf_PDB.txt", "wb")
writeLines(names_pocket_ligand_unique,output.file)
close(output.file)

write(paste(unique(names_ligand), collapse = ','), "../data/list_lig_unique.txt")

t = table(names_ligand)
which(t>1)
#TEST TANIMOTO
library(vegan)
data(varespec)
vare.dist <- vegdist(varespec, method = "jaccard")

#### SDF to MACCS ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rcpi")
BiocManager::install("ChemmineOB")

library(Rcpi)
library(ChemmineOB)
mol = readChar(system.file('../4A3H_DCB_A_1.sdf', package = 'Rcpi'),
                nchars = 1e+6)
mol = 'C1CCC1CC(CN(C)(C))CC(=O)CC'
smif = extractDrugOBMACCS(mol, type = 'smile')
## FAMILY LIG 190 regad 
lig_names = read.table("../../results/lig_classification/lig_names_regad.txt")
lig_names[,1]
length(names_ligand) - length(which(is.element(names_ligand,lig_names[,1]) == FALSE))
setdiff(names_ligand,lig_names[,1])
length(names_prot)
prot_exclude_lig = names_prot[which(is.element(names_ligand,lig_names[,1]) == FALSE)]

### H CLUST matrice ###
mat_dist = read.csv("../../results/lig_classification/mat_dist_tanimoto_regad.csv",row.names = 1, header = T)

hc_ward = hclust(as.dist(1-mat_dist), method = "ward.D2")
hc_average = hclust(as.dist(1-mat_dist), method = "average")
SEUIL = 0.1 #(1-0.9)
groups = cutree(hc_average, h=SEUIL)
                
pdf("../../results/lig_classification/hclust_dist_tanimoto_m_w.pdf")
plot(hc_ward, hang = -1, cex = 0.6)
dev.off()

pdf("../../results/lig_classification/hclust_dist_tanimoto_m_average.pdf")
plot(hc_average, hang = -1, cex = 0.6)
dev.off()

pdf("../../results/lig_classification/hclust_dist_tanimoto_m_average_table_ctree.pdf")
plot(table(cutree(hc, h=SEUIL)))
dev.off()

###Combien de mL par rapport à ligand ###
#mat_sL = data.frame(names_pockets = rownames(dt), lig = groups[names_ligand])
#mat_sL = write.table(mat_sL, file = "../data/lig_classification/mat_pocket_lig.txt",sep =",", quote = F, row.names=F)
mat_sL = read.table("../data/lig_classification/mat_pocket_lig.txt",sep =",",header = T, row.names = 1)

plot(table(mat_sL$lig))

length(which(table(mat_sL$lig) >= 2))
nrow(mat_sL[which(mat_sL$lig == 23),])
####
class_lig = read.table("../results/lig_classification/lig_class_0.9.txt",sep =",",header = T,row.names = 1)
class_lig["17X",1]
summary(class_lig)
groups = class_lig[,1]
names(groups) = rownames(class_lig)

groups["17X"]

mat_sL = data.frame(names_pockets = rownames(dt), lig = groups[names_ligand])

list_na = which(is.na(mat_sL$lig))
list_na
list_class_na = as.vector(Na_groups[sapply(strsplit(as.character(mat_sL[list_na,1]), "_"), "[", 2)])
list_class_na
length(list_class_na)
length(list_na)
mat_sL[list_na,2] = list_class_na

summary(mat_sL$lig)
Na_groups = unique(sapply(strsplit(as.character(mat_sL[is.na(mat_sL$lig),1]), "_"), "[", 2))
Na_groups = (max(groups)+1):(max(groups)+length(Na_groups))
names(Na_groups) = unique(sapply(strsplit(as.character(mat_sL[is.na(mat_sL$lig),1]), "_"), "[", 2))

groups["HEM"]
##
mat_sL[which(mat_sL$lig == 0),]

summary(mat_sL$lig)
barplot(table(mat_sL$lig))
table(mat_sL$lig)[1]
#plot distribution par categories
barplot(table(mat_sL$lig)[which(table(mat_sL$lig)==1)])

barplot(table(mat_sL$lig)[which(table(mat_sL$lig)>100)])

df <- data.frame(
  group = c("1","]1;5]",  "]5;10]", "]10;1000]", "> 1000"),
  value = c(sum(table(mat_sL$lig)[which(table(mat_sL$lig) == 1)]),
            sum(table(mat_sL$lig)[which(table(mat_sL$lig) <= 5 & table(mat_sL$lig) > 1)]),
            sum(table(mat_sL$lig)[which(table(mat_sL$lig) <= 10 & table(mat_sL$lig) > 5)]),
            sum(table(mat_sL$lig)[which(table(mat_sL$lig) <= 1000 & table(mat_sL$lig) > 10)]),
            sum(table(mat_sL$lig)[which(table(mat_sL$lig) > 1000)])
  ),
  size = c(length(table(mat_sL$lig)[which(table(mat_sL$lig) == 1)]),
           length(table(mat_sL$lig)[which(table(mat_sL$lig) <= 5 & table(mat_sL$lig) > 1)]),
           length(table(mat_sL$lig)[which(table(mat_sL$lig) <= 10 & table(mat_sL$lig) > 5)]),
           length(table(mat_sL$lig)[which(table(mat_sL$lig) <= 1000 & table(mat_sL$lig) > 10)]),
           length(table(mat_sL$lig)[which(table(mat_sL$lig) > 1000)])
  )
)
library(ggplot2)
# Bar plot
bp<- ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1,stat="identity")+
  geom_text(aes(label=size))

bp + theme(legend.text = element_text(colour="black", size=5, 
                                      face="bold"))


#
ggplot(data=df, aes(x=group, y=value)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=size), vjust=1.6, color="black", size=3.5)+
  xlab("categorie de ligands liés aux poches")+ylab("Nombre de poches associés")
  theme_minimal()

##### Creation des paires de poches ####
#Pckets with overlaping score > Seuil
Seuil_so = 0.1
SO_fpocket = read.table("../data/SO_desc_fpocket_all.txt",sep=",",header = F)
summary(SO_fpocket)
SO_fpocket$V2 = as.numeric(levels(SO_fpocket$V2))[SO_fpocket$V2]

length(SO_fpocket[which(SO_fpocket$V2 > Seuil_so),1])
so_pocket_names = as.character( SO_fpocket[which(SO_fpocket$V2 > Seuil_so),1])
length(so_pocket_names)
# SO VS fpocket
dt_dt = dt
dt_overlap = dt
#
dt_overlap = na.omit(dt_overlap[,features])
#
y_predict = sqrt((dt_dt[rownames(dt_overlap),features] - dt_overlap[rownames(dt_overlap),features])**2)
names(y_predict) = rownames(dt_overlap)
#
SO_50 = SO_fpocket$V2
names(SO_50) = as.character(SO_fpocket$V1)
#
sLsP_50 = cbind(SO_50[rownames(dt_overlap)],y_predict[rownames(dt_overlap)])
sLsP_50 = as.data.frame(sLsP_50)
colnames(sLsP_50) = c("SO","ypred")
cor.test(sLsP_50$SO,sLsP_50$ypred)
#
plot(sLsP_50$SO[which(sLsP_50$ypred>1)],sLsP_50$ypred[which(sLsP_50$ypred>1)], xlab = "score overlap", ylab = "distance", ylim=c(0,6))
cor.test(sLsP_50$SO[which(sLsP_50$ypred>1)],sLsP_50$ypred[which(sLsP_50$ypred>1)])
#
boxplot(sLsP_50$SO)
mean(sLsP_50$SO)
sd(sLsP_50$SO)
hist(sLsP_50$SO)
#
boxplot(sLsP_50$ypred)
mean(sLsP_50$ypred)
sd(sLsP_50$ypred)
hist(sLsP_50$ypred)
###
library(ggplot2)
library(ggExtra)
g <- ggplot(sLsP_50, aes(SO, ypred)) +
  geom_count() +
  geom_smooth(method="lm", se=F) + xlim(0,0.6)+ylim(0,1)


ggMarginal(g, type = "histogram", fill="transparent")

p1 <- ggplot(sP_50, aes(SO,ypred)) +
  geom_smooth(method="lm", se=F) + xlim(0,0.8)+ylim(0,1)
ggMarginal(p1 + geom_point(), type = "histogram", fill="transparent")

# POCHES A RETIRER POUR TEST
dt[which(names_prot  == "1HOP"),1:2]
names_prot_so = sapply(strsplit(so_pocket_names, "_"), "[", 1)

length(which(names_prot_so == "1HOP"))
so_pocket_names = so_pocket_names[which(names_prot_so != "1HOP")]
#DIST POCK same prot same lig
dist_lig_sameP_sameL = function(dt, names_ligand, mat_sL, so_pocket_names) {
  dist_ligs = rep(0,length(which(table(names_ligand) > 1)))
  names(dist_ligs) = names(which(table(names_ligand) > 1))
  names_ligand_unique = unique(mat_sL[so_pocket_names,1])
  flag = 1
  dist_features_pock_all = NULL
  dist_features = matrix(data = 0,
                         nrow = length(which(table(names_ligand) > 1)),
                         ncol = ncol(dt))
  dist_features = as.data.frame(dist_features)
  rownames(dist_features) = names(which(table(names_ligand) > 1))
  data_samePsameL = NULL
  rnames = so_pocket_names
  for (lig_name in names_ligand_unique){#length(names_ligand_unique)
    #print(i)
    index = which(mat_sL[so_pocket_names,1] == lig_name)
    if(length(index) > 1) {
      print(flag)
      if(length(index) > 100) {
        index = sample(index,100)
      }
      dist_pock = NULL
      dist_features_pock = NULL
      for (j in 1:(length(index)-1)) {
        for (k in (j+1):length(index)) {
          if(sapply(strsplit(rnames[index[j]], "_"), "[", 1) != sapply(strsplit(rnames[index[k]], "_"), "[", 1)) { #CHANGE TO != FOR DIFFERENT PROTEINS
            data_samePsameL = rbind(data_samePsameL, c(rnames[index[j]],rnames[index[k]]))
            #dist_pock = c(dist_pock, dist(rbind(dt[index[j],],
            #                                    dt[index[k],])))
            #dist_features_pock = rbind(dist_features_pock,  apply(dt[c(index[j],index[k]),],2,dist))
          }
        }
      }
      flag = flag+1
      if(!is.null(dist_pock)) {
        dist_ligs[lig_name] = mean(dist_pock)
        #dist_features_pock_all = rbind(dist_features_pock_all,dist_features_pock)
        #dist_features[lig_name,] = apply(dist_features_pock,2,mean)
        ##flag = flag+1
      }
    }
  }
  data_samePsameL = as.data.frame(data_samePsameL)
  colnames(data_samePsameL) = c("name_pock1","name_pock2")
  #write.csv(data_samePsameL, "../data/data_structure_comparaison/data_samePsameL_dt72_pharmacophores.csv")
  #save(dist_features_pock_all, file = "../results/protocol_features/dist_features_sameP_sameL_pock.Rdata")
  #save(dist_features, file = "../results/protocol_features/dist_features_sameP_sameL.Rdata")
  return(data_samePsameL)
}
data_samePsameL = dist_lig_sameP_sameL(dt,names_ligand, mat_sL, so_pocket_names)
write.csv(data_samePsameL, "../data/data_structure_comparaison/data_samePsameL_dt72_pharmacophores_positive.csv")

nrow(data_samePsameL)

#negative pairs
dist_lig_diffP_diffL = function(dt, names_ligand,mat_sL,so_pocket_names, n = 50000) {
  unique_names_ligand = unique(mat_sL[so_pocket_names,1])
  dist_ligs_random = rep(0,n)#length(unique(names_ligand)))
  dist_features = NULL
  dist_features_pock_all = NULL
  data_diffPdiffL = NULL
  rnames = so_pocket_names
  for (i in 1:n){#length(unique(names_ligand))){
    print(i)
    lig = sample(unique_names_ligand, 2)
    i_pock1 = sample(which(mat_sL[so_pocket_names,1]==lig[1]),1)
    i_pock2 = sample(which(mat_sL[so_pocket_names,1]==lig[2]),1)
    #dist_ligs_random[i] = dist(rbind(dt[i_pock1,],
    #                                 dt[i_pock2,]))
    data_diffPdiffL = rbind(data_diffPdiffL, c(rnames[i_pock1],rnames[i_pock2]))
    #dist_features = rbind(dist_features,  apply(dt[c(i_pock1,i_pock2),],2,dist))
    
  }
  data_diffPdiffL = as.data.frame(data_diffPdiffL)
  colnames(data_diffPdiffL) = c("name_pock1","name_pock2")
  #write.csv(data_diffPdiffL, "../data/data_structure_comparaison/data_diffPdiffL_dt72_pharmacophores_50000.csv")
  #dist_features_pock_all = rbind(dist_features_pock_all,dist_features)
  #save(dist_features_pock_all, file = "../results/protocol_features/dist_features_diffP_diffL_pock.Rdata")
  #save(dist_features, file = "../results/protocol_features/dist_features_diffP_diffL.Rdata")
  return(data_diffPdiffL)
}
data_diffPdiffL = dist_lig_diffP_diffL(dt,names_ligand, mat_sL, so_pocket_names)
write.csv(data_diffPdiffL, "../data/data_structure_comparaison/data_diffPdiffL_dt72_pharmacophores_negative.csv")
# to linux format


###PLOT DIST ###
data_samePsameL_positive = read.csv("../data/data_structure_comparaison/data_samePsameL_dt72_pharmacophores_positive.csv", colClasses = "character")
data_diffPdiffL_negative = read.csv("../data/data_structure_comparaison/data_diffPdiffL_dt72_pharmacophores_negative.csv", colClasses = "character")
#
output.f = file("../data/data_structure_comparaison/pairs_positive.txt","wb")
write.table(data_samePsameL_positive,file = output.f,row.names = F,col.names = F,quote = F,append = T,sep=",")
close(output.f)
output.f = file("../data/data_structure_comparaison/pairs_negative.txt","wb")
write.table(data_diffPdiffL_negative,file = output.f,row.names = F,col.names = F,quote = F,append = T,sep=",")
close(output.f)
##samePsameL
dist_lig_sameP_sameL = NULL
for (i in 1:nrow(data_samePsameL_positive)) {
  print(i)
  dist_lig_sameP_sameL = c(dist_lig_sameP_sameL,
                           dist(rbind(dt[data_samePsameL_positive[i,"name_pock1"],features],
                                      dt[data_samePsameL_positive[i,"name_pock2"],features])))
}
##diffPdiffL
dist_lig_diffP_diffL = NULL
for (i in 1:nrow(data_diffPdiffL_negative)) {
  print(i)
  dist_lig_diffP_diffL = c(dist_lig_diffP_diffL,
                           dist(rbind(dt[data_diffPdiffL_negative[i,"name_pock1"],features],
                                      dt[data_diffPdiffL_negative[i,"name_pock2"],features])))
}
#library(sm)
sm.density.compare(c(dist_lig_diffP_diffL,
                     dist_lig_sameP_sameL),
                   c(rep(1,length(dist_lig_diffP_diffL)),
                     rep(2,length(dist_lig_sameP_sameL))),
                   model = "none",xlim = c(0,100)
                   , xlab = " distance bewteen pockets")

## COMPARE TM SCORE
###MAT positive prox
mat_tmscsore = read.table("../data/tmscore_pock_PDB/tmscore_positive.txt",sep=",",header = F)
mat_tmscsore_geom = read.table("../data/tmscore_pock_PDB/tmscore_positive_geom.txt",sep=",",header = F)
#
mat_tmscsore$V1 = as.character(mat_tmscsore$V1)
mat_tmscsore$V2 = as.character(mat_tmscsore$V2)
mat_tmscsore_geom$V1 = as.character(mat_tmscsore_geom$V1)
mat_tmscsore_geom$V2 = as.character(mat_tmscsore_geom$V2)
#
mat_res = data_samePsameL_positive
mat_res = cbind(mat_res, rep(NA,nrow(data_samePsameL_positive)))
mat_res = cbind(mat_res, rep(NA,nrow(data_samePsameL_positive)))
mat_res = cbind(mat_res, rep(NA,nrow(data_samePsameL_positive)))
mat_res = cbind(mat_res, rep(NA,nrow(data_samePsameL_positive)))
colnames(mat_res)
c_names = colnames(mat_res)
c_names[4:7] = c("dist_prox","tmscore_prox","dist_geom","tmscore_geom")
colnames(mat_res) = c_names
#
#dt_dt = dt
#dt_o = dt
#
tmscore_lig_sameP_sameL = NULL
dist_lig_sameP_sameL = NULL
for (i in 1:nrow(data_samePsameL_positive)) {
  #print(i)
  mat_res[i,"dist_prox"]= dist(rbind(dt_dt[data_samePsameL_positive[i,"name_pock1"],features],
                                     dt_dt[data_samePsameL_positive[i,"name_pock2"],features]))
}
print("done")
for (i in 1:nrow(mat_tmscsore)) {
  #print(i)
  N = which(mat_tmscsore[i,1] == mat_res$name_pock1 & mat_tmscsore[i,2] == mat_res$name_pock2)
  mat_res[N,"tmscore_prox"] = mat_tmscsore[i,3]
}
print("done")
for (i in 1:nrow(data_samePsameL_positive)) {
  #print(i)
  mat_res[i,"dist_geom"]= dist(rbind(dt_o[data_samePsameL_positive[i,"name_pock1"],features],
                                     dt_o[data_samePsameL_positive[i,"name_pock2"],features]))
}
print("done")
for (i in 1:nrow(mat_tmscsore_geom)) {
  #print(i)
  N = which(mat_tmscsore_geom[i,1] == mat_res$name_pock1 & mat_tmscsore_geom[i,2] == mat_res$name_pock2)
  mat_res[N,"tmscore_geom"] = mat_tmscsore_geom[i,3]
}
#
plot(mat_res$dist,mat_res$tmscore)
cor.test(mat_res$dist,mat_res$tmscore)
#mat_res = mat_res[,2:5]
output.f = file("../data/data_structure_comparaison/mat_pairs_infos_positive.txt","wb")
write.table(mat_res,file = output.f,row.names = F,col.names = T,quote = F,append = T,sep=",")
close(output.f)

#
print("positive done")
###MAT negative prox
 
sm.density.compare(c(mat_infos_neg$tmscore_prox,
                     mat_infos_pos$tmscore_prox),
                   c(rep(1,length(which(!is.na(mat_infos_neg$tmscore_prox)))),
                     rep(2,length(which(!is.na(mat_infos_pos$tmscore_prox))))),
                   model = "none",xlim = c(0,1)
                   , xlab = " distance bewteen pockets")
                   
###### JEU TEST DE 100 POCHES
# paires positives
mat_infos_pos = read.table("../data/data_structure_comparaison/mat_pairs_infos_positive.txt",sep=",",header=T)
# paires negatives
mat_infos_neg = read.table("../data/data_structure_comparaison/mat_pairs_infos_negative.txt",sep=",",header=T)
#### RETIRER 100 poches
#Retirer 4opa et 1qjc
#100
set.seed(83)
pocket_pos_100 = mat_infos_pos[sample(nrow(mat_infos_pos),size = 100),2]
length(pocket_pos_100)
pocket_neg_100 = mat_infos_neg[sample(nrow(mat_infos_neg),size = 100),2]
length(pocket_neg_100)
#
length(unique(as.character(pocket_pos_100)))
length(unique(as.character(pocket_neg_100)))
#save(pocket_pos_100, file = "../data/data_structure_comparaison/pocket_pos_100.Rdata")
#save(pocket_neg_100, file = "../data/data_structure_comparaison/pocket_neg_100.Rdata")
#load(file = "../data/data_structure_comparaison/pocket_pos_100.Rdata")
#load(file = "../data/data_structure_comparaison/pocket_neg_100.Rdata")

##mat_infos_pos mat_infos_neg - 100
#pos
prot_pos_100 = sapply(strsplit(as.character(pocket_pos_100), "_"), "[", 1)
length(prot_pos_100)
length(unique(prot_pos_100))
#neg
prot_neg_100 = sapply(strsplit(as.character(pocket_neg_100), "_"), "[", 1)
length(prot_neg_100)
length(unique(prot_neg_100))
#
### retirer matrice
#pos
name_mat_pos_prot_1 = sapply(strsplit(as.character(mat_infos_pos$name_pock1), "_"), "[", 1)
name_mat_pos_prot_2 = sapply(strsplit(as.character(mat_infos_pos$name_pock2), "_"), "[", 1)
#neg
name_mat_neg_prot_1 = sapply(strsplit(as.character(mat_infos_neg$name_pock1), "_"), "[", 1)
name_mat_neg_prot_2 = sapply(strsplit(as.character(mat_infos_neg$name_pock2), "_"), "[", 1)
#
is.element("1H1T",c(name_mat_pos_prot_1,name_mat_pos_prot_2,name_mat_neg_prot_1,name_mat_neg_prot_2))
which("1QJC" == names_prot)
dt[10794,1:3]
rownames(dt)[10794]

is.element("2XCS",c(name_mat_pos_prot_1,name_mat_pos_prot_2,name_mat_neg_prot_1,name_mat_neg_prot_2))
which("2XCS" == names_prot)
dt[29915,1:3]
rownames(dt)[29915]
#p
index_p_100 = which(is.element(name_mat_pos_prot_1,prot_pos_100) | is.element(name_mat_pos_prot_2,prot_pos_100) |
                    is.element(name_mat_pos_prot_1,prot_neg_100) | is.element(name_mat_pos_prot_2,prot_neg_100) |
                    is.element(name_mat_pos_prot_1,"1QJC") | is.element(name_mat_pos_prot_2,"1QJC") |
                    is.element(name_mat_pos_prot_1,"2XCS") | is.element(name_mat_pos_prot_2,"2XCS"))
mat_res = mat_infos_pos[-index_p_100,]
output.f = file("../data/data_structure_comparaison/mat_pairs_infos_positive_100.txt","wb")
write.table(mat_res,file = output.f,row.names = F,col.names = T,quote = F,append = T,sep=",")
close(output.f)
#n
index_n_100 = which(is.element(name_mat_neg_prot_1,prot_pos_100) | is.element(name_mat_neg_prot_2,prot_pos_100) |
                    is.element(name_mat_neg_prot_1,prot_neg_100) | is.element(name_mat_neg_prot_2,prot_neg_100) |
                    is.element(name_mat_neg_prot_1,"1QJC") | is.element(name_mat_neg_prot_2,"1QJC") |
                    is.element(name_mat_neg_prot_1,"2XCS") | is.element(name_mat_neg_prot_2,"2XCS"))
mat_res = mat_infos_neg[-index_n_100,]
output.f = file("../data/data_structure_comparaison/mat_pairs_infos_negative_100.txt","wb")
write.table(mat_res,file = output.f,row.names = F,col.names = T,quote = F,append = T,sep=",")
close(output.f)
##CSV
name_prot_pos_100 = rownames(dt)[which(is.element(names_prot,prot_pos_100))]
name_prot_neg_100 = rownames(dt)[which(is.element(names_prot,prot_neg_100))]
write.csv(dt[setdiff(rownames(dt),c(as.character(pocket_pos_100),as.character(pocket_neg_100),"1QJC_PNS_B_1","2XTO_GDP_B_2")),],
          file = "../data/dt_72clean-100.csv")


### REGRESSION LOGISTIQUE ###
features = c(descriptors_hydrophobicity,
             descriptors_aromatic,
             descriptors_polarity,
             descriptors_physicochemical,
             descriptors_geometrical)
features = c(features,"A","C","E","D","G","F","I","H","K","M","L","N","Q","P","S","R","T","W","V","Y")
# DATA
# paires positives
mat_infos_pos = read.table("../data/data_structure_comparaison/mat_pairs_infos_positive_100.txt",sep=",",header=T)
# paires negatives
mat_infos_neg = read.table("../data/data_structure_comparaison/mat_pairs_infos_negative_100.txt",sep=",",header=T)
## PROXIMITE
dt_pock_pos = sqrt((dt[as.character(mat_infos_pos[,"name_pock1"]),features] - dt[as.character(mat_infos_pos[,"name_pock2"]),features])**2)
names(dt_pock_pos) = sapply(strsplit(as.character(mat_infos_pos$name_pock1), "_"), "[", 2)
dt_pock_neg = sqrt((dt[as.character(mat_infos_neg[,"name_pock1"]),features] - dt[as.character(mat_infos_neg[,"name_pock2"]),features])**2)
#
#FOR OVERLAP
#withoutNA
try(na.fail(dt[,features]))
summary(dt[,features])
is.element(names(which(is.na(dt[,"RADIUS_CYLINDER"]))),as.character(mat_infos_pos[,"name_pock1"]))
is.element(names(which(is.na(dt[,"RADIUS_CYLINDER"]))),as.character(mat_infos_pos[,"name_pock2"]))
is.element(names(which(is.na(dt[,"RADIUS_CYLINDER"]))),as.character(mat_infos_neg[,"name_pock1"]))
is.element(names(which(is.na(dt[,"RADIUS_CYLINDER"]))),as.character(mat_infos_neg[,"name_pock2"]))
#
which(is.na(dt[,"hydrophobic_kyte"]))
is.element(names(which(is.na(dt[,"hydrophobic_kyte"]))),as.character(mat_infos_pos[,"name_pock1"]))
is.element(names(which(is.na(dt[,"hydrophobic_kyte"]))),as.character(mat_infos_pos[,"name_pock2"]))
is.element(names(which(is.na(dt[,"hydrophobic_kyte"]))),as.character(mat_infos_neg[,"name_pock1"]))
is.element(names(which(is.na(dt[,"hydrophobic_kyte"]))),as.character(mat_infos_neg[,"name_pock2"]))

is.element("2XCS_RXV_F_1",as.character(mat_infos_neg[,"name_pock1"]))
#
dt = na.omit(dt[,features])
na.action(dt)
#
index_mat_infos_pos = which(is.element(mat_infos_pos[,"name_pock1"],rownames(dt)) == TRUE &
                            is.element(mat_infos_pos[,"name_pock2"],rownames(dt)) == TRUE)
#
index_mat_infos_neg = which(is.element(mat_infos_neg[,"name_pock1"],rownames(dt)) == TRUE &
                            is.element(mat_infos_neg[,"name_pock2"],rownames(dt)) == TRUE)

dt_pock_pos = sqrt((dt[mat_infos_pos[index_mat_infos_pos,"name_pock1"],features] - dt[mat_infos_pos[index_mat_infos_pos,"name_pock2"],features])**2)
names(dt_pock_pos) = sapply(strsplit(as.character(mat_infos_pos[index_mat_infos_pos,]$name_pock1), "_"), "[", 2)
dt_pock_neg = sqrt((dt[mat_infos_neg[index_mat_infos_neg,"name_pock1"],features] - dt[mat_infos_neg[index_mat_infos_neg,"name_pock2"],features])**2)

#data preparation
dt_predict = NULL
dt_predict = rbind(dt_pock_pos,dt_pock_neg)
y_true = c(rep(1,nrow(dt_pock_pos)),
           rep(0,nrow(dt_pock_neg))) 
dt_predict = cbind(abs(dt_predict),y_true)
dt_predict = as.data.frame(dt_predict)
dt_predict$y_true = as.factor(dt_predict$y_true)
### data preparation ###
set.seed(83)
Index = sample(1:nrow(dt_predict),size = nrow(dt_predict)*2/3 )
dt_predict_app = dt_predict[Index,]
dt_predict_test = dt_predict[-Index,]

y_true_app = y_true[Index]
y_true_test = y_true[-Index]
#### Linear regression#### 
dt_predict.glm = glm(y_true~.,data = dt_predict_app, family = "binomial")
attributes(dt_predict.glm)
summary(dt_predict.glm)

dt_predict.glm.step = step(dt_predict.glm, direction = "both")
#save(dt_predict.glm.step, file = "../results/kmeans_results_reglog/model_article_prox.glm.step.Rdata")
#save(dt_predict.glm.step, file = "../results/kmeans_results_reglog/model_article_geom.glm.step.Rdata")

#load(file = "../results/kmeans_results_reglog/model_article_prox.glm.step.Rdata")
#load(file = "../results/kmeans_results_reglog/model_article_geom.glm.step.Rdata")
summary(dt_predict.glm.step)

### PERFORMANCES ###
library(ROCR)
perf_auc = function(y_predict, y_true_test) {
  dt.pred = prediction(y_predict, y_true_test) #y_predict[-Index]#1-y_predict#y_true_test#y_predict[,2]
  dt.auc = performance(dt.pred, "auc")
  print("AUC")
  print(attr(dt.auc, "y.values"))
}
#
y_predict = predict.glm(dt_predict.glm, newdata=dt_predict_app[,features],type = "response" )
y_predict = predict.glm(dt_predict.glm, newdata=dt_predict_test[,features],type = "response" )
#
perf_auc(y_predict, y_true_app)
perf_auc(y_predict, y_true_test)

## perf tm score
#dist
m = rbind(mat_infos_pos, mat_infos_neg)
perf_auc(m[Index,"dist_prox"],y_true_app)
perf_auc(m[-Index,"dist_prox"],y_true_test)

dist_app = m[Index,"dist_geom"]
dist_app_na = dist_app[-which(is.na(dist_app))]
y_true_app_na = y_true_app[-which(is.na(dist_app))]

perf_auc(dist_app_na,y_true_app_na)

dist_test = m[-Index,"dist_geom"]
dist_test_na = dist_test[-which(is.na(dist_test))]
y_true_test_na = y_true_test[-which(is.na(dist_test))]
perf_auc(dist_test_na,y_true_test_na)

dist = 
#tmscore
nrow(m)
nrow(dt_predict)

tmscore_app = m[Index,"tmscore_prox"]
tmscore_app_na = tmscore_app[-which(is.na(tmscore_app))]
y_true_app_na = y_true_app[-which(is.na(tmscore_app))]

perf_auc(tmscore_app_na,y_true_app_na)

tmscore_test = m[-Index,"tmscore_prox"]
tmscore_test_na = tmscore_test[-which(is.na(tmscore_test))]
y_true_test_na = y_true_test[-which(is.na(tmscore_test))]

perf_auc(tmscore_test_na,y_true_test_na)

length(y_true_test_na)
#
nrow(m)
nrow(dt_predict)

tmscore_app = m[Index,"tmscore_geom"]
tmscore_app_na = tmscore_app[-which(is.na(tmscore_app))]
y_true_app_na = y_true_app[-which(is.na(tmscore_app))]

perf_auc(tmscore_app_na,y_true_app_na)

tmscore_test = m[-Index,"tmscore_geom"]
tmscore_test_na = tmscore_test[-which(is.na(tmscore_test))]
y_true_test_na = y_true_test[-which(is.na(tmscore_test))]

perf_auc(tmscore_test_na,y_true_test_na)

length(y_true_test)
length(tmscore_test)

### PLOT DENSITY ###
dist_lig_diffP_diffL = mat_infos_neg$tmscore_prox
dist_lig_sameP_sameL = mat_infos_pos$tmscore_prox
#
dist_lig_diffP_diffL = mat_infos_neg$tmscore_geom
dist_lig_sameP_sameL = mat_infos_pos$tmscore_geom
#
dist_lig_diffP_diffL = mat_infos_neg$dist_prox
dist_lig_sameP_sameL = mat_infos_pos$dist_prox
#
dist_lig_diffP_diffL = mat_infos_neg$dist_geom
dist_lig_sameP_sameL = mat_infos_pos$dist_geom
#
dist_lig_diffP_diffL = predict.glm(dt_predict.glm.step, newdata=dt_predict[which(dt_predict$y_true == 0),features],type = "response" )
dist_lig_sameP_sameL = predict.glm(dt_predict.glm.step, newdata=dt_predict[which(dt_predict$y_true == 1),features],type = "response" )
#
y_predict = c(na.omit(dist_lig_sameP_sameL),na.omit(dist_lig_diffP_diffL))
y_true = c(rep(1,length(na.omit(dist_lig_sameP_sameL))),rep(0,length(na.omit(dist_lig_diffP_diffL))))
perf_auc(y_predict,y_true)
#
library(sm)
sm.density.compare(c(dist_lig_diffP_diffL,
                     dist_lig_sameP_sameL),
                   c(rep(1,length(dist_lig_diffP_diffL)),
                     rep(2,length(dist_lig_sameP_sameL))),
                   model = "none",xlim = c(0,max(na.omit(c(dist_lig_diffP_diffL,dist_lig_sameP_sameL))))
                   , xlab = " distance bewteen pockets")
### KMEANS RESULT ###

