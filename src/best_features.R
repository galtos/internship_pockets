####select best features####
#Libraries
library(caret)
library(corrplot)
library(car)
library(FactoMineR)
library(randomForest)
####load data####
dt_12descriptors = read.table("../data/data_desc.csv", header = T, sep = ",", row.names = 1)
dt_72descriptors = read.table("../data/data_72desc.csv", header = T, sep = "\t", row.names = 1, fill=TRUE)
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
####
##CLEAN DATA
dt = dt_12descriptors[,]
dt = delete_clean_data(dt)
#Management dt_72
r_names_dt_72 = toupper(gsub("_prox5_5.desR","",rownames(dt_72descriptors)))
dt_72 = dt_72descriptors[which(is.element(r_names_dt_72, toupper(rownames(dt)))),]
rownames(dt_72) = toupper(gsub("_prox5_5.desR","",rownames(dt_72)))
summary(dt_72)
dt_72[,"hydrophobicity"] = NULL
dt_72[,"polarity"] = NULL
dt_72[,"RADIUS_MIN_CYLINDER"] = NULL
dt_72[,"drugg"] = NULL
dt_72[,"HEIGHT_MIN_CYLINDER"] = NULL
nrow(dt_72)
dt_72 = na.omit(dt_72)
nrow(dt_72)
ncol(dt_72)
####
dt = dt_72[,c("SURFACE_HULL","VOLUME_HULL","RADIUS_HULL",
              "RADIUS_CYLINDER","DIAMETER_HULL","C_ATOM",
              "p_side_chain_atom","p_O_atom","p_main_chain_atom",
              "PCI","X._ATOM_CONVEXE","p_carbone_atom","p_nitrogen_atom",
              "C_RESIDUES")]
##SCALE DT
dt = scale(dt)
##LIGAND PROT SELECTION
names_prot = sapply(strsplit(rownames(dt), "_"), "[", 1)
names_ligand = sapply(strsplit(rownames(dt), "_"), "[", 2)
unique_names_prot = unique(names_prot)
unique_names_ligand = unique(names_ligand)
table_names_ligand = table(names_ligand)
length(which(names_ligand == "HEM"))
##
n_data_lig = 2000
n_data_random = 2000
##
dt_72_bind = NULL
dt_72_random = NULL
## SAME LIGANDS
unique_names_ligand_sup1 = names(which(table_names_ligand > 1))
dt_predict = NULL

for (i in 1:n_data_lig) {
  lig = sample(unique_names_ligand_sup1, 1)
  lig = sample(unique_names_ligand_sup1, 1)
  pocket_same_lig = sample(which(names_ligand == lig),2)
  while(names_prot[pocket_same_lig][1] == names_prot[pocket_same_lig][2]) {
    lig = sample(unique_names_ligand_sup1, 1)
    pocket_same_lig = sample(which(names_ligand == lig),2)
  }
  print(lig)
  dt_72_bind = rbind(dt_72_bind,c(dt[pocket_same_lig[1],],dt[pocket_same_lig[2],]))
}

## RANDOM LIGAND
for (i in 1:n_data_random) {
  pocket_random_lig = sample(nrow(dt),2)
  dt_72_random = rbind(dt_72_random,c(dt[pocket_random_lig[1],],dt[pocket_random_lig[2],]))
}

##
nrow(dt_72_bind)
ncol(dt_72_bind)
###ADD col prediction
dt_72_bind = cbind(dt_72_bind,rep(1,nrow(dt_72_bind)))
dt_72_random = cbind(dt_72_random,rep(0,nrow(dt_72_random)))
###
dt_72_predict = rbind(dt_72_bind, dt_72_random)
dt_72_predict = as.data.frame(dt_72_predict)

colnames(dt_72_predict) = c(paste0(colnames(dt_72_bind)[1:70],"_1"),
                            paste0(colnames(dt_72_bind)[71:140],"_2"),'prediction')
colnames(dt_72_predict)
nrow(dt_72_predict)
ncol(dt_72_predict)
dt_72_predict$prediction <- as.factor(dt_72_predict$prediction)


####PREDICTION
set.seed(71)
rf <-randomForest(prediction~.,data=dt_72_predict, ntree=500) 
print(rf)

 
rf$predicted
importance(rf)

varImpPlot(rf)

###ROC
pred1=predict(rf,type = "prob")
#library(ROCR)
perf = prediction(pred1[,2], dt_72_predict$prediction)
# 1. Area under curve
auc = performance(perf, "auc")
auc
# 2. True Positive and Negative Rate
pred3 = performance(perf, "tpr","fpr")
# 3. Plot the ROC curve
plot(pred3,main="ROC Curve for Random Forest",col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")


