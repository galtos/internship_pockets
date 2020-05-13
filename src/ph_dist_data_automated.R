##AUTOMATED PHARMACOPHORES ANALYSIS ###
print("START")
print("Rcpp")
library("Rcpp")
sourceCpp("C_code/dist_fuzcav.cpp")
sourceCpp("C_code/dist_fuzcav_norm.cpp")
print("READ  PHARMACOPHORE FILE")
dt_pharmacophores = read.table("../data/FPCount_save_all_inter_dt72.txt",sep
                               =";",row.names=1)
print("to data frame")
dt = as.data.frame(dt_pharmacophores)
print("Get lig prot names")
names_ligand = sapply(strsplit(rownames(dt),"_"),"[",2)
names_prot = sapply(strsplit(rownames(dt),"_"),"[",1)
#
rnames_dt = rownames(dt)
names_ligand_unique = unique(names_ligand)
#
data_samePsameL_dt72_pharmacophores = read.csv("../data/data_structure_comparaison/data_samePsameL_dt72_pharmacophores.csv", colClasses = "character")
data_diffPsameL_dt72_pharmacophores = read.csv("../data/data_structure_comparaison/data_diffPsameL_dt72_pharmacophores.csv", colClasses = "character")
data_diffPdiffL_dt72_pharmacophores = read.csv("../data/data_structure_comparaison/data_diffPdiffL_dt72_pharmacophores.csv", colClasses = "character")
#
print("###Distance enclidean ###")
##samePsameL
dist_lig_sameP_sameL = NULL
for (i in 1:nrow(data_samePsameL_dt72_pharmacophores)) {
  dist_lig_sameP_sameL = c(dist_lig_sameP_sameL,
                           dist(rbind(dt[data_samePsameL_dt72_pharmacophores[i,"name_pock1"],],
                                      dt[data_samePsameL_dt72_pharmacophores[i,"name_pock2"],])))
}
save(dist_lig_sameP_sameL, file = "../results/pharmacophores_results/dist_data_intersect72/dist_lig_sameP_sameL_euclidean.Rdata")
print("dist_lig_sameP_sameL dist euclidean")
##diffPsameL
dist_lig_diffP_sameL = NULL
for (i in 1:nrow(data_diffPsameL_dt72_pharmacophores)) {
  dist_lig_diffP_sameL = c(dist_lig_diffP_sameL,
                           dist(rbind(dt[data_diffPsameL_dt72_pharmacophores[i,"name_pock1"],],
                                      dt[data_diffPsameL_dt72_pharmacophores[i,"name_pock2"],])))
}
save(dist_lig_diffP_sameL, file = "../results/pharmacophores_results/dist_data_intersect72/dist_lig_diffP_sameL_euclidean.Rdata")
print("dist_lig_diffP_sameL dist euclidean")
##diffPdiffL
dist_lig_diffP_diffL = NULL
for (i in 1:nrow(data_diffPdiffL_dt72_pharmacophores)) {
  dist_lig_diffP_diffL = c(dist_lig_diffP_diffL,
                           dist(rbind(dt[data_diffPdiffL_dt72_pharmacophores[i,"name_pock1"],],
                                      dt[data_diffPdiffL_dt72_pharmacophores[i,"name_pock2"],])))
}
save(dist_lig_diffP_diffL, file = "../results/pharmacophores_results/dist_data_intersect72/dist_lig_diffP_diffL_euclidean.Rdata")
print("dist_lig_diffP_diffL dist euclidean")
#################################
print("###Distance manhattan ###")
##samePsameL
dist_lig_sameP_sameL = NULL
for (i in 1:nrow(data_samePsameL_dt72_pharmacophores)) {
  dist_lig_sameP_sameL = c(dist_lig_sameP_sameL,
                           dist(rbind(dt[data_samePsameL_dt72_pharmacophores[i,"name_pock1"],],
                                      dt[data_samePsameL_dt72_pharmacophores[i,"name_pock2"],]), 
                                method = "manhattan"))
}
save(dist_lig_sameP_sameL, file = "../results/pharmacophores_results/dist_data_intersect72/dist_lig_sameP_sameL_manhattan.Rdata")
print("dist_lig_sameP_sameL dist manhattan")
##diffPsameL
dist_lig_diffP_sameL = NULL
for (i in 1:nrow(data_diffPsameL_dt72_pharmacophores)) {
  dist_lig_diffP_sameL = c(dist_lig_diffP_sameL,
                           dist(rbind(dt[data_diffPsameL_dt72_pharmacophores[i,"name_pock1"],],
                                      dt[data_diffPsameL_dt72_pharmacophores[i,"name_pock2"],]), 
                                method = "manhattan"))
}
save(dist_lig_diffP_sameL, file = "../results/pharmacophores_results/dist_data_intersect72/dist_lig_diffP_sameL_manhattan.Rdata")
print("dist_lig_diffP_sameL dist manhattan")
##diffPdiffL
dist_lig_diffP_diffL = NULL
for (i in 1:nrow(data_diffPdiffL_dt72_pharmacophores)) {
  dist_lig_diffP_diffL = c(dist_lig_diffP_diffL,
                           dist(rbind(dt[data_diffPdiffL_dt72_pharmacophores[i,"name_pock1"],],
                                      dt[data_diffPdiffL_dt72_pharmacophores[i,"name_pock2"],]), 
                                method = "manhattan"))
}
save(dist_lig_diffP_diffL, file = "../results/pharmacophores_results/dist_data_intersect72/dist_lig_diffP_diffL_manhattan.Rdata")
print("dist_lig_diffP_diffL dist manhattan")
##############################"
print("###Distance fuzcav perso n cerisier ###")
##samePsameL
dist_lig_sameP_sameL = NULL
for (i in 1:nrow(data_samePsameL_dt72_pharmacophores)) {
  dist_lig_sameP_sameL = c(dist_lig_sameP_sameL,
                           dist_fuzcav_ph(as.integer(dt[data_samePsameL_dt72_pharmacophores[i,"name_pock1"],]),
                                          as.integer(dt[data_samePsameL_dt72_pharmacophores[i,"name_pock2"],])))
}
save(dist_lig_sameP_sameL, file = "../results/pharmacophores_results/dist_data_intersect72/dist_lig_sameP_sameL_fuzcavPerso.Rdata")
print("dist_lig_sameP_sameL dist fuzcav cerisier")
##diffPsameL
dist_lig_diffP_sameL = NULL
for (i in 1:nrow(data_diffPsameL_dt72_pharmacophores)) {
  dist_lig_diffP_sameL = c(dist_lig_diffP_sameL,
                           dist_fuzcav_ph(as.integer(dt[data_diffPsameL_dt72_pharmacophores[i,"name_pock1"],]),
                                          as.integer(dt[data_diffPsameL_dt72_pharmacophores[i,"name_pock2"],])))
  
}
save(dist_lig_diffP_sameL, file = "../results/pharmacophores_results/dist_data_intersect72/dist_lig_diffP_sameL_fuzcavPerso.Rdata")
print("dist_lig_diffP_sameL dist fuzcav cerisier")
##diffPdiffL
dist_lig_diffP_diffL = NULL
for (i in 1:nrow(data_diffPdiffL_dt72_pharmacophores)) {
  dist_lig_diffP_diffL = c(dist_lig_diffP_diffL,
                           dist_fuzcav_ph(as.integer(dt[data_diffPdiffL_dt72_pharmacophores[i,"name_pock1"],]),
                                          as.integer(dt[data_diffPdiffL_dt72_pharmacophores[i,"name_pock2"],])))
  
}
save(dist_lig_diffP_diffL, file = "../results/pharmacophores_results/dist_data_intersect72/dist_lig_diffP_diffL_fuzcavPerso.Rdata")
print("dist_lig_diffP_diffL dist fuzcav cerisier")
##############################"
print("###Distance fuzcav la vrai ###")
##samePsameL
dist_lig_sameP_sameL = NULL
for (i in 1:nrow(data_samePsameL_dt72_pharmacophores)) {
  dist_lig_sameP_sameL = c(dist_lig_sameP_sameL,
                           dist_fuzcav_ph_norm(as.integer(dt[data_samePsameL_dt72_pharmacophores[i,"name_pock1"],]),
                                          as.integer(dt[data_samePsameL_dt72_pharmacophores[i,"name_pock2"],])))
}
save(dist_lig_sameP_sameL, file = "../results/pharmacophores_results/dist_data_intersect72/dist_lig_sameP_sameL_fuzcaVrai.Rdata")
print("dist_lig_sameP_sameL dist fuzcav vrai")
##diffPsameL
dist_lig_diffP_sameL = NULL
for (i in 1:nrow(data_diffPsameL_dt72_pharmacophores)) {
  dist_lig_diffP_sameL = c(dist_lig_diffP_sameL,
                           dist_fuzcav_ph_norm(as.integer(dt[data_diffPsameL_dt72_pharmacophores[i,"name_pock1"],]),
                                          as.integer(dt[data_diffPsameL_dt72_pharmacophores[i,"name_pock2"],])))
  
}
save(dist_lig_diffP_sameL, file = "../results/pharmacophores_results/dist_data_intersect72/dist_lig_diffP_sameL_fuzcaVrai.Rdata")
print("dist_lig_diffP_sameL dist fuzcav vrai")
##diffPdiffL
dist_lig_diffP_diffL = NULL
for (i in 1:nrow(data_diffPdiffL_dt72_pharmacophores)) {
  dist_lig_diffP_diffL = c(dist_lig_diffP_diffL,
                           dist_fuzcav_ph_norm(as.integer(dt[data_diffPdiffL_dt72_pharmacophores[i,"name_pock1"],]),
                                          as.integer(dt[data_diffPdiffL_dt72_pharmacophores[i,"name_pock2"],])))
  
}
save(dist_lig_diffP_diffL, file = "../results/pharmacophores_results/dist_data_intersect72/dist_lig_diffP_diffL_fuzcaVrai.Rdata")
print("dist_lig_diffP_diffL dist fuzcav vrai")