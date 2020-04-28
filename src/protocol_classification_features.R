### PROTOCOL CLASSIFICATION ###
## DATA PREPARATION 72 desc ##
#Libraries
library(caret)
library(corrplot)
library(car)
library(FactoMineR)
library(sm)
library(colorspace)
####load data####
dt_72descriptors = read.table("../data/data_PDB_72desc.txt", header = T, sep = "", row.names = 1, fill = TRUE)
dt_12descriptors = read.table("../data/data_desc.csv", header = T, sep = ",", row.names = 1, nrow = 2)
##DT cleaning##
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
###
#Management dt_72
rownames(dt_72descriptors)
dt_72descriptors[1:2,1:10]
dt = dt_72descriptors
### 1: supprime features NULL
nrow(dt_72descriptors)
#summary(dt_72descriptors)
#dt[,"HEIGHT_MIN_CYLINDER"] = NULL
#dt[,"RADIUS_MIN_CYLINDER"] = NULL

#dt[,"VOLUME_HULL"] = NULL
#dt[,"drugg"] = NULL
### 2: selection 14 60 residues poches  minimum_size_pocket = 60
maximum_size_pocket = 14
minimum_size_pocket = 60
inf_60 = which(dt[,"C_RESIDUES"] <= minimum_size_pocket)
print("number of pockets superior 60:")
print(nrow(dt)-length(inf_60))
sup_14 = which(dt[,"C_RESIDUES"] >= maximum_size_pocket)
print("number of pockets inferior 14:")
print(nrow(dt)-length(sup_14))

print("Number of pockets >=14 <=60")
print(length(intersect(inf_60,sup_14)))
dt = dt[intersect(inf_60,sup_14),]
#ligands dans jeux données
names_prot = sapply(strsplit(rownames(dt), "_"), "[", 1)
names_ligand = sapply(strsplit(rownames(dt), "_"), "[", 2)
### 3: supprimer poches lié aux ligand pas interressants
excluded_ligands = scan("../data/PLIC_lig/excluded_ligands.txt", sep ="\n", what = character())
dt = dt[!is.element(names_ligand, excluded_ligands),]
nrow(dt)
summary(dt)
### 4: supprimer NA
#print("number pockets NA")
#print(nrow(dt) - nrow(na.omit(dt)))  
#dt = na.omit(dt)
### 5: scale
dt = scale(dt)
#again
names_prot = sapply(strsplit(rownames(dt), "_"), "[", 1)
names_ligand = sapply(strsplit(rownames(dt), "_"), "[", 2)
#
length(unique(names_ligand))
### INTERSECT PHARMACOPHORES ###
dt_pharmacophores = read.table("../data/FPCount_save_all_inter_dt72.txt", 
                               sep = ";", row.names = 1, colClasses = c(rep("charachter", 1), rep("NULL", 6)))
dt_pharmacophores[,1:20]

###########################################
descriptors_hydrophobicity =  c(
                                "p_hydrophobic_residues",
                                "p_hyd_atom",
                                #"hydrophobicity_pocket", not found
                                "hydrophobic_kyte"
                                )
descriptors_aromatic =        c(
                                "p_aromatic_residues",
                                "p_Car_atom"
                                )
descriptors_polarity =        c(
                                "p_polar_residues"
                                #"polarity"#with NA
                                )
descriptors_physicochemical = c(
                                "p_charged_residues",
                                "p_negative_residues",
                                "p_positive_residues",
                                "charge",
                                "p_aliphatic_residues",
                                "p_Nlys_atom",
                                "p_Ntrp_atom",
                                "p_Ocoo_atom",
                                "p_Cgln_atom",
                                "p_Ccoo_atom",
                                "p_Carg_atom",
                                "p_S_atom",
                                "p_Otyr_atom",
                                "p_Ooh_atom",
                                "p_O_atom",
                                "p_N_atom",
                                #"p_ND1_atom",
                                "p_NE2_atom",
                                #"p_pro_residues",
                                "p_small_residues",
                                "p_tiny_residues",
                                "p_main_chain_atom",
                                "p_side_chain_atom",
                                "p_C_atom",
                                "p_sulfur_atom",
                                "p_carbone_atom",
                                "p_oxygen_atom",
                                "p_nitrogen_atom"
                                )
descriptors_geometrical =     c(
                                "RADIUS_HULL",
                                "SURFACE_HULL",
                                "DIAMETER_HULL",
                                "VOLUME_HULL",
                                "FACE",
                                "SMALLEST_SIZE",
                                "RADIUS_CYLINDER",
                                "X._ATOM_CONVEXE",
                                "C_ATOM",
                                "C_RESIDUES",
                                "PSI",
                                "PCI",
                                "CONVEX.SHAPE_COEFFICIENT",
                                "INERTIA_1",
                                "INERTIA_2",
                                "INERTIA_3"
                                )
#
summary(dt[,descriptors_geometrical])
length(descriptors_physicochemical)
setdiff(colnames(dt),c(descriptors_physicochemical,
                       descriptors_geometrical,
                       descriptors_polarity,
                       descriptors_aromatic,
                       descriptors_hydrophobicity))
setdiff(descriptors_geometrical,colnames(dt))

#### FUNCTION###
#DIST POCK same prot same lig
dist_lig_sameP_sameL = function(dt, names_ligand) {
  dist_ligs = rep(0,length(which(table(names_ligand) > 1)))
  names(dist_ligs) = names(which(table(names_ligand) > 1))
  names_ligand_unique = unique(names_ligand)
  flag = 1
  dist_features_pock_all = NULL
  dist_features = matrix(data = 0,
                         nrow = length(which(table(names_ligand) > 1)),
                         ncol = ncol(dt))
  dist_features = as.data.frame(dist_features)
  rownames(dist_features) = names(which(table(names_ligand) > 1))
  data_samePsameL = NULL
  rnames = rownames(dt)
  for (lig_name in names_ligand_unique){#length(names_ligand_unique)
    #print(i)
    index = grep(paste0(paste0("_",lig_name),"_"),rownames(dt))
    if(length(index) > 1) {
      print(flag)
      if(length(index) > 10) {
        index = sample(index,10)
      }
      dist_pock = NULL
      dist_features_pock = NULL
      for (j in 1:(length(index)-1)) {
        for (k in (j+1):length(index)) {
          if(sapply(strsplit(rnames[index[j]], "_"), "[", 1) == sapply(strsplit(rnames[index[k]], "_"), "[", 1)) { #CHANGE TO != FOR DIFFERENT PROTEINS
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
  write.csv(data_samePsameL, "../data/data_structure_comparaison/data_samePsameL.csv")
  #save(dist_features_pock_all, file = "../results/protocol_features/dist_features_sameP_sameL_pock.Rdata")
  #save(dist_features, file = "../results/protocol_features/dist_features_sameP_sameL.Rdata")
  return(dist_ligs)
}
#DIST POCK diff prot same lig
dist_lig_diffP_sameL = function(dt, names_ligand) {
  dist_ligs = rep(0,length(which(table(names_ligand) > 1)))
  names(dist_ligs) = names(which(table(names_ligand) > 1))
  names_ligand_unique = unique(names_ligand)
  flag = 1
  dist_features_pock_all = NULL
  dist_features = matrix(data = 0,
                         nrow = length(which(table(names_ligand) > 1)),
                         ncol = ncol(dt))
  dist_features = as.data.frame(dist_features)
  rownames(dist_features) = names(which(table(names_ligand) > 1))
  #save file active samePsameL
  data_diffPsameL = NULL
  rnames = rownames(dt)
  for (lig_name in names_ligand_unique){
    #print(i)
    index = grep(paste0(paste0("_",lig_name),"_"),rownames(dt))
    if(length(index) > 1) {
      print(flag)
      if(length(index) > 10) {
        index = sample(index,10)
      }
      dist_pock = NULL
      dist_features_pock = NULL
      for (j in 1:(length(index)-1)) {
        for (k in (j+1):length(index)) {
          if(sapply(strsplit(rnames[index[j]], "_"), "[", 1) != sapply(strsplit(rnames[index[k]], "_"), "[", 1)) { #CHANGE TO != FOR DIFFERENT PROTEINS
            data_diffPsameL = rbind(data_diffPsameL, c(rnames[index[j]],rnames[index[k]]))
            #dist_pock = c(dist_pock, dist(rbind(dt[index[j],],
            #                                    dt[index[k],])))
            #dist_features_pock = rbind(dist_features_pock,  apply(dt[c(index[j],index[k]),],2,dist))
          }
        }
      }
      flag = flag+1
      if(!is.null(dist_pock)) {
        #dist_ligs[lig_name] = mean(dist_pock)
        #dist_features_pock_all = rbind(dist_features_pock_all,dist_features_pock)
        #dist_features[lig_name,] = apply(dist_features_pock,2,mean)
        ##flag = flag+1      
      }
    }
  }
  data_diffPsameL = as.data.frame(data_diffPsameL)
  colnames(data_diffPsameL) = c("name_pock1","name_pock2")
  write.csv(data_diffPsameL, "../data/data_structure_comparaison/data_diffPsameL.csv")
  #save(dist_features_pock_all, file = "../results/protocol_features/dist_features_diffP_sameL_pock.Rdata")
  #save(dist_features, file = "../results/protocol_features/dist_features_diffP_sameL.Rdata")
  return(dist_ligs)
}
#DIST POCK diff prot diff lig
dist_lig_diffP_diffL = function(dt, names_ligand, n = 5000) {
  unique_names_ligand = unique(names_ligand)
  dist_ligs_random = rep(0,n)#length(unique(names_ligand)))
  dist_features = NULL
  dist_features_pock_all = NULL
  data_diffPdiffL = NULL
  rnames = rownames(dt)
  for (i in 1:n){#length(unique(names_ligand))){
    print(i)
    lig = sample(unique_names_ligand, 2)
    i_pock1 = sample(grep(paste0(paste0("_",lig[1]),"_"),rnames),1)
    i_pock2 = sample(grep(paste0(paste0("_",lig[2]),"_"),rnames),1)
    #dist_ligs_random[i] = dist(rbind(dt[i_pock1,],
    #                                 dt[i_pock2,]))
    data_diffPdiffL = rbind(data_diffPdiffL, c(rnames[i_pock1],rnames[i_pock2]))
    #dist_features = rbind(dist_features,  apply(dt[c(i_pock1,i_pock2),],2,dist))
     
  }
  data_diffPdiffL = as.data.frame(data_diffPdiffL)
  colnames(data_diffPdiffL) = c("name_pock1","name_pock2")
  write.csv(data_diffPdiffL, "../data/data_structure_comparaison/data_diffPdiffL.csv")
  #dist_features_pock_all = rbind(dist_features_pock_all,dist_features)
  #save(dist_features_pock_all, file = "../results/protocol_features/dist_features_diffP_diffL_pock.Rdata")
  #save(dist_features, file = "../results/protocol_features/dist_features_diffP_diffL.Rdata")
  return(dist_ligs_random)
}
###########################
names_prot = sapply(strsplit(rownames(dt), "_"), "[", 1)
names_ligand = sapply(strsplit(rownames(dt), "_"), "[", 2)

### Distance plot ###
plot_distance = function(dt, names_ligand, features_name) {
  dist_lig_sameP_sameL = dist_lig_sameP_sameL(dt, names_ligand)
  #save(dist_lig_sameP_sameL, file = paste0(paste0("../results/protocol_features/", features_name),
  #                                         "/dist_lig_sameP_sameL.Rdata"))
  dist_lig_diffP_sameL = dist_lig_diffP_sameL(dt, names_ligand)
  #save(dist_lig_diffP_sameL, file = paste0(paste0("../results/protocol_features/", features_name),
  #                                         "/dist_lig_diffP_sameL.Rdata"))  
  dist_lig_diffP_diffL = dist_lig_diffP_diffL(dt, names_ligand)
  #save(dist_lig_diffP_diffL, file = paste0(paste0("../results/protocol_features/", features_name),
  #                                         "/dist_lig_diffP_diffL.Rdata"))
  png(paste0(paste0("../results/protocol_features/", features_name),
             "/distance_pockets.png"))
  sm.density.compare(c(dist_lig_diffP_diffL,
                       dist_lig_sameP_sameL[which(dist_lig_sameP_sameL>0)],
                       dist_lig_diffP_sameL[which(dist_lig_diffP_sameL>0)]),
                     c(rep(1,length(dist_lig_diffP_diffL)),
                       rep(2,length(dist_lig_sameP_sameL[which(dist_lig_sameP_sameL>0)])),
                       rep(3,length(dist_lig_diffP_sameL[which(dist_lig_diffP_sameL>0)]))
                     ),
                     model = "none"
                     , xlab = "Mean distance bewteen pockets")
  dev.off()
}

#############################################################
#save distance
nrow(dt)
features_random_forest = c(
            "SURFACE_HULL",
            "VOLUME_HULL",
            "RADIUS_HULL",
            "RADIUS_CYLINDER",
            "DIAMETER_HULL",
            "C_ATOM",
            "p_side_chain_atom",
            "p_main_chain_atom",
            "p_O_atom",
            "p_main_chain_atom",
            "PCI",
            "PSI",
            "X._ATOM_CONVEXE",
            "p_carbone_atom",
            "p_nitrogen_atom"
            )
features_selection_interet = c(
  "RADIUS_HULL",
  "SURFACE_HULL",
  "DIAMETER_HULL",
  "VOLUME_HULL",
  "RADIUS_CYLINDER",
  "C_RESIDUES",
  "hydrophobic_kyte",
  "SMALLEST_SIZE",
  "X._ATOM_CONVEXE",
  "p_aliphatic_residues",
  "p_main_chain_atom",
  "p_side_chain_atom",
  "p_polar_residues",
  "p_aromatic_residues",
  "p_Car_atom",
  "p_hydrophobic_residues",
  "p_hyd_atom"
)
features_selection_interet_cor9 = c(
  "SURFACE_HULL",
  "VOLUME_HULL",
  "C_RESIDUES",
  "hydrophobic_kyte",
  "SMALLEST_SIZE",
  "X._ATOM_CONVEXE",
  "p_aliphatic_residues",
  "p_main_chain_atom",
  "p_side_chain_atom",
  "p_polar_residues",
  "p_aromatic_residues",
  "p_Car_atom",
  "p_hydrophobic_residues",
  "p_hyd_atom"
)
features_best_geometrique =  c(
  "RADIUS_HULL",
  "SURFACE_HULL",
  "DIAMETER_HULL",
  "VOLUME_HULL",
  "RADIUS_CYLINDER",
  "C_RESIDUES",
  "SMALLEST_SIZE",
  "X._ATOM_CONVEXE"
)
features_best_physicochimique = c(
  "hydrophobic_kyte",
  "p_aliphatic_residues",
  "p_main_chain_atom",
  "p_side_chain_atom",
  "p_polar_residues",
  "p_aromatic_residues",
  "p_Car_atom",
  "p_hydrophobic_residues",
  "p_hyd_atom",
  "p_charged_residues"
)
features_features_tandencieux = c(features_selection_interet,
                                  "PSI",
                                  "INERTIA_1",
                                  "p_charged_residues",
                                  "charge",
                                  "p_small_residues",
                                  "p_tiny_residues",
                                  "p_nitrogen_atom",
                                  "p_negative_residues",
                                  "C_ATOM"
                               )
c(descriptors_hydrophobicity,
  descriptors_aromatic,
  descriptors_polarity,
  descriptors_physicochemical,
  descriptors_geometrical)
features_name = c("features_selection_interet_cor9_radiushull")
features = list(c("RADIUS_HULL",features_selection_interet_cor9))
flag = 1
for (f in features) {
  print(features_name[flag])
  plot_distance(dt[,f],names_ligand = names_ligand, features_name = features_name[flag])
  flag = flag+1
}

#plot
load("../results/protocol_features/features_selection_interet//dist_lig_sameP_sameL.Rdata")
load("../results/protocol_features/features_selection_interet/dist_lig_diffP_sameL.Rdata")
load("../results/protocol_features/features_selection_interet/dist_lig_diffP_diffL.Rdata")
png(paste0(paste0("../results/protocol_features/", features_name[2]),
           "/distance_pockets.png"))
sm.density.compare(c(dist_lig_diffP_diffL,
                     dist_lig_sameP_sameL[which(dist_lig_sameP_sameL>0)],
                     dist_lig_diffP_sameL[which(dist_lig_diffP_sameL>0)]),
                   c(rep(1,length(dist_lig_diffP_diffL)),
                     rep(2,length(dist_lig_sameP_sameL[which(dist_lig_sameP_sameL>0)])),
                     rep(3,length(dist_lig_diffP_sameL[which(dist_lig_diffP_sameL>0)]))
                   ),
                   model = "none"
                   , xlab = "Mean distance bewteen pockets")

dev.off()
#physico + geo
#P
load("../results/protocol_features/features_physicochimiques/dist_lig_sameP_sameL.Rdata")
load("../results/protocol_features/features_physicochimiques/dist_lig_diffP_sameL.Rdata")
load("../results/protocol_features/features_physicochimiques/dist_lig_diffP_diffL.Rdata")
dist_lig_sameP_sameL_P = dist_lig_sameP_sameL
dist_lig_diffP_sameL_P = dist_lig_diffP_sameL
dist_lig_diffP_diffL_P = dist_lig_diffP_diffL
#G
load("../results/protocol_features/features_geometriques/dist_lig_sameP_sameL.Rdata")
load("../results/protocol_features/features_geometriques/dist_lig_diffP_sameL.Rdata")
load("../results/protocol_features/features_geometriques/dist_lig_diffP_diffL.Rdata")
dist_lig_sameP_sameL_G = dist_lig_sameP_sameL
dist_lig_diffP_sameL_G = dist_lig_diffP_sameL
dist_lig_diffP_diffL_G = dist_lig_diffP_diffL
#P+G
load("../results/protocol_features/features_geo+phys/dist_lig_sameP_sameL.Rdata")
load("../results/protocol_features/features_geo+phys/dist_lig_diffP_sameL.Rdata")
load("../results/protocol_features/features_geo+phys/dist_lig_diffP_diffL.Rdata")
#12
load("../results/protocol_features/features_selection_12/dist_lig_sameP_sameL.Rdata")
load("../results/protocol_features/features_selection_12/dist_lig_diffP_sameL.Rdata")
load("../results/protocol_features/features_selection_12/dist_lig_diffP_diffL.Rdata")
#best
load("../results/protocol_features/features_best_physicochimique/dist_lig_sameP_sameL.Rdata")
load("../results/protocol_features/features_best_physicochimique/dist_lig_diffP_sameL.Rdata")
load("../results/protocol_features/features_best_physicochimique/dist_lig_diffP_diffL.Rdata")
#
#all
load("../results/protocol_features/features_all/dist_lig_sameP_sameL.Rdata")
load("../results/protocol_features/features_all/dist_lig_diffP_sameL.Rdata")
load("../results/protocol_features/features_all/dist_lig_diffP_diffL.Rdata")

sm.density.compare(c(dist_lig_diffP_diffL,
                     dist_lig_sameP_sameL[which(dist_lig_sameP_sameL>0)],
                     dist_lig_diffP_sameL[which(dist_lig_diffP_sameL>0)]),
                   c(rep(1,length(dist_lig_diffP_diffL)),
                     rep(2,length(dist_lig_sameP_sameL[which(dist_lig_sameP_sameL>0)])),
                     rep(3,length(dist_lig_diffP_sameL[which(dist_lig_diffP_sameL>0)]))
                   ),
                   model = "none", xlim = c(0,10)
                   , xlab = "Mean distance bewteen pockets")

### BAR PLOT FEATURES ###
f = c(descriptors_hydrophobicity,
      descriptors_aromatic,
      descriptors_polarity,
      descriptors_physicochemical,
      descriptors_geometrical)
load("../results/protocol_features/dist_features_sameP_sameL.Rdata")
dist_featuresD_sameP_sameL = dist_features
colnames(dist_featuresD_sameP_sameL) = f
dist_features_sameP_sameL = apply(dist_featuresD_sameP_sameL[which(apply(dist_featuresD_sameP_sameL[,1:10],1,sum) != 0),],2,mean)
dist_features_sameP_sameL_sd = apply(dist_featuresD_sameP_sameL[which(apply(dist_featuresD_sameP_sameL[,1:10],1,sum) != 0),],2,sd)
load("../results/protocol_features/dist_features_diffP_sameL.Rdata")
dist_featuresD_diffP_sameL = dist_features
colnames(dist_featuresD_diffP_sameL) = f
dist_features_diffP_sameL = apply(dist_featuresD_diffP_sameL[which(apply(dist_featuresD_diffP_sameL[,1:10],1,sum) != 0),],2,mean)
dist_features_diffP_sameL_sd = apply(dist_featuresD_diffP_sameL[which(apply(dist_featuresD_diffP_sameL[,1:10],1,sum) != 0),],2,sd)
load("../results/protocol_features/dist_features_diffP_diffL.Rdata")
dist_featuresD_diffP_diffL = dist_features
dist_features_diffP_diffL = apply(dist_features,2,mean)
dist_features_diffP_diffL_sd = apply(dist_features,2,sd)

dt_dist_features_mean = apply(dist_features,2,mean)
names(dt_dist_features_mean) = tolower(names(dt_dist_features_mean))
#
par(cex=1, mar=c(12,4,4,2))
#
barplot(dt_dist_features_mean, las=2)
#
library(ggplot2)
dt_dist_features_mean = data.frame(features = rep(names(dist_features_diffP_sameL),3),
                                   dist_features = c(dist_features_sameP_sameL,
                                                     dist_features_diffP_sameL,
                                                     dist_features_diffP_diffL),
                                   sd_features = c(dist_features_sameP_sameL_sd,
                                                   dist_features_diffP_sameL_sd,
                                                   dist_features_diffP_diffL_sd),
                                   dist_type = c(rep("dist_features_sameP_sameL",length(dist_features_sameP_sameL)),
                                                 rep("dist_features_diffP_sameL",length(dist_features_diffP_sameL)),
                                                 rep("dist_features_diffP_diffL",length(dist_features_diffP_diffL))))

ggplot(dt_dist_features_mean, aes(x = features, y = dist_features, fill=dist_type)) +
  geom_bar(stat = "identity",  position=position_dodge()) +
  geom_errorbar(aes(ymin=dist_features,
                    ymax=dist_features+sd_features),
                width=.2, position=position_dodge(0.9)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

barplot(dist_features_diffP_diffL - dist_features_diffP_sameL, las = 2)
dev.off()
### bar plot dist > 6 ###
load("../results/protocol_features/features_best_physicochimique/dist_lig_sameP_sameL.Rdata")
which(dist_lig_sameP_sameL > 2.5)
load("../results/protocol_features/features_best_physicochimique/dist_lig_diffP_sameL.Rdata")
which(dist_lig_diffP_sameL > 2.5)
#features
#mean
dist_features_sameP_sameL = apply(dist_featuresD_sameP_sameL[which(dist_lig_sameP_sameL > 6),],2,mean)
dist_features_diffP_sameL = apply(dist_featuresD_diffP_sameL[which(dist_lig_diffP_sameL > 6),],2,mean)
dist_features_diffP_diffL = apply(dist_featuresD_diffP_diffL,2,mean)
#sd
dist_features_sameP_sameL_sd = apply(dist_featuresD_sameP_sameL[which(dist_lig_sameP_sameL > 6),],2,sd)
dist_features_diffP_sameL_sd = apply(dist_featuresD_diffP_sameL[which(dist_lig_diffP_sameL > 6),],2,sd)
dist_features_diffP_diffL_sd = apply(dist_featuresD_diffP_diffL,2,sd)
#
dt_dist_features_mean = data.frame(features = rep(names(dist_features_diffP_sameL),3),
                                   dist_features = c(dist_features_sameP_sameL,
                                                     dist_features_diffP_sameL,
                                                     dist_features_diffP_diffL),
                                   sd_features = c(dist_features_sameP_sameL_sd,
                                                   dist_features_diffP_sameL_sd,
                                                  dist_features_diffP_diffL_sd),
                                   dist_type = c(rep("dist_features_sameP_sameL",length(dist_features_sameP_sameL)),
                                                 rep("dist_features_diffP_sameL",length(dist_features_diffP_sameL)),
                                                 rep("dist_features_diffP_diffL",length(dist_features_diffP_diffL))))

ggplot(dt_dist_features_mean, aes(x = features, y = dist_features, fill=dist_type)) +
  geom_bar(stat = "identity",  position=position_dodge()) +
  geom_errorbar(aes(ymin=dist_features,
                    ymax=dist_features+sd_features),
                width=.2, position=position_dodge(0.9)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

### PLOT BIAIS if n pocket attached to ligand > 100 ###
load("../results/protocol_features/dist_features_sameP_sameL_pock.Rdata")
dist_featuresD_sameP_sameL = dist_features_pock_all
dist_features_sameP_sameL = apply(dist_features_pock_all,2,mean)
dist_features_sameP_sameL_sd = apply(dist_features_pock_all,2,sd)
load("../results/protocol_features/dist_features_diffP_sameL_pock.Rdata")
dist_featuresD_diffP_sameL = dist_features_pock_all
dist_features_diffP_sameL = apply(dist_features_pock_all,2,mean)
dist_features_diffP_sameL_sd = apply(dist_features_pock_all,2,sd)
load("../results/protocol_features/dist_features_diffP_diffL_pock.Rdata")
dist_featuresD_diffP_diffL = dist_features_pock_all
dist_features_diffP_diffL = apply(dist_features_pock_all,2,mean)
dist_features_diffP_diffL_sd = apply(dist_features_pock_all,2,sd)
#
#
library(ggplot2)
dt_dist_features_mean = data.frame(features = rep(names(dist_features_diffP_sameL),3),
                                   dist_features = c(dist_features_sameP_sameL,
                                                     dist_features_diffP_sameL,
                                                     dist_features_diffP_diffL),
                                   sd_features = c(dist_features_sameP_sameL_sd,
                                                   dist_features_diffP_sameL_sd,
                                                   dist_features_diffP_diffL_sd),
                                   dist_type = c(rep("dist_features_sameP_sameL",length(dist_features_sameP_sameL)),
                                                 rep("dist_features_diffP_sameL",length(dist_features_diffP_sameL)),
                                                 rep("dist_features_diffP_diffL",length(dist_features_diffP_diffL))))
#
ggplot(dt_dist_features_mean, aes(x = features, y = dist_features, fill=dist_type)) +
  geom_bar(stat = "identity",  position=position_dodge()) +
  geom_errorbar(aes(ymin=dist_features-sd_features,
                    ymax=dist_features+sd_features),
                width=.2, position=position_dodge(0.9)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

#
par(cex=1, mar=c(12,4,4,2))
### T.TEST ###
mydataframe_dist_feature = NULL
mydataframe_dist_feature = rbind(dist_featuresD_sameP_sameL,
                                 rbind(dist_featuresD_diffP_sameL,
                                       dist_featuresD_diffP_diffL))
mydataframe_dist_feature = cbind(mydataframe_dist_feature,c(rep(1,nrow(dist_featuresD_sameP_sameL)),
                                                            rep(1,nrow(dist_featuresD_diffP_sameL)),
                                                            rep(0,nrow(dist_featuresD_diffP_diffL))))
mydataframe_dist_feature = as.data.frame(mydataframe_dist_feature)
##diffPdiffL VS samePsameL avec ou sans BIAIS
p_value_diffPdiffL_samePsameL_biais100 = rep(0,length(colnames(dist_featuresD_sameP_sameL)))
names(p_value_diffPdiffL_samePsameL_biais100) = colnames(dist_featuresD_sameP_sameL)
for (descriptor in colnames(dist_featuresD_sameP_sameL)) {
  p_value_diffPdiffL_samePsameL_biais100[descriptor] = t.test(dist_featuresD_sameP_sameL[,descriptor],
                                                              dist_featuresD_diffP_diffL[,descriptor])$p.value
}
sort(p_value_diffPdiffL_samePsameL_biais100)
##diffPdiffL VS samePsameL avec ou sans BIAIS
p_value_diffPdiffL_diffPsameL_biais100 = rep(0,length(colnames(dist_featuresD_diffP_sameL)))
names(p_value_diffPdiffL_samePsameL_biais100) = colnames(dist_featuresD_sameP_sameL)
for (descriptor in colnames(dist_featuresD_diffP_sameL)) {
  p_value_diffPdiffL_diffPsameL_biais100[descriptor] = t.test(dist_featuresD_diffP_sameL[,descriptor],
                                                              dist_featuresD_diffP_diffL[,descriptor])$p.value
}
which(p_value_diffPdiffL_diffPsameL_biais100 > 0)
sort(p_value_diffPdiffL_diffPsameL_biais100)
names(p_value_diffPdiffL_diffPsameL_biais100)
#
par(cex=1, mar=c(12,4,4,2))
#
plot(scale(p_value_diffPdiffL_diffPsameL_biais100))
dev.off()
#COR
res <- apply(mydataframe_dist_feature[,1:48],2, function(x) {x1 <-cor(x,mydataframe_dist_feature$V49)})
par(cex=1, mar=c(12,4,4,2))
barplot(res, las = 2)
dev.off()
###ROC CURVE ###
load("../results/protocol_features/features_selection_interet_cor9_radiushull/dist_lig_sameP_sameL.Rdata")
load("../results/protocol_features/features_selection_interet_cor9_radiushull/dist_lig_diffP_sameL.Rdata")
load("../results/protocol_features/features_selection_interet_cor9_radiushull/dist_lig_diffP_diffL.Rdata")
#
##SCALE
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
#
min(c(dist_lig_sameP_sameL[which(dist_lig_sameP_sameL>0)],
      dist_lig_diffP_sameL[which(dist_lig_diffP_sameL>0)],
      dist_lig_diffP_diffL))
max(c(dist_lig_sameP_sameL,dist_lig_diffP_sameL,dist_lig_diffP_diffL))
#
dist_lig_diffP_sameL = dist_ligs[which(dist_ligs > 0)]
dist_lig_diffP_diffL = dist_ligs_random
# AUC  diffLdiffP vs sameLsameP
y_true = c(rep(0,length(dist_lig_sameP_sameL[which(dist_lig_sameP_sameL>0)])),
           rep(1,length(dist_lig_diffP_diffL)))
y_predict = c(range01(dist_lig_sameP_sameL[which(dist_lig_sameP_sameL>0)]),
              range01(dist_lig_diffP_diffL))
#AUC diffLdiffP vs sameLdiffP
y_true = c(rep(0,length(dist_lig_diffP_sameL[which(dist_lig_diffP_sameL>0)])),
           rep(1,length(dist_lig_diffP_diffL)))
y_predict = c(range01(dist_lig_diffP_sameL[which(dist_lig_diffP_sameL>0)]),
              range01(dist_lig_diffP_diffL))
#AUC  diffLdiffP vs (sameLsameP|sameLdiffP)
y_true = c(rep(0,length(dist_lig_sameP_sameL[which(dist_lig_sameP_sameL>0)])+
                 length(dist_lig_diffP_sameL[which(dist_lig_diffP_sameL>0)])),
           rep(1,length(dist_lig_diffP_diffL)))
y_predict = c(range01(dist_lig_sameP_sameL[which(dist_lig_sameP_sameL>0)]),
              range01(dist_lig_diffP_sameL[which(dist_lig_diffP_sameL>0)]),
              range01(dist_lig_diffP_diffL))
#library(ROCR)
dt.pred = prediction(y_predict, y_true)
dt.perf = performance(dt.pred, "tpr", "fpr")
plot(dt.perf)
dt.auc = performance(dt.pred, "auc")
attr(dt.auc, "y.values")
#

##load data BENCHMARK ##
#dataset : review_structures
dt_review_structures = read.table("../data/benchmark/review_structures.tar/review_structures/review_structures/data_review_structures.txt", header = TRUE, row.names = 1, sep = "")
#
review_structures = read.csv("../data/benchmark/review_structures.tar/review_structures/review_structures/review_structures.csv", header = FALSE)
review_structures[,1] = substr(toupper(review_structures[,1]),1,4)
review_structures[,2] = substr(toupper(review_structures[,2]),1,4)
#
dt_benchmark = review_structures##to load before
dt_structures = dt_review_structures
#dataset : kahraman_structures
dt_kahraman_structures = read.table("../data/benchmark/kahraman_structures.tar/kahraman_structures/kahraman_structures/data_kahraman_structures.txt", header = TRUE, row.names = 1, sep = "")
#
kahraman_structures = read.csv("../data/benchmark/kahraman_structures.tar/kahraman_structures/kahraman_structures/kahraman_structures.csv", header = FALSE)
kahraman_structures[,1] = substr(toupper(kahraman_structures[,1]),1,4)
kahraman_structures[,2] = substr(toupper(kahraman_structures[,2]),1,4)
#
dt_benchmark = kahraman_structures##to load before
dt_structures = dt_kahraman_structures
## PREDICTION ##
names_prot_structures = sapply(strsplit(rownames(dt_structures), "_"), "[", 1)
#SCALE
dt_structures_scale = scale(dt_structures, center = attr(dt, "scaled:center"),
                                           scale = attr(dt, "scaled:scale"))
#
#FEATURES
features = descriptors_geometrical
#
y_true = NULL
y_predict = NULL
for (i in 1:10000) {
  print(i)
  if(is.element(dt_benchmark[i,1], names_prot_structures) &
     is.element(dt_benchmark[i,2], names_prot_structures) & 
     dt_benchmark[i,1] != dt_benchmark[i,2]) {
    if(dt_benchmark[i,3] == "active") {
      y_true = c(y_true,1)
    } else {
      y_true = c(y_true,0)
    }
    pock1 = which(names_prot_structures == dt_benchmark[i,1]) 
    pock2 = which(names_prot_structures == dt_benchmark[i,2]) 
    y_predict = c(y_predict,dist(rbind(dt_structures_scale[pock1,features],
                                       dt_structures_scale[pock2,features])))
  }
}


## correlation descripteurs ##
cor(dt[,c(descriptors_hydrophobicity,
          descriptors_aromatic,
          descriptors_polarity,
          descriptors_physicochemical,
          descriptors_geometrical)])
corrplot(abs(cor(dt[,c(descriptors_hydrophobicity,
                   descriptors_aromatic,
                   descriptors_polarity,
                   descriptors_physicochemical)])))
dev.off()


cor(dt[,c(descriptors_geometrical)])
findCorrelation(abs(cor(dt[,c(descriptors_hydrophobicity,
                            descriptors_aromatic,
                            descriptors_polarity,
                            descriptors_physicochemical)])), cutoff = 0.9, verbose = T, names = T)


