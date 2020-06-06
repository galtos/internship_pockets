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
#dt overlap
dt_overlap = read.table("../data/data_PDB_72desc_overlap_all.txt", header = T, sep = "", row.names = 1, fill = TRUE)
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
###---OVERLAP-----###dt = dt_overlap
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
length(unique(names_prot))
#save scale infos
#write.table(attr(dt, "scaled:center"), file = "../results/scale/scaled_center_dt72clean.Rdata")
#write.table(attr(dt, "scaled:scale"), file = "../results/scale/scaled_scale_dt72clean.Rdata")
#save min max
#write.table(rbind(apply(dt,2,max),apply(dt,2,min)), file = "../results/scale/minmax_value_dt72clean.Rdata")
#
### Management scale OVERLAP ###
scaled_center_dt_t = read.table(file = "../results/scale/scaled_center_dt72clean.Rdata", col.names = F, row.names = 1)
scaled_scale_dt_t = read.table(file = "../results/scale/scaled_scale_dt72clean.Rdata", col.names = F, row.names = 1)
scaled_center_dt = scaled_center_dt_t[,1]
names(scaled_center_dt) = rownames(scaled_center_dt_t)
scaled_scale_dt = scaled_scale_dt_t[,1]
names(scaled_scale_dt) = rownames(scaled_scale_dt_t)
#
dt = scale(dt, scaled_center_dt, scaled_scale_dt)
### INTERSECT PHARMACOPHORES ###
#dt_pharmacophores = read.table("../data/FPCount_save_all_inter_dt72.txt",sep = ";", row.names = 1)
load("../data/intersect_dt72_pharmacophores.Rdata")
dt = dt[intersect_dt72_pharmacophores,]
#again
names_prot = sapply(strsplit(rownames(dt), "_"), "[", 1)
names_ligand = sapply(strsplit(rownames(dt), "_"), "[", 2)
#
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
  write.csv(data_samePsameL, "../data/data_structure_comparaison/data_samePsameL_dt72_pharmacophores.csv")
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
  write.csv(data_diffPsameL, "../data/data_structure_comparaison/data_diffPsameL_dt72_pharmacophores.csv")
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
  write.csv(data_diffPdiffL, "../data/data_structure_comparaison/data_diffPdiffL_dt72_pharmacophores_50000.csv")
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
features_12 = colnames(dt_12descriptors)
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
### DISTANCE AVEC DATASET ###
data_samePsameL_dt72_pharmacophores = read.csv("../data/data_structure_comparaison/data_samePsameL.csv")
data_diffPsameL_dt72_pharmacophores = read.csv("../data/data_structure_comparaison/data_diffPsameL.csv")
data_diffPdiffL_dt72_pharmacophores = read.csv("../data/data_structure_comparaison/data_diffPdiffL.csv")
#
data_samePsameL_dt72_pharmacophores = read.csv("../data/data_structure_comparaison/data_samePsameL_dt72_pharmacophores.csv", colClasses = "character")
data_diffPsameL_dt72_pharmacophores = read.csv("../data/data_structure_comparaison/data_diffPsameL_dt72_pharmacophores.csv", colClasses = "character")
data_diffPdiffL_dt72_pharmacophores = read.csv("../data/data_structure_comparaison/data_diffPdiffL_dt72_pharmacophores.csv", colClasses = "character")
data_samePsameL_dt72_pharmacophores$name_pock1
#
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
#
features = c("RADIUS_HULL",features_selection_interet_cor9)
features = c(descriptors_hydrophobicity,
              descriptors_aromatic,
              descriptors_polarity,
              descriptors_physicochemical,
              descriptors_geometrical)
length(features)
features = features_best_geometrique
features = c("SURFACE_HULL","SMALLEST_SIZE","VOLUME_HULL","X._ATOM_CONVEXE","C_RESIDUES","RADIUS_HULL")
features = c(  "SURFACE_HULL",
               "VOLUME_HULL",
               "C_RESIDUES",
               "hydrophobic_kyte",
               "SMALLEST_SIZE",
               "X._ATOM_CONVEXE")

features = features_best_physicochimique
features = features_selection_interet
features = features_12
##samePsameL
dist_lig_sameP_sameL = NULL
for (i in 1:nrow(data_samePsameL_dt72_pharmacophores)) {
  print(i)
  dist_lig_sameP_sameL = c(dist_lig_sameP_sameL,
                           dist(rbind(dt[data_samePsameL_dt72_pharmacophores[i,"name_pock1"],features],
                                      dt[data_samePsameL_dt72_pharmacophores[i,"name_pock2"],features])))
}
names(dist_lig_sameP_sameL) = sapply(strsplit(data_samePsameL_dt72_pharmacophores$name_pock1, "_"), "[", 2)
which(dist_lig_sameP_sameL == 0)
length(which(dist_lig_sameP_sameL == 0))
##diffPsameL
dist_lig_diffP_sameL = NULL
for (i in 1:nrow(data_diffPsameL_dt72_pharmacophores)) {
  print(i)
  dist_lig_diffP_sameL = c(dist_lig_diffP_sameL,
                           dist(rbind(dt[data_diffPsameL_dt72_pharmacophores[i,"name_pock1"],features],
                                      dt[data_diffPsameL_dt72_pharmacophores[i,"name_pock2"],features])))
}
names(dist_lig_diffP_sameL) = sapply(strsplit(data_diffPsameL_dt72_pharmacophores$name_pock1, "_"), "[", 2)
dist_lig_diffP_sameL
mean(dist_lig_diffP_sameL)
##diffPdiffL
dist_lig_diffP_diffL = NULL
for (i in 1:nrow(data_diffPdiffL_dt72_pharmacophores)) {
  print(i)
  dist_lig_diffP_diffL = c(dist_lig_diffP_diffL,
                           dist(rbind(dt[data_diffPdiffL_dt72_pharmacophores[i,"name_pock1"],features],
                                      dt[data_diffPdiffL_dt72_pharmacophores[i,"name_pock2"],features])))
}
dist_lig_diffP_diffL

mean(dist_lig_diffP_diffL)
#### distance sans biais en faisant la moyenne ####
#samePsameL
dist_lig_sameP_sameL_mean = rep(0,length(unique(names(dist_lig_sameP_sameL))))
names(dist_lig_sameP_sameL_mean) = unique(names(dist_lig_sameP_sameL))

for (name_lig in unique(names(dist_lig_sameP_sameL))) {
  #print(dist_lig_sameP_sameL[which(names(dist_lig_sameP_sameL) == name_lig)])
  dist_lig_sameP_sameL_mean[name_lig] = mean(dist_lig_sameP_sameL[which(names(dist_lig_sameP_sameL) == name_lig)])
}
boxplot(dist_lig_sameP_sameL_mean)
#diffPsameL
dist_lig_diffP_sameL_mean = rep(0,length(unique(names(dist_lig_diffP_sameL))))
names(dist_lig_diffP_sameL_mean) = unique(names(dist_lig_diffP_sameL))

for (name_lig in unique(names(dist_lig_diffP_sameL))) {
  #print(dist_lig_diffP_sameL[which(names(dist_lig_diffP_sameL) == name_lig)])
  dist_lig_diffP_sameL_mean[name_lig] = mean(dist_lig_diffP_sameL[which(names(dist_lig_diffP_sameL) == name_lig)])
}
boxplot(dist_lig_diffP_sameL_mean)

#
dist_lig_sameP_sameL = dist_lig_sameP_sameL_mean
dist_lig_diffP_sameL = dist_lig_diffP_sameL_mean
dist_lig_diffP_diffL = dist_lig_diffP_diffL
#
### resultat distance pharmacophores ###
#euclidean
load("../results/pharmacophores_results/data/dist_data_intersect72/dist_lig_sameP_sameL_euclidean.Rdata")
load("../results/pharmacophores_results/data/dist_data_intersect72/dist_lig_diffP_sameL_euclidean.Rdata")
load("../results/pharmacophores_results/data/dist_data_intersect72/dist_lig_diffP_diffL_euclidean.Rdata")
#manhattan
load("../results/pharmacophores_results/data/dist_data_intersect72/dist_lig_sameP_sameL_manhattan.Rdata")
load("../results/pharmacophores_results/data/dist_data_intersect72/dist_lig_diffP_sameL_manhattan.Rdata")
load("../results/pharmacophores_results/data/dist_data_intersect72/dist_lig_diffP_diffL_manhattan.Rdata")
#fuzcav ncerisier
load("../results/pharmacophores_results/data/dist_data_intersect72/dist_lig_sameP_sameL_fuzcavPerso.Rdata")
load("../results/pharmacophores_results/data/dist_data_intersect72/dist_lig_diffP_sameL_fuzcavPerso.Rdata")
load("../results/pharmacophores_results/data/dist_data_intersect72/dist_lig_diffP_diffL_fuzcavPerso.Rdata")
#fuzcav vrai
load("../results/pharmacophores_results/data/dist_data_intersect72/dist_lig_sameP_sameL_fuzcaVrai.Rdata")
load("../results/pharmacophores_results/data/dist_data_intersect72/dist_lig_diffP_sameL_fuzcaVrai.Rdata")
load("../results/pharmacophores_results/data/dist_data_intersect72/dist_lig_diffP_diffL_fuzcaVrai.Rdata")
#
names(dist_lig_sameP_sameL) = sapply(strsplit(data_samePsameL_dt72_pharmacophores$name_pock1, "_"), "[", 2)
names(dist_lig_diffP_sameL) = sapply(strsplit(data_diffPsameL_dt72_pharmacophores$name_pock1, "_"), "[", 2)
#
length(which(dist_lig_sameP_sameL == 0))
###ROC CURVE ###
##SCALE
range01 <- function(x,min_scale,max_scale){(x-min_scale)/(max_scale-min_scale)}
#
min_scale = min(c(dist_lig_sameP_sameL,dist_lig_diffP_sameL,dist_lig_diffP_diffL))
max_scale = max(c(dist_lig_sameP_sameL,dist_lig_diffP_sameL,dist_lig_diffP_diffL))
# AUC  diffLdiffP vs sameLsameP
y_true = c(rep(0,length(dist_lig_sameP_sameL)),
           rep(1,length(dist_lig_diffP_diffL)))
y_predict = c(range01(dist_lig_sameP_sameL,min_scale,max_scale),
              range01(dist_lig_diffP_diffL,min_scale,max_scale))
#AUC diffLdiffP vs sameLdiffP
y_true = c(rep(0,length(dist_lig_diffP_sameL)),
           rep(1,length(dist_lig_diffP_diffL)))
y_predict = c(range01(dist_lig_diffP_sameL,min_scale,max_scale),
              range01(dist_lig_diffP_diffL,min_scale,max_scale))
#AUC  diffLdiffP vs (sameLsameP|sameLdiffP)
y_true = c(rep(0,length(dist_lig_sameP_sameL)+
                 length(dist_lig_diffP_sameL)),
           rep(1,length(dist_lig_diffP_diffL)))
y_predict = c(range01(dist_lig_sameP_sameL,min_scale,max_scale),
              range01(dist_lig_diffP_sameL,min_scale,max_scale),
              range01(dist_lig_diffP_diffL,min_scale,max_scale))
#library(ROCR)
dt.pred = prediction(y_predict, y_true_test) #y_predict[-Index]#1-y_predict#y_true_test#y_predict[,2]
dt.perf = performance(dt.pred, "tpr", "fpr")
plot(dt.perf)
dt.auc = performance(dt.pred, "auc")
attr(dt.auc, "y.values")
#
#plot
sm.density.compare(c(range01(dist_lig_diffP_diffL,min_scale,max_scale),
                     range01(dist_lig_sameP_sameL,min_scale,max_scale),
                     range01(dist_lig_diffP_sameL,min_scale,max_scale)),
                   c(rep(1,length(dist_lig_diffP_diffL)),
                     rep(2,length(dist_lig_sameP_sameL)),
                     rep(3,length(dist_lig_diffP_sameL))
                   ),
                   model = "none",xlim=c(0,1)
                   , xlab = "Mean distance bewteen pockets")
sm.density.compare(c(dist_lig_diffP_diffL,
                     dist_lig_sameP_sameL,
                     dist_lig_diffP_sameL),
                   c(rep(1,length(dist_lig_diffP_diffL)),
                     rep(2,length(dist_lig_sameP_sameL)),
                     rep(3,length(dist_lig_diffP_sameL))
                   ),
                   model = "none", xlim = c(0,500),
                   , xlab = "Mean distance bewteen pockets")

### distance euclideanne + distance fuzcav ###
dist_lig_sameP_sameL_desc = range01(dist_lig_sameP_sameL,min_scale,max_scale)
dist_lig_diffP_sameL_desc = range01(dist_lig_diffP_sameL,min_scale,max_scale)
dist_lig_diffP_diffL_desc = range01(dist_lig_diffP_diffL,min_scale,max_scale)
#
dist_lig_sameP_sameL_desc_fuzcav = range01(dist_lig_sameP_sameL,min_scale,max_scale) #1-
dist_lig_diffP_sameL_desc_fuzcav = range01(dist_lig_diffP_sameL,min_scale,max_scale)
dist_lig_diffP_diffL_desc_fuzcav = range01(dist_lig_diffP_diffL,min_scale,max_scale)
#additionner les deux
dist_lig_sameP_sameL = dist_lig_sameP_sameL_desc + dist_lig_sameP_sameL_desc_fuzcav
dist_lig_diffP_sameL = dist_lig_diffP_sameL_desc + dist_lig_diffP_sameL_desc_fuzcav
dist_lig_diffP_diffL = dist_lig_diffP_diffL_desc + dist_lig_diffP_diffL_desc_fuzcav

### linear regression on difference between descriptors ###
dt_pock_samePsameL = dt[data_samePsameL_dt72_pharmacophores[,"name_pock1"],features] - dt[data_samePsameL_dt72_pharmacophores[,"name_pock2"],features]
dt_pock_diffPsameL = dt[data_diffPsameL_dt72_pharmacophores[,"name_pock1"],features] - dt[data_diffPsameL_dt72_pharmacophores[,"name_pock2"],features]
dt_pock_diffPdiffL = dt[data_diffPdiffL_dt72_pharmacophores[,"name_pock1"],features] - dt[data_diffPdiffL_dt72_pharmacophores[,"name_pock2"],features]

y_true = c(rep(1,nrow(dt_pock_samePsameL)),
           rep(1,nrow(dt_pock_diffPsameL)),
           rep(0,nrow(dt_pock_diffPdiffL)))  
dt_pock = rbind(dt_pock_samePsameL,
                rbind(dt_pock_diffPsameL,
                      dt_pock_diffPdiffL))

dt_pock = cbind(abs(dt_pock),y_true)
dt_pock = as.data.frame(dt_pock)
similarity.modele = glm(y_true~.,data = dt_pock, family = "binomial")

attributes(similarity.modele)

summary(similarity.modele)
similarity.modele.step = step(similarity.modele, direction = "both")
summary(similarity.modele.step)

y_predict = predict.glm(similarity.modele.step, newdata=dt_pock[,features],type = "response" )
which(y_predict == 0)
mean(y_predict)
boxplot(y_predict)
#y_predict = range01(y_predict,min(y_predict),max(y_predict))
## PLOT pharmacophhore / descripteur
plot(dist_lig_sameP_sameL_desc, 1-dist_lig_sameP_sameL_desc_fuzcav,xlim = c(0,1),ylim = c(0,1))
cor(dist_lig_sameP_sameL_desc, 1-dist_lig_sameP_sameL_desc_fuzcav)
### check nombre total interaction dataset
names_ligand
table_names_ligand = table(names_ligand)
N_inter = 0
for (i in 1:length(table_names_ligand)) {
  if(table_names_ligand[i] > 1) {
    N_inter = N_inter + (table_names_ligand[i]*table_names_ligand[i])/2 -(table_names_ligand[i]/2)
  }
}
max(table_names_ligand)

#### PREDICTION ####
features = c(descriptors_hydrophobicity,
             descriptors_aromatic,
             descriptors_polarity,
             descriptors_physicochemical,
             descriptors_geometrical)
features = c("SURFACE_HULL","hydrophobic_kyte")
features = descriptors_geometrical
features = c(descriptors_hydrophobicity,
             descriptors_aromatic,
             descriptors_polarity,
             descriptors_physicochemical)
features = features_selection_interet_cor9
colnames(dt[,35:54])
features = c(features,"A","C","E","D","G","F","I","H","K","M","L","N","Q","P","S","R","T","W","V","Y")
#
#data_diffPdiffL_dt72_pharmacophores = read.csv("../data/data_structure_comparaison/data_diffPdiffL_dt72_pharmacophores_50000.csv", colClasses = "character")
data_diffPdiffL_dt72_pharmacophores = read.csv("../data/data_structure_comparaison/data_diffPdiffL_dt72_pharmacophores.csv", colClasses = "character")
#data_prediction
dt_pock_samePsameL = sqrt((dt[data_samePsameL_dt72_pharmacophores[,"name_pock1"],features] - dt[data_samePsameL_dt72_pharmacophores[,"name_pock2"],features])**2)
names(dt_pock_samePsameL) = sapply(strsplit(data_samePsameL_dt72_pharmacophores$name_pock1, "_"), "[", 2)
dt_pock_diffPsameL = sqrt((dt[data_diffPsameL_dt72_pharmacophores[,"name_pock1"],features] - dt[data_diffPsameL_dt72_pharmacophores[,"name_pock2"],features])**2)
names(dt_pock_diffPsameL) = sapply(strsplit(data_diffPsameL_dt72_pharmacophores$name_pock1, "_"), "[", 2)
dt_pock_diffPdiffL = sqrt((dt[data_diffPdiffL_dt72_pharmacophores[,"name_pock1"],features] - dt[data_diffPdiffL_dt72_pharmacophores[,"name_pock2"],features])**2)
#FOR OVERLAP
#withoutNA
dt = na.omit(dt[,features])
#
index_data_samePsameL_dt72_pharmacophores = which(is.element(data_samePsameL_dt72_pharmacophores[,"name_pock1"],rownames(dt)) == TRUE &
                                                    is.element(data_samePsameL_dt72_pharmacophores[,"name_pock2"],rownames(dt)) == TRUE)
#
index_data_diffPsameL_dt72_pharmacophores = which(is.element(data_diffPsameL_dt72_pharmacophores[,"name_pock1"],rownames(dt)) == TRUE &
                                                    is.element(data_diffPsameL_dt72_pharmacophores[,"name_pock2"],rownames(dt)) == TRUE)
#
index_data_diffPdiffL_dt72_pharmacophores = which(is.element(data_diffPdiffL_dt72_pharmacophores[,"name_pock1"],rownames(dt)) == TRUE &
                                                    is.element(data_diffPdiffL_dt72_pharmacophores[,"name_pock2"],rownames(dt)) == TRUE)

dt_pock_samePsameL = sqrt((dt[data_samePsameL_dt72_pharmacophores[index_data_samePsameL_dt72_pharmacophores,"name_pock1"],features] - dt[data_samePsameL_dt72_pharmacophores[index_data_samePsameL_dt72_pharmacophores,"name_pock2"],features])**2)
names(dt_pock_samePsameL) = sapply(strsplit(data_samePsameL_dt72_pharmacophores[index_data_samePsameL_dt72_pharmacophores,]$name_pock1, "_"), "[", 2)
dt_pock_diffPsameL = sqrt((dt[data_diffPsameL_dt72_pharmacophores[index_data_diffPsameL_dt72_pharmacophores,"name_pock1"],features] - dt[data_diffPsameL_dt72_pharmacophores[index_data_diffPsameL_dt72_pharmacophores,"name_pock2"],features])**2)
names(dt_pock_diffPsameL) = sapply(strsplit(data_diffPsameL_dt72_pharmacophores[index_data_diffPsameL_dt72_pharmacophores,]$name_pock1, "_"), "[", 2)
dt_pock_diffPdiffL = sqrt((dt[data_diffPdiffL_dt72_pharmacophores[index_data_diffPdiffL_dt72_pharmacophores,"name_pock1"],features] - dt[data_diffPdiffL_dt72_pharmacophores[index_data_diffPdiffL_dt72_pharmacophores,"name_pock2"],features])**2)

# MEAN
#### distance sans biais en faisant la moyenne ####
#samePsameL
dist_lig_sameP_sameL_mean = NULL
unique_names_dist_lig_sameP_sameL = unique(names(dt_pock_samePsameL))
for (name_lig in 1:length(unique_names_dist_lig_sameP_sameL)) {
  samePsameL_pock = which(names(dt_pock_samePsameL) == unique_names_dist_lig_sameP_sameL[name_lig])
  if(length(samePsameL_pock) > 1) {
    dist_lig_sameP_sameL_mean= rbind(dist_lig_sameP_sameL_mean, apply(dt_pock_samePsameL[samePsameL_pock,],2,mean))
  } else {
    dist_lig_sameP_sameL_mean= rbind(dist_lig_sameP_sameL_mean, dt_pock_samePsameL[samePsameL_pock,])
  }
}
names(dist_lig_sameP_sameL_mean) = unique_names_dist_lig_sameP_sameL
#diffPsameL
dist_lig_diffP_sameL_mean = NULL
unique_names_dist_lig_diffP_sameL = unique(names(dt_pock_diffPsameL))
for (name_lig in 1:length(unique_names_dist_lig_diffP_sameL)) {
  diffPsameL_pock = which(names(dt_pock_diffPsameL) == unique_names_dist_lig_diffP_sameL[name_lig])
  if(length(diffPsameL_pock) > 1) {
    dist_lig_diffP_sameL_mean= rbind(dist_lig_diffP_sameL_mean, apply(dt_pock_diffPsameL[diffPsameL_pock,],2,mean))
  } else {
    dist_lig_diffP_sameL_mean= rbind(dist_lig_diffP_sameL_mean, dt_pock_diffPsameL[diffPsameL_pock,])
  }
}
names(dist_lig_diffP_sameL_mean) = unique_names_dist_lig_diffP_sameL
#
dt_pock_samePsameL = dist_lig_sameP_sameL_mean
dt_pock_diffPsameL = dist_lig_diffP_sameL_mean
#data preparation
dt_predict = rbind(dt_pock_samePsameL,
             rbind(dt_pock_diffPsameL,
                   dt_pock_diffPdiffL))
y_true = c(rep(1,nrow(dt_pock_samePsameL)),
           rep(1,nrow(dt_pock_diffPsameL)),
           rep(0,nrow(dt_pock_diffPdiffL))) 
dt_predict = cbind(abs(dt_predict),y_true)
dt_predict = as.data.frame(dt_predict)
dt_predict$y_true = as.factor(dt_predict$y_true)
##ACP
dt_predict.apc = PCA(dt_predict, scale.unit = T, quali.sup=49)
dt_predict.apc$eig
dt_predict.apc$ind$contrib

sort(dt_predict.apc$var$contrib[,"Dim.1"])

barplot(dt_predict.apc$eig[,2])
plot(dt_predict.apc)
#reprÃ©ssentation similar or not
plot(dt_predict.apc,habillage = 49, col.hab = c("green","blue"), label="none")
plot(dt_predict.apc,choix = "ind", col.ind = dt_predict$y_true)

### data preparation ###
set.seed(84)
Index = sample(1:nrow(dt_predict),size = nrow(dt_predict)*2/3 )
dt_predict_app = dt_predict[Index,]
dt_predict_test = dt_predict[-Index,]
dt_predict_test_samePsameL = dt_predict[-unique(c(Index,nrow(dt_pock_samePsameL):(nrow(dt_pock_samePsameL)+nrow(dt_pock_diffPsameL)))),]
dt_predict_test_diffPsameL = dt_predict[-unique(c(Index,1:nrow(dt_pock_samePsameL))),]
#
y_true_app = y_true[Index]
y_true_test = y_true[-Index]
y_true_test_samePsameL = y_true[-unique(c(Index,nrow(dt_pock_samePsameL):(nrow(dt_pock_samePsameL)+nrow(dt_pock_diffPsameL))))]
y_true_test_diffPsameL = y_true[-unique(c(Index,1:nrow(dt_pock_samePsameL)))]

#### Linear regression#### 
dt_predict.glm = glm(y_true~.,data = dt_predict_app, family = "binomial")
attributes(dt_predict.glm)
summary(dt_predict.glm)

dt_predict.glm.step = step(dt_predict.glm, direction = "both")
summary(dt_predict.glm.step)
#save(dt_predict.glm.step, file = "../results/kmeans_results_reglog/model.glm.step.Rdata")
pock_test = sqrt((dt["101M_HEM_A_1",features] - dt["102M_HEM_A_1",features])**2)
pock_test = as.data.frame(rbind(pock_test,pock_test))

predict.glm(dt_predict.glm.step, newdata=pock_test[,features], type = "response" )
exp(sim)/(1+exp(sim))

similarity_regol = function(dt_pock1,dt_pock2) {
  Intercept = 6.61543
  hydrophobic_kyte = -0.52657#vec1[56]-vec2[56]
  p_aromatic_residues = -0.26214#vec1[69]-vec2[69]
  p_charged_residues = -0.42272#vec1[6]-vec2[6]
  p_negative_residues = 0.33210#vec1[23]-vec2[23]
  p_positive_residues = 0.59123#vec1[10]-vec2[10]
  charge = -0.73920#vec1[24]-vec2[24]
  p_Nlys_atom = -0.26117#vec1[15]-vec2[15]
  p_Ocoo_atom = -0.16479#vec1[67]-vec2[67]
  p_small_residues = -0.16480#vec1[7]-vec2[7]
  p_C_atom = 0.19381#vec1[64]-vec2[64]
  p_nitrogen_atom = -0.36986#vec1[9]-vec2[9]
  RADIUS_HULL = -3.33198#vec1[31]-vec2[31]
  SURFACE_HULL = -3.57145#vec1[73]-vec2[73] 
  DIAMETER_HULL = 1.10726#vec1[49]-vec2[49]
  VOLUME_HULL = 2.46588#vec1[28]-vec2[28]
  RADIUS_CYLINDER = 1.07161#vec1[3]-vec2[3]
  C_RESIDUES = -0.65951#vec1[63]-vec2[63]
  PSI = -0.18916#vec1[58]-vec2[58]
  PCI = 0.16018#vec1[12]-vec2[12]
  INERTIA_2 = -0.32137#vec1[18]-vec2[18]
  INERTIA_3 = -0.15445#vec1[70]-vec2[70]
  A = -0.16218#vec1[32]-vec2[32]
  C = -0.34735#vec1[34]-vec2[34]
  D = -0.15076#vec1[36]-vec2[36]
  G = -0.43486#vec1[37]-vec2[37]
  F = -0.18449#vec1[38]-vec2[38]
  I = -0.30728#vec1[39]-vec2[39]
  H = -0.27257#vec1[40]-vec2[40]
  K = -0.27611#vec1[41]-vec2[41]
  M = -0.27610#vec1[42]-vec2[42]
  L = -0.34013#vec1[43]-vec2[43]
  Q = -0.11223#vec1[45]-vec2[45]
  P = -0.11730#vec1[46]-vec2[46]
  S = -0.30554#vec1[47]-vec2[47]
  R = -0.22685#vec1[48]-vec2[48]
  T = -0.30320#vec1[50]-vec2[50]
  W = -0.15142#vec1[51]-vec2[51]
  Y = -0.11486#vec1[53]-vec2[53]
  sim = (Intercept + 
           hydrophobic_kyte * sqrt((dt[dt_pock1,"hydrophobic_kyte"] - dt[dt_pock2,"hydrophobic_kyte"])**2) + 
           p_aromatic_residues * sqrt((dt[dt_pock1,"p_aromatic_residues"] - dt[dt_pock2,"p_aromatic_residues"])**2) +
           p_charged_residues * sqrt((dt[dt_pock1,"p_charged_residues"] - dt[dt_pock2,"p_charged_residues"])**2) +
           p_negative_residues * sqrt((dt[dt_pock1,"p_negative_residues"] - dt[dt_pock2,"p_negative_residues"])**2) +
           p_positive_residues * sqrt((dt[dt_pock1,"p_positive_residues"] - dt[dt_pock2,"p_positive_residues"])**2) +
           charge * sqrt((dt[dt_pock1,"charge"] - dt[dt_pock2,"charge"])**2) +
           p_Nlys_atom *sqrt((dt[dt_pock1,"p_Nlys_atom"] - dt[dt_pock2,"p_Nlys_atom"])**2) +
           p_Ocoo_atom *sqrt((dt[dt_pock1,"p_Ocoo_atom"] - dt[dt_pock2,"p_Ocoo_atom"])**2) +
           p_small_residues *sqrt((dt[dt_pock1,"p_small_residues"] - dt[dt_pock2,"p_small_residues"])**2) +
           p_C_atom *sqrt((dt[dt_pock1,"p_C_atom"] - dt[dt_pock2,"p_C_atom"])**2) +
           p_nitrogen_atom *sqrt((dt[dt_pock1,"p_nitrogen_atom"] - dt[dt_pock2,"p_nitrogen_atom"])**2) +
           RADIUS_HULL *sqrt((dt[dt_pock1,"RADIUS_HULL"] - dt[dt_pock2,"RADIUS_HULL"])**2) +
           SURFACE_HULL *sqrt((dt[dt_pock1,"SURFACE_HULL"] - dt[dt_pock2,"SURFACE_HULL"])**2) +
           DIAMETER_HULL *sqrt((dt[dt_pock1,"DIAMETER_HULL"] - dt[dt_pock2,"DIAMETER_HULL"])**2) +
           VOLUME_HULL *sqrt((dt[dt_pock1,"VOLUME_HULL"] - dt[dt_pock2,"VOLUME_HULL"])**2) +
           RADIUS_CYLINDER *sqrt((dt[dt_pock1,"RADIUS_CYLINDER"] - dt[dt_pock2,"RADIUS_CYLINDER"])**2) +
           C_RESIDUES *sqrt((dt[dt_pock1,"C_RESIDUES"] - dt[dt_pock2,"C_RESIDUES"])**2) +
           PSI *sqrt((dt[dt_pock1,"PSI"] - dt[dt_pock2,"PSI"])**2) +
           PCI *sqrt((dt[dt_pock1,"PCI"] - dt[dt_pock2,"PCI"])**2) +
           INERTIA_2 *sqrt((dt[dt_pock1,"INERTIA_2"] - dt[dt_pock2,"INERTIA_2"])**2) +
           INERTIA_3 *sqrt((dt[dt_pock1,"INERTIA_3"] - dt[dt_pock2,"INERTIA_3"])**2) +
           A *sqrt((dt[dt_pock1,"A"] - dt[dt_pock2,"A"])**2) +
           C *sqrt((dt[dt_pock1,"C"] - dt[dt_pock2,"C"])**2) +
           D *sqrt((dt[dt_pock1,"D"] - dt[dt_pock2,"D"])**2) +
           G *sqrt((dt[dt_pock1,"G"] - dt[dt_pock2,"G"])**2) +
           F *sqrt((dt[dt_pock1,"F"] - dt[dt_pock2,"F"])**2) +
           I *sqrt((dt[dt_pock1,"I"] - dt[dt_pock2,"I"])**2) +
           H *sqrt((dt[dt_pock1,"H"] - dt[dt_pock2,"H"])**2) +
           K *sqrt((dt[dt_pock1,"K"] - dt[dt_pock2,"K"])**2) +
           M *sqrt((dt[dt_pock1,"M"] - dt[dt_pock2,"M"])**2) +
           L *sqrt((dt[dt_pock1,"L"] - dt[dt_pock2,"L"])**2) +
           Q *sqrt((dt[dt_pock1,"Q"] - dt[dt_pock2,"Q"])**2) +
           P *sqrt((dt[dt_pock1,"P"] - dt[dt_pock2,"P"])**2) +
           S *sqrt((dt[dt_pock1,"S"] - dt[dt_pock2,"S"])**2) +
           R *sqrt((dt[dt_pock1,"R"] - dt[dt_pock2,"R"])**2) +
           T *sqrt((dt[dt_pock1,"T"] - dt[dt_pock2,"T"])**2) +
           W *sqrt((dt[dt_pock1,"W"] - dt[dt_pock2,"W"])**2) +
           Y *sqrt((dt[dt_pock1,"Y"] - dt[dt_pock2,"Y"])**2) 
  )
  return(1-exp(sim)/(1+exp(sim)))
}
#
y_predict = predict.glm(dt_predict.glm.step, newdata=dt_predict_app[,features],type = "response" )
y_predict = predict.glm(dt_predict.glm.step, newdata=dt_predict_test[,features],type = "response" )
y_predict = predict.glm(dt_predict.glm.step, newdata=dt_predict_test_samePsameL[,features],type = "response" )
y_predict = predict.glm(dt_predict.glm.step, newdata=dt_predict_test_diffPsameL[,features],type = "response" )
#
perf_auc(y_predict, y_true_app)
perf_auc(y_predict, y_true_test)
perf_auc(y_predict, y_true_test_samePsameL)
perf_auc(y_predict, y_true_test_diffPsameL)
#
y_predict[which(y_predict > 0.5)] = 1
y_predict[which(y_predict < 0.5)] = 0
#MCC
#library(mltools)
mcc(preds = y_predict, actuals = y_true_test_diffPsameL)
#
glm.table <- table(y_predict, y_true_test)
Se = glm.table[2,2] / (glm.table[2,2] + glm.table[1,2])
Sp = glm.table[1,1] / (glm.table[1,1] + glm.table[2,1])
Se
Sp
#
nrow(dt_predict)
### log + pharmacophores fuzcav
y_predict_desc_test = y_predict
y_predict_desc_samePsameL = y_predict
y_predict_desc_diffPsameL = y_predict
#
y_predict_fuzcav_test = y_predict[-Index]
y_predict_fuzcav_samePsameL = y_predict[-unique(c(Index,nrow(dt_pock_samePsameL):(nrow(dt_pock_samePsameL)+nrow(dt_pock_diffPsameL))))]
y_predict_fuzcav_diffPsameL = y_predict[-unique(c(Index,1:nrow(dt_pock_samePsameL)))]
#
y_predict = y_predict_desc_test + y_predict_fuzcav_test
y_predict = y_predict_desc_samePsameL + y_predict_fuzcav_samePsameL
y_predict = y_predict_desc_diffPsameL + y_predict_fuzcav_diffPsameL
#
plot(y_predict_desc_test,y_predict_fuzcav_test)
plot(y_predict_desc_samePsameL,y_predict_fuzcav_samePsameL)
plot(y_predict_desc_diffPsameL,y_predict_fuzcav_diffPsameL)

#### CROSS VALIDATION REGRESSION LOGISTIQUE ###
fct.CV = function(Mat,nbG){
  ###Mat : matrice a echantillonner
  ##nbG : nombre de groupe a faire
  Ech.sample=sample(dim(Mat)[1])
  vect=ceiling(seq(1,dim(Mat)[1],length=(nbG+1)))
  
  echlist=as.list(rep(0,nbG))
  
  for(i in 1:length(echlist)){
    if (i!=nbG){ echlist[[i]]=Ech.sample[vect[i]:((vect[i+1])-1)] }
    if (i==nbG){ echlist[[i]]=Ech.sample[vect[i]:length(Ech.sample)] }
  }
  return (echlist)
}
CV.GLM= function(X,nm.cl,iteration, nb.grp){
  ###iteration : nombre de fois que la CV est rÃ©aliser	
  ###nb.grp :  nbr de fold pour la CV
  ###X : matrice de descripteurs 
  ###nm.cl :  vecteur de classe Ã predire
  ###Pouc chq CV, on va calculer AUC sur l'ech a predire.
  ###Ces taux seront stockes dans le vecteur res.tx
  AUC.tx <- NULL
  for(it in 1:iteration){
    perf_ssjeu <- rep(NA, length = nb.grp )
    ech.list <- fct.CV(Mat = X, nbG = nb.grp)
    for(i in 1:nb.grp){
      #choisi les 3 echantillons pour l apprentissage
      num.ech <- seq(1,nb.grp)[-i]
      list.mol <- NULL
      for (ind in num.ech){
        list.mol <- c(list.mol, ech.list[[ind]])     ###training set
      }
      X1 <- X[list.mol,]
      Y <- as.factor(nm.cl[list.mol])
      ###creation de la matrice a predire
      Matpred <- cbind(X1,Y)
      colnames(Matpred) <- c(colnames(X1), "pred")
      #calcul du modele
      modele <- glm(pred ~ ., data = Matpred,family = "binomial")
      #modele = step(modele, direction = "both")
      ###creation of validation set
      new.data <- X[ech.list[[i]],]
      Ypred <- predict(modele,newdata=new.data, type = "response")
      Yobs <- nm.cl[ech.list[[i]]]
      dt.pred = prediction(Ypred, Yobs)
      dt.auc = performance(dt.pred, "auc")
      perf_ssjeu[i] <- attr(dt.auc, "y.values")[[1]]
    }   ###fin de la boucle de cross-validation
    #compute well-predicted rate for the CV test set
    AUC.tx = c(AUC.tx,mean(perf_ssjeu))
  }
  return(AUC.tx)
}
AUC.taux = CV.GLM(dt_predict_app[,-ncol(dt_predict_app)],dt_predict_app$y_true,200,3) #grand nombre de iterations en cross validation 500,1000...
t.test(AUC.taux)
summary(dt_predict_app)
sd(AUC.taux)

#### CART ####
library(MASS)
library(rpart)
library(rpart.plot)
#MODELE
dt_predict.cart = rpart(y_true~., data = dt_predict_app,
                        method = "class",
                        parms = list(split = 'information'), 
                        minsplit = 2, 
                        minbucket = 1,
                        cp = -1)
rpart.plot(dt_predict.cart)

plotcp(dt_predict.cart)
print(dt_predict.cart$cptable[which.min(dt_predict.cart$cptable[,4]),1])
print(dt_predict.cart)
summary(dt_predict.cart)
prp(dt_predict.cart, extra = 2)
#otp tree
dt_predict.cart.Tree_Opt <- prune(dt_predict.cart,cp=dt_predict.cart$cptable[which.min(dt_predict.cart$cptable[,4]),1])
#dt_predict.cart.Tree_Opt <- prune(dt_predict.cart,cp = 0.0072)

plotcp(dt_predict.cart.Tree_Opt)
rpart.plot(dt_predict.cart.Tree_Opt)

y_predict = predict(dt_predict.cart.Tree_Opt, newdata=dt_predict_test,type = "prob")
#
y_predict = predict(dt_predict.cart.Tree_Opt, newdata=dt_predict_app[,features],type = "prob" )
y_predict = predict(dt_predict.cart.Tree_Opt, newdata=dt_predict_test[,features],type = "prob" )
y_predict = predict(dt_predict.cart.Tree_Opt, newdata=dt_predict_test_samePsameL[,features],type = "prob" )
y_predict = predict(dt_predict.cart.Tree_Opt, newdata=dt_predict_test_diffPsameL[,features],type = "prob" )
#
perf_auc(y_predict[,2], y_true_app)
perf_auc(y_predict[,2], y_true_test)
perf_auc(y_predict[,2], y_true_test_samePsameL)
perf_auc(y_predict[,2], y_true_test_diffPsameL)
#
y_predict[which(y_predict[,2] > 0.5),2] = 1
y_predict[which(y_predict[,2] < 0.5),2] = 0
glm.table <- table(y_predict[,2], y_true_test_diffPsameL)
Se = glm.table[2,2] / (glm.table[2,2] + glm.table[1,2])
Sp = glm.table[1,1] / (glm.table[1,1] + glm.table[2,1])
Se
Sp
#
### CART + pharmacophores fuzcav
y_predict_desc_test = y_predict[,2]
y_predict_desc_samePsameL = y_predict[,2]
y_predict_desc_diffPsameL = y_predict[,2]
#
y_predict_fuzcav_test = y_predict[-Index]
y_predict_fuzcav_samePsameL = y_predict[-unique(c(Index,nrow(dt_pock_samePsameL):(nrow(dt_pock_samePsameL)+nrow(dt_pock_diffPsameL))))]
y_predict_fuzcav_diffPsameL = y_predict[-unique(c(Index,1:nrow(dt_pock_samePsameL)))]
#
y_predict = y_predict_desc_test + y_predict_fuzcav_test
y_predict = y_predict_desc_samePsameL + y_predict_fuzcav_samePsameL
y_predict = y_predict_desc_diffPsameL + y_predict_fuzcav_diffPsameL
#
plot(y_predict_desc_test,y_predict_fuzcav_test)
plot(y_predict_desc_samePsameL,y_predict_fuzcav_samePsameL)
plot(y_predict_desc_diffPsameL,y_predict_fuzcav_diffPsameL)
###Comparaison  avec dist normale
#### CROSS VALIDATION CART ###
CV.CART = function(X,nm.cl,iteration, nb.grp){
  ###iteration : nombre de fois que la CV est rÃ©aliser	
  ###nb.grp :  nbr de fold pour la CV
  ###X : matrice de descripteurs 
  ###nm.cl :  vecteur de classe Ã predire
  ###Pouc chq CV, on va calculer AUC sur l'ech a predire.
  ###Ces taux seront stockes dans le vecteur res.tx
  AUC.tx <- NULL
  for(it in 1:iteration){
    perf_ssjeu <- rep(NA, length = nb.grp )
    ech.list <- fct.CV(Mat = X, nbG = nb.grp)
    for(i in 1:nb.grp){
      #choisi les 3 echantillons pour l apprentissage
      num.ech <- seq(1,nb.grp)[-i]
      list.mol <- NULL
      for (ind in num.ech){
        list.mol <- c(list.mol, ech.list[[ind]])     ###training set
      }
      X1 <- X[list.mol,]
      Y <- as.factor(nm.cl[list.mol])
      ###creation de la matrice a predire
      Matpred <- cbind(X1,Y)
      colnames(Matpred) <- c(colnames(X1), "pred")
      #calcul du modele
      modele <- rpart(y_true~., data = dt_predict_app,
                      method = "class",
                      parms = list(split = 'information'), 
                      minsplit = 2, 
                      minbucket = 1,
                      cp = -1)
      modele <- prune(modele,cp=modele$cptable[which.min(modele$cptable[,4]),1])
      
      #modele = step(modele, direction = "both")
      ###creation of validation set
      new.data <- X[ech.list[[i]],]
      Ypred <- predict(modele,newdata=new.data, type = "prob")
      Yobs <- nm.cl[ech.list[[i]]]
      dt.pred = prediction(Ypred[,2], Yobs)
      dt.auc = performance(dt.pred, "auc")
      perf_ssjeu[i] <- attr(dt.auc, "y.values")[[1]]
    }   ###fin de la boucle de cross-validation
    #compute well-predicted rate for the CV test set
    AUC.tx = c(AUC.tx,mean(perf_ssjeu))
  }
  return(AUC.tx)
}
AUC.taux = CV.CART(dt_predict_app[,-ncol(dt_predict_app)],dt_predict_app$y_true,20,3) #grand nombre de iterations en cross validation 500,1000...
t.test(AUC.taux)
summary(dt_predict_app)
sd(AUC.taux)

#### LDA ####
dt_predict.lda = lda(y_true~.,data = dt_predict_app)
dt_predict.lda$prior
attributes(dt_predict.lda)
summary(dt_predict.lda)

y_predict = predict(dt_predict.lda, newdata=dt_predict_test[,features],type = "response" )

F = as.matrix(dt_predict_app[,-ncol(dt_predict_app)]) %*% as.matrix(dt_predict.lda$scaling)
#library(knitr)
kable(cor(F, dt_predict_app[,-ncol(dt_predict_app)]), digits = 2)
nrow(dt_predict)
#
y_predict = predict(dt_predict.lda, newdata=dt_predict_app[,features],type = "response" )
y_predict = predict(dt_predict.lda, newdata=dt_predict_test[,features],type = "response" )
y_predict = predict(dt_predict.lda, newdata=dt_predict_test_samePsameL[,features],type = "response" )
y_predict = predict(dt_predict.lda, newdata=dt_predict_test_diffPsameL[,features],type = "response" )
#
perf_auc(y_predict$x, y_true_app)
perf_auc(y_predict$x, y_true_test)
perf_auc(y_predict$x, y_true_test_samePsameL)
perf_auc(y_predict$x, y_true_test_diffPsameL)
#
y_predict$x[which(y_predict$x > 0.5)] = 1
y_predict$x[which(y_predict$x < 0.5)] = 0
glm.table <- table(y_predict$x, y_true_test_diffPsameL)
Se = glm.table[2,2] / (glm.table[2,2] + glm.table[1,2])
Sp = glm.table[1,1] / (glm.table[1,1] + glm.table[2,1])
Se
Sp
#
### LDA + pharmacophores fuzcav
y_predict_desc_test = y_predict$x
y_predict_desc_samePsameL = y_predict$x
y_predict_desc_diffPsameL = y_predict$x
#
y_predict_fuzcav_test = y_predict[-Index]
y_predict_fuzcav_samePsameL = y_predict[-unique(c(Index,nrow(dt_pock_samePsameL):(nrow(dt_pock_samePsameL)+nrow(dt_pock_diffPsameL))))]
y_predict_fuzcav_diffPsameL = y_predict[-unique(c(Index,1:nrow(dt_pock_samePsameL)))]
#
y_predict = y_predict_desc_test + y_predict_fuzcav_test
y_predict = y_predict_desc_samePsameL + y_predict_fuzcav_samePsameL
y_predict = y_predict_desc_diffPsameL + y_predict_fuzcav_diffPsameL
#
plot(y_predict_desc_test,y_predict_fuzcav_test)
plot(y_predict_desc_samePsameL,y_predict_fuzcav_samePsameL)
plot(y_predict_desc_diffPsameL,y_predict_fuzcav_diffPsameL)
#### Random Forest ####
library(randomForest)
dt_predict.rf = randomForest(y_true~., data = dt_predict_app, ntree = 200, mtry = 3, importance = T)

plot(dt_predict.rf)
dt_predict.rf
attributes(dt_predict.rf)

y_predict = predict(dt_predict.rf, newdata=dt_predict_test[,features],type = "response" )
table(y_predict,y_true_test)
#
y_predict = predict(dt_predict.rf, newdata=dt_predict_app[,features],type = "prob" )
y_predict = predict(dt_predict.rf, newdata=dt_predict_test[,features],type = "prob" )
y_predict = predict(dt_predict.rf, newdata=dt_predict_test_samePsameL[,features],type = "prob" )
y_predict = predict(dt_predict.rf, newdata=dt_predict_test_diffPsameL[,features],type = "prob" )
#
perf_auc(y_predict[,2], y_true_app)
perf_auc(y_predict[,2], y_true_test)
perf_auc(y_predict[,2], y_true_test_samePsameL)
perf_auc(y_predict[,2], y_true_test_diffPsameL)
#
y_predict[which(y_predict[,2] > 0.5),2] = 1
y_predict[which(y_predict[,2] <= 0.5),2] = 0
glm.table <- table(y_predict[,2], y_true_test_diffPsameL)
Se = glm.table[2,2] / (glm.table[2,2] + glm.table[1,2])
Sp = glm.table[1,1] / (glm.table[1,1] + glm.table[2,1])
Se
Sp
### rf + pharmacophores fuzcav
y_predict_desc_test = y_predict[,2]
y_predict_desc_samePsameL = y_predict[,2]
y_predict_desc_diffPsameL = y_predict[,2]
#
y_predict_fuzcav_test = y_predict[-Index]
y_predict_fuzcav_samePsameL = y_predict[-unique(c(Index,nrow(dt_pock_samePsameL):(nrow(dt_pock_samePsameL)+nrow(dt_pock_diffPsameL))))]
y_predict_fuzcav_diffPsameL = y_predict[-unique(c(Index,1:nrow(dt_pock_samePsameL)))]
#
y_predict = y_predict_desc_test + y_predict_fuzcav_test
y_predict = y_predict_desc_samePsameL + y_predict_fuzcav_samePsameL
y_predict = y_predict_desc_diffPsameL + y_predict_fuzcav_diffPsameL
#
plot(y_predict_desc_test,y_predict_fuzcav_test)
plot(y_predict_desc_samePsameL,y_predict_fuzcav_samePsameL)
plot(y_predict_desc_diffPsameL,y_predict_fuzcav_diffPsameL)
#
nrow(dt_predict_test)
### Performances ###
perf_auc = function(y_predict, y_true_test) {
  dt.pred = prediction(y_predict, y_true_test) #y_predict[-Index]#1-y_predict#y_true_test#y_predict[,2]
  dt.auc = performance(dt.pred, "auc")
  print("AUC")
  print(attr(dt.auc, "y.values"))
}

##### DENSITY CURVE #####

#reg log
dist_lig_sameP_sameL = predict(dt_predict.glm.step, newdata=dt_predict[1:nrow(dt_pock_samePsameL),features],type = "response" )
dist_lig_diffP_sameL = predict(dt_predict.glm.step, newdata=dt_predict[nrow(dt_pock_samePsameL):(nrow(dt_pock_samePsameL)+nrow(dt_pock_diffPsameL)),features],type = "response" )
dist_lig_diffP_diffL = predict(dt_predict.glm.step, newdata=dt_predict[(nrow(dt_pock_samePsameL)+nrow(dt_pock_diffPsameL)):nrow(dt_predict),features],type = "response" )
#reg CART
dist_lig_sameP_sameL = predict(dt_predict.cart.Tree_Opt, newdata=dt_predict[1:nrow(dt_pock_samePsameL),features],type = "prob" )[,2]
dist_lig_diffP_sameL = predict(dt_predict.cart.Tree_Opt, newdata=dt_predict[nrow(dt_pock_samePsameL):(nrow(dt_pock_samePsameL)+nrow(dt_pock_diffPsameL)),features],type = "prob" )[,2]
dist_lig_diffP_diffL = predict(dt_predict.cart.Tree_Opt, newdata=dt_predict[(nrow(dt_pock_samePsameL)+nrow(dt_pock_diffPsameL)):nrow(dt_predict),features],type = "prob" )[,2]
#reg LDA
dist_lig_sameP_sameL = predict(dt_predict.lda, newdata=dt_predict[1:nrow(dt_pock_samePsameL),features],type = "prob" )$x
dist_lig_diffP_sameL = predict(dt_predict.lda, newdata=dt_predict[nrow(dt_pock_samePsameL):(nrow(dt_pock_samePsameL)+nrow(dt_pock_diffPsameL)),features],type = "prob" )$x
dist_lig_diffP_diffL = predict(dt_predict.lda, newdata=dt_predict[(nrow(dt_pock_samePsameL)+nrow(dt_pock_diffPsameL)):nrow(dt_predict),features],type = "prob" )$x
#reg RF
dist_lig_sameP_sameL = predict(dt_predict.rf, newdata=dt_predict[1:nrow(dt_pock_samePsameL),features],type = "prob" )[,2]
dist_lig_diffP_sameL = predict(dt_predict.rf, newdata=dt_predict[nrow(dt_pock_samePsameL):(nrow(dt_pock_samePsameL)+nrow(dt_pock_diffPsameL)),features],type = "prob" )[,2]
dist_lig_diffP_diffL = predict(dt_predict.rf, newdata=dt_predict[(nrow(dt_pock_samePsameL)+nrow(dt_pock_diffPsameL)):nrow(dt_predict),features],type = "prob" )[,2]

#+pharmacophores
dist_lig_sameP_sameL = dist_lig_sameP_sameL + y_predict[1:nrow(dt_pock_samePsameL)]
dist_lig_diffP_sameL = dist_lig_diffP_sameL + y_predict[nrow(dt_pock_samePsameL):(nrow(dt_pock_samePsameL)+nrow(dt_pock_diffPsameL))]
dist_lig_diffP_diffL = dist_lig_diffP_diffL + y_predict[(nrow(dt_pock_samePsameL)+nrow(dt_pock_diffPsameL)):nrow(dt_predict)]
##plot
sm.density.compare(c(dist_lig_diffP_diffL,
                     dist_lig_sameP_sameL,
                     dist_lig_diffP_sameL),
                   c(rep(1,length(dist_lig_diffP_diffL)),
                     rep(2,length(dist_lig_sameP_sameL)),
                     rep(3,length(dist_lig_diffP_sameL))
                   ),
                   model = "none"
                   , xlab = "Mean distance bewteen pockets")
plot(dist_lig_sameP_sameL,y_predict[1:nrow(dt_pock_samePsameL)],
     ylim = c(0,1))
plot(dist_lig_diffP_sameL,y_predict[nrow(dt_pock_samePsameL):(nrow(dt_pock_samePsameL)+nrow(dt_pock_diffPsameL))],
     ylim = c(0,1))
plot(dist_lig_diffP_diffL,y_predict[(nrow(dt_pock_samePsameL)+nrow(dt_pock_diffPsameL)):nrow(dt_predict)],
     ylim = c(0,1))

###POCHE à RETIRER ###
data_samePsameL_dt72_pharmacophores
data_diffPsameL_dt72_pharmacophores
data_diffPdiffL_dt72_pharmacophores

#50
set.seed(83)
pocket_samePsameL_50 = data_samePsameL_dt72_pharmacophores[sample(nrow(data_samePsameL_dt72_pharmacophores),
                                                      size = 50),2]
length(pocket_samePsameL_50)
pocket_diffPsameL_50 = data_diffPsameL_dt72_pharmacophores[sample(nrow(data_diffPsameL_dt72_pharmacophores),
                                                                 size = 50),2]
length(pocket_diffPsameL_50)
pocket_diffPdiffL_50 = data_diffPdiffL_dt72_pharmacophores[sample(nrow(data_diffPdiffL_dt72_pharmacophores),
                                                                 size = 50),2]
length(pocket_diffPdiffL_50)

save(pocket_samePsameL_50, file = "../results/kmedoids_results_reglog/pocket_samePsameL_50.Rdata")
save(pocket_diffPsameL_50, file = "../results/kmedoids_results_reglog/pocket_diffPsameL_50.Rdata")
save(pocket_diffPdiffL_50, file = "../results/kmedoids_results_reglog/pocket_diffPdiffL_50.Rdata")
##dt-50pockets.
dt = dt[setdiff(rownames(dt),c(pocket_samePsameL_50,pocket_diffPsameL_50,pocket_diffPdiffL_50)),]
nrow(dt)
write.csv(dt[setdiff(rownames(dt),c(pocket_samePsameL_50,pocket_diffPsameL_50,pocket_diffPdiffL_50)),],
          file = "../data/dt_72clean_overlap-50.csv")

### REGLOG PLOT Fscore Precision Recall ###
y_predict = predict.glm(dt_predict.glm.step, newdata=dt_predict_test[,features],type = "response" )

perf_precision = NULL
perf_recall = NULL
perf_F1 = NULL
for (i in seq(0,1,0.05)) {
  y_predict = predict.glm(dt_predict.glm.step, newdata=dt_predict_test[,features],type = "response" )
  y_predict[which(y_predict > i)] = 1
  y_predict[which(y_predict <= i)] = 0
  glm.table <- table(factor(y_predict,levels = 0:1), factor(y_true_test,levels = 0:1))

  perf_precision = c(perf_precision,glm.table[2,2] / (glm.table[2,2] + glm.table[2,1]))
  perf_recall = c(perf_recall,glm.table[2,2] / (glm.table[2,2] + glm.table[1,2]))
  perf_F1 = c(perf_F1, (2*glm.table[2,2])/(2*glm.table[2,2]+glm.table[1,2]+glm.table[2,1]))
}

plot(seq(0,1,0.05),perf_precision, type = 'l',col = 2,ylim=c(0,1))
points(seq(0,1,0.05),perf_recall, type = 'l', col = 3)
points(seq(0,1,0.05),perf_F1, type = 'l', col = 4)
legend("bottom", legend = c("precision","recall","score F1"), col = c(2,3,4), lty=1)
