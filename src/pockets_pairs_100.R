### paires Kmeans test ###
#pocket_diffPsameL_50
load("../data/data_structure_comparaison/pocket_pos_100.Rdata")
#pocket_neg_100
load("../data/data_structure_comparaison/pocket_neg_100.Rdata")

#same
dist_lig_sameP_sameL = function(dt, names_ligand, mat_sL, pocket_pos_100) {
  dist_ligs = rep(0,length(which(table(names_ligand) > 1)))
  names(dist_ligs) = names(which(table(names_ligand) > 1))
  flag = 1
  dist_features_pock_all = NULL
  dist_features = matrix(data = 0,
                         nrow = length(which(table(names_ligand) > 1)),
                         ncol = ncol(dt))
  dist_features = as.data.frame(dist_features)
  rownames(dist_features) = names(which(table(names_ligand) > 1))
  data_samePsameL = NULL
  rnames = as.character(mat_sL$names_pockets)
  for (pock in pocket_pos_100){#length(names_ligand_unique)
    print(pock)
    lig_name = mat_sL[which(mat_sL$names_pockets == pock),2]
    index = as.character(mat_sL$names_pockets[which(mat_sL[,2] == lig_name)])
    if(length(index) > 5) {
      index = sample(index,2)
    }
    dist_pock = NULL
    dist_features_pock = NULL
    for (j in index) {
      if(pock != j) {
        if(sapply(strsplit(pock, "_"), "[", 1) != sapply(strsplit(j, "_"), "[", 1)) { #CHANGE TO != FOR DIFFERENT PROTEINS
          data_samePsameL = rbind(data_samePsameL, c(pock,j)) 
          #dist_pock = c(dist_pock, dist(rbind(dt[index[j],],
          #                                    dt[index[k],])))
          #dist_features_pock = rbind(dist_features_pock,  apply(dt[c(index[j],index[k]),],2,dist))
        }
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

data_samePsameL = dist_lig_sameP_sameL(dt,names_ligand, mat_sL, pocket_pos_100)
nrow(data_samePsameL)

#write.csv(data_samePsameL, "../data/data_structure_comparaison/data_samePsameL_100_pos.csv")

## diff
dist_lig_diffP_diffL = function(dt, names_ligand, mat_sL, pocket_pos_100) {
  dist_ligs = rep(0,length(which(table(names_ligand) > 1)))
  names(dist_ligs) = names(which(table(names_ligand) > 1))
  flag = 1
  dist_features_pock_all = NULL
  dist_features = matrix(data = 0,
                         nrow = length(which(table(names_ligand) > 1)),
                         ncol = ncol(dt))
  dist_features = as.data.frame(dist_features)
  rownames(dist_features) = names(which(table(names_ligand) > 1))
  data_samePsameL = NULL
  rnames = as.character(mat_sL$names_pockets)
  for (pock in pocket_pos_100){#length(names_ligand_unique)
    print(pock)
    lig_name = mat_sL[which(mat_sL$names_pockets == pock),2]
    index = as.character(mat_sL$names_pockets[which(mat_sL[,2] != lig_name)])
    if(length(index) > 2) {
      index = sample(index,2)
    }
    dist_pock = NULL
    dist_features_pock = NULL
    for (j in index) {
      if(pock != j) {
        if(sapply(strsplit(pock, "_"), "[", 1) != sapply(strsplit(j, "_"), "[", 1)) { #CHANGE TO != FOR DIFFERENT PROTEINS
          data_samePsameL = rbind(data_samePsameL, c(pock,j)) 
          #dist_pock = c(dist_pock, dist(rbind(dt[index[j],],
          #                                    dt[index[k],])))
          #dist_features_pock = rbind(dist_features_pock,  apply(dt[c(index[j],index[k]),],2,dist))
        }
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

data_diffPdiffL = dist_lig_diffP_diffL(dt,names_ligand, mat_sL, pocket_pos_100)
nrow(data_diffPdiffL)

#write.csv(data_diffPdiffL, "../data/data_structure_comparaison/data_samePsameL_100_neg.csv")
###
data_samePsameL_positive = read.csv("../data/data_structure_comparaison/data_samePsameL_100_pos.csv", colClasses = "character")
data_diffPdiffL_negative = read.csv("../data/data_structure_comparaison/data_samePsameL_100_neg.csv", colClasses = "character")
#
data_samePsameL_positive = data_samePsameL_positive[1:20,]
data_diffPdiffL_negative = data_diffPdiffL_negative[1:20,]
