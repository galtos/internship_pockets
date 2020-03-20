#### IMPORT LIBRARIES ####
library(data.tree)
library(treemap)
library(fmsb)
radarchart(langues.means)
####LOADING TREE AND INFORMATIONS FOR THE SCALE ####
load(file = "../results/dt_12clean_tree_mcqueen_seeds100.Rdata")
load(file = "../results/dt_12clean_tree_mcqueen_seeds800.Rdata")
#load(file = "../results/res_12desc/dt_12dsc_tree_withinss400_mcqueen.Rdata")
#load(file = "../results/dt_72dsc_tree_withinss400_mcqueen.Rdata")
#load(file = "../results/dt_72dsc_tree_withinss400_mcqueen_physicochemicals.Rdata")
#load(file = "../results/dt_72dsc_tree_withinss400_mcqueen_geometrical.Rdata")
#
scaled_center_dt_t = read.table(file = "../results/scaled:center_dt12clean.Rdata", col.names = F, row.names = 1)
scaled_scale_dt_t = read.table(file = "../results/scaled:scale_dt12clean.Rdata", col.names = F, row.names = 1)
scaled_center_dt = scaled_center_dt_t[,1]
names(scaled_center_dt) = rownames(scaled_center_dt_t)
scaled_scale_dt = scaled_scale_dt_t[,1]
names(scaled_scale_dt) = rownames(scaled_scale_dt_t)
####################################################
####LOADING INFORMATION NEW POCKET####
##DES file
dt_new_pocket_des = read.table("../data/pockets_MD_NS1/Res_pocketConf0101-p0_atm/pocketConf101-0_atm.des", row.names = 1)
dt_new_pocket_des = read.table("../data/pockets_MD_NS1/Res_pocketConf047-0_atm/pocketConf047-0_atm.des", row.names = 1)
dt_new_pocket_des = read.table("../data/pockets_MD_NS1/Res_pocketConf137-0_atm/pocket-Conf137-0_atm.des", row.names = 1)
dt_new_pocket_des = read.table("../data/pockets_MD_NS1/Res_pocketConf0093-2_atm/pocket-Conf0093-2_atm.des", row.names = 1)
##TXT after fpocket file
#covid19:
dt_new_pocket_txt = read.table("../data/pockets_COVID19/pocket_PPE/pocket_PPE.txt", header = T, sep = "\t", row.names = 1, fill=TRUE)
dt_new_pocket_txt = dt_new_pocket_txt["pocket2_atm",]

dt_new_pocket_desR = 
## DATA 12 file TXT ##
new_pocket = data.frame(centers.C_RESIDUES = dt_new_pocket_txt[1,"Nb.RES"],
                        centers.DIAMETER_HULL = dt_new_pocket_txt[1,"Diameter.hull"],
                        centers.hydrophobic_kyte = dt_new_pocket_txt[1,"Hydrophobic.kyte"],
                        centers.p_aliphatic_residues = dt_new_pocket_txt[1,"Aliphatic.residues"],
                        centers.p_aromatic_residues = dt_new_pocket_txt[1,"Aromatic.residues"],
                        centers.p_hydrophobic_residues = dt_new_pocket_txt[1,"Hydrophobic.residues"],
                        centers.p_Nlys_atom = dt_new_pocket_txt[1,"Nlys.atom"],
                        centers.p_Ntrp_atom = dt_new_pocket_txt[1,"Ntrp.atom"],
                        centers.p_Ooh_atom = dt_new_pocket_txt[1,"Ooh.atom"],
                        centers.p_Otyr_atom = dt_new_pocket_txt[1,"Otyr.atom"],
                        centers.p_polar_residues = dt_new_pocket_txt[1,"Polar.residues"],
                        centers.VOLUME_HULL = dt_new_pocket_txt[1,"Volume.hull"]
                        )
## DATA 12 ##
new_pocket = data.frame(centers.C_RESIDUES = dt_new_pocket_des["pocket_C_RESIDUES",1],
                        centers.DIAMETER_HULL = dt_new_pocket_des["pocket_DIAMETER_HULL",1],
                        centers.hydrophobic_kyte = dt_new_pocket_des["pocket_hydrophobic_kyte",1],
                        centers.p_aliphatic_residues = dt_new_pocket_des["pocket_p_aliphatic_residues",1],
                        centers.p_aromatic_residues = dt_new_pocket_des["pocket_p_aromatic_residues",1],
                        centers.p_hydrophobic_residues = dt_new_pocket_des["pocket_p_hydrophobic_residues",1],
                        centers.p_Nlys_atom = dt_new_pocket_des["pocket_p_Nlys_atom",1],
                        centers.p_Ntrp_atom = dt_new_pocket_des["pocket_p_Ntrp_atom",1],
                        centers.p_Ooh_atom = dt_new_pocket_des["pocket_p_Ooh_atom",1],
                        centers.p_Otyr_atom = dt_new_pocket_des["pocket_p_Otyr_atom",1],
                        centers.p_polar_residues = dt_new_pocket_des["pocket_p_polar_residues",1],
                        centers.VOLUME_HULL = dt_new_pocket_des["pocket_VOLUME_HULL",1]
                        )
## DATA 72 ##
new_pocket = data.frame(centers.A = dt_new_pocket_des["pocket_A",1],
                        centers.C = dt_new_pocket_des["pocket_C",1],
                        centers.C_ATOM = dt_new_pocket_des["pocket_C_ATOM",1],
                        centers.C_RESIDUES = dt_new_pocket_des["pocket_C_RESIDUES",1],
                        centers.charge = dt_new_pocket_des["pocket_charge",1],
                        centers.CONVEX.SHAPE_COEFFICIENT = dt_new_pocket_des["pocket_CONVEX-SHAPE_COEFFICIENT",1],
                        centers.D = dt_new_pocket_des["pocket_D",1],
                        centers.DIAMETER_HULL = dt_new_pocket_des["pocket_DIAMETER_HULL",1],
                        centers.E = dt_new_pocket_des["pocket_E",1],
                        centers.F = dt_new_pocket_des["pocket_F",1],
                        centers.FACE = dt_new_pocket_des["pocket_FACE",1],
                        centers.G = dt_new_pocket_des["pocket_G",1],
                        centers.H = dt_new_pocket_des["pocket_H",1],
                        centers.hydrophobic_kyte = dt_new_pocket_des["pocket_hydrophobic_kyte",1],
                        centers.hydrophobicity = dt_new_pocket_des["pocket_hydrophobicity",1],
                        centers.I = dt_new_pocket_des["pocket_I",1],
                        centers.INERTIA_1 = dt_new_pocket_des["pocket_INERTIA_1",1],
                        centers.INERTIA_2 = dt_new_pocket_des["pocket_INERTIA_2",1],
                        centers.INERTIA_3 = dt_new_pocket_des["pocket_INERTIA_3",1],
                        centers.K = dt_new_pocket_des["pocket_K",1],
                        centers.L = dt_new_pocket_des["pocket_L",1],
                        centers.M = dt_new_pocket_des["pocket_M",1],
                        centers.N = dt_new_pocket_des["pocket_N",1],
                        centers.P = dt_new_pocket_des["pocket_P",1],
                        centers.p_aliphatic_residues = dt_new_pocket_des["pocket_p_aliphatic_residues",1],
                        centers.p_aromatic_residues = dt_new_pocket_des["pocket_p_aromatic_residues",1],
                        centers.p_C_atom = dt_new_pocket_des["pocket_p_C_atom",1],
                        centers.p_Car_atom = dt_new_pocket_des["pocket_p_Car_atom",1],
                        centers.p_carbone_atom = dt_new_pocket_des["pocket_p_carbone_atom",1],
                        centers.p_Carg_atom = dt_new_pocket_des["pocket_p_Carg_atom",1],
                        centers.p_Ccoo_atom = dt_new_pocket_des["pocket_p_Ccoo_atom",1],
                        centers.p_Cgln_atom = dt_new_pocket_des["pocket_p_Cgln_atom",1],
                        centers.p_charged_residues = dt_new_pocket_des["pocket_p_charged_residues",1],
                        centers.p_hyd_atom = dt_new_pocket_des["pocket_p_hyd_atom",1],
                        centers.p_hydrophobic_atom = dt_new_pocket_des["pocket_p_hydrophobic_atom",1],
                        centers.p_hydrophobic_residues = dt_new_pocket_des["pocket_p_hydrophobic_residues",1],
                        centers.p_main_chain_atom = dt_new_pocket_des["pocket_p_main_chain_atom",1],
                        centers.p_N_atom = dt_new_pocket_des["pocket_p_N_atom",1],
                        centers.p_ND1_atom = dt_new_pocket_des["pocket_p_ND1_atom",1],
                        centers.p_NE2_atom = dt_new_pocket_des["pocket_p_NE2_atom",1],
                        centers.p_negative_residues = dt_new_pocket_des["pocket_p_negative_residues",1],
                        centers.p_nitrogen_atom = dt_new_pocket_des["pocket_p_nitrogen_atom",1],
                        centers.p_Nlys_atom = dt_new_pocket_des["pocket_p_Nlys_atom",1],
                        centers.p_Ntrp_atom = dt_new_pocket_des["pocket_p_Ntrp_atom",1],
                        centers.p_O_atom = dt_new_pocket_des["pocket_p_O_atom",1],
                        centers.p_Ocoo_atom = dt_new_pocket_des["pocket_p_Ocoo_atom",1],
                        centers.p_Ooh_atom = dt_new_pocket_des["pocket_p_Ooh_atom",1],
                        centers.p_Otyr_atom = dt_new_pocket_des["pocket_p_Otyr_atom",1],
                        centers.p_oxygen_atom = dt_new_pocket_des["pocket_p_oxygen_atom",1],
                        centers.p_polar_residues = dt_new_pocket_des["pocket_p_polar_residues",1],
                        centers.p_positive_residues = dt_new_pocket_des["pocket_p_positive_residues",1],
                        centers.p_S_atom = dt_new_pocket_des["pocket_p_S_atom",1],
                        centers.p_side_chain_atom = dt_new_pocket_des["pocket_p_side_chain_atom",1],
                        centers.p_small_residues = dt_new_pocket_des["pocket_p_small_residues",1],
                        centers.p_sulfur_atom = dt_new_pocket_des["pocket_p_sulfur_atom",1],
                        centers.p_tiny_residues = dt_new_pocket_des["pocket_p_tiny_residues",1],
                        centers.PCI = dt_new_pocket_des["pocket_PCI",1],
                        centers.polarity = dt_new_pocket_des["pocket_polarity",1],
                        centers.PSI = dt_new_pocket_des["pocket_PSI",1],
                        centers.Q = dt_new_pocket_des["pocket_Q",1],
                        centers.R = dt_new_pocket_des["pocket_R",1],
                        centers.RADIUS_CYLINDER = dt_new_pocket_des["pocket_RADIUS_CYLINDER",1],
                        centers.RADIUS_HULL = dt_new_pocket_des["pocket_RADIUS_HULL",1],
                        centers.S = dt_new_pocket_des["pocket_S",1],
                        centers.SMALLEST_SIZE = dt_new_pocket_des["pocket_SMALLEST_SIZE",1],
                        centers.SURFACE_HULL = dt_new_pocket_des["pocket_SURFACE_HULL",1],
                        centers.T = dt_new_pocket_des["pocket_T",1],
                        centers.V = dt_new_pocket_des["pocket_V",1],
                        centers.W = dt_new_pocket_des["pocket_W",1],
                        centers.X._ATOM_CONVEXE = dt_new_pocket_des["pocket_%_ATOM_CONVEXE",1],
                        centers.Y = dt_new_pocket_des["pocket_Y",1]
)

## DT PHYSICOCHYMIQUES ##
new_pocket = data.frame(centers.p_aliphatic_residues = dt_new_pocket_des["pocket_p_aliphatic_residues",1],
                        centers.p_ND1_atom = dt_new_pocket_des["pocket_p_ND1_atom",1],
                        centers.p_NE2_atom = dt_new_pocket_des["pocket_p_NE2_atom",1],
                        centers.p_Nlys_atom = dt_new_pocket_des["pocket_p_Nlys_atom",1],
                        centers.p_Ntrp_atom = dt_new_pocket_des["pocket_p_Ntrp_atom",1],
                        centers.p_Ooh_atom = dt_new_pocket_des["pocket_p_Ooh_atom",1],
                        centers.p_Otyr_atom = dt_new_pocket_des["pocket_p_Otyr_atom",1]
)
## DT GEOMETRIQUES ##
new_pocket = data.frame(centers.C_ATOM = dt_new_pocket_des["pocket_C_ATOM",1],
                          centers.C_RESIDUES = dt_new_pocket_des["pocket_C_RESIDUES",1],
                          centers.CONVEX.SHAPE_COEFFICIENT = dt_new_pocket_des["pocket_CONVEX-SHAPE_COEFFICIENT",1],
                          centers.RADIUS_CYLINDER = dt_new_pocket_des["pocket_RADIUS_CYLINDER",1],
                          centers.RADIUS_HULL = dt_new_pocket_des["pocket_RADIUS_HULL",1],
                          centers.SURFACE_HULL = dt_new_pocket_des["pocket_SURFACE_HULL",1]
)


#new_pocket_1 = scale(new_pocket_1, scaled_center_dt[sort(names(scaled_center_dt))], scaled_scale_dt[sort(names(scaled_scale_dt))])
#new_pocket_2 = scale(new_pocket_2, scaled_center_dt[sort(names(scaled_center_dt))], scaled_scale_dt[sort(names(scaled_scale_dt))])
#dist(rbind(new_pocket_1,new_pocket_2))

new_pocket = scale(new_pocket, scaled_center_dt[sort(names(scaled_center_dt))], scaled_scale_dt[sort(names(scaled_scale_dt))])
##when scaled names need to be sorted (dt 72)
#new_pocket = scale(new_pocket, scaled_center_dt[sort(names(scaled_center_dt))], scaled_scale_dt[sort(names(scaled_scale_dt))])

## DISTANCE 12 ##
alltree$Do(function(node) node$dist <- dist(rbind(c(
  node$centers.p_polar_residues,
  node$centers.p_Nlys_atom,
  node$centers.p_aliphatic_residues,
  node$centers.VOLUME_HULL,
  node$centers.DIAMETER_HULL,
  node$centers.p_Ooh_atom,
  node$centers.hydrophobic_kyte,
  node$centers.p_Ntrp_atom,
  node$centers.p_Otyr_atom,
  node$centers.C_RESIDUES,
  node$centers.p_aromatic_residues,
  node$centers.p_hydrophobic_residues
  ),
  new_pocket)))
## DISTANCE 72 ##
alltree$Do(function(node) node$dist <- dist(rbind(c(
  node$centers.A,
  node$centers.C,
  node$centers.C_ATOM,
  node$centers.C_RESIDUES,
  node$centers.charge,
  node$centers.CONVEX.SHAPE_COEFFICIENT,
  node$centers.D,
  node$centers.DIAMETER_HULL,
  node$centers.E,
  node$centers.F,
  node$centers.FACE,
  node$centers.G,
  node$centers.H,
  node$centers.hydrophobic_kyte,
  node$centers.hydrophobicity,
  node$centers.I,
  node$centers.INERTIA_1,
  node$centers.INERTIA_2,
  node$centers.INERTIA_3,
  node$centers.K,
  node$centers.L,
  node$centers.M,
  node$centers.N,
  node$centers.P,
  node$centers.p_aliphatic_residues,
  node$centers.p_aromatic_residues,
  node$centers.p_C_atom,
  node$centers.p_Car_atom,
  node$centers.p_carbone_atom,
  node$centers.p_Carg_atom,
  node$centers.p_Ccoo_atom,
  node$centers.p_Cgln_atom,
  node$centers.p_charged_residues,
  node$centers.p_hyd_atom,
  node$centers.p_hydrophobic_atom,
  node$centers.p_hydrophobic_residues,
  node$centers.p_main_chain_atom,
  node$centers.p_N_atom,
  node$centers.p_ND1_atom,
  node$centers.p_NE2_atom,
  node$centers.p_negative_residues,
  node$centers.p_nitrogen_atom,
  node$centers.p_Nlys_atom,
  node$centers.p_Ntrp_atom,
  node$centers.p_O_atom,
  node$centers.p_Ocoo_atom,
  node$centers.p_Ooh_atom,
  node$centers.p_Otyr_atom,
  node$centers.p_oxygen_atom,
  node$centers.p_polar_residues,
  node$centers.p_positive_residues,
  node$centers.p_S_atom,
  node$centers.p_side_chain_atom,
  node$centers.p_small_residues,
  node$centers.p_sulfur_atom,
  node$centers.p_tiny_residues,
  node$centers.PCI,
  node$centers.polarity,
  node$centers.PSI,
  node$centers.Q,
  node$centers.R,
  node$centers.RADIUS_CYLINDER,
  node$centers.RADIUS_HULL,
  node$centers.S,
  node$centers.SMALLEST_SIZE,
  node$centers.SURFACE_HULL,
  node$centers.T,
  node$centers.V,
  node$centers.W,
  node$centers.X._ATOM_CONVEXE,
  node$centers.Y
  ),
  new_pocket)))

## DISTANCE DT PHYSICOCHYMIQUES ##
alltree$Do(function(node) node$dist <- dist(rbind(c(
  node$centers.p_aliphatic_residues,
  node$centers.p_ND1_atom,
  node$centers.p_NE2_atom,
  node$centers.p_Nlys_atom,
  node$centers.p_Ntrp_atom,
  node$centers.p_Ooh_atom,
  node$centers.p_Otyr_atom
  ),
  new_pocket)))
## DISTANCE DT GEOMETRIQUES ##
alltree$Do(function(node) node$dist <- dist(rbind(c(
  node$centers.C_ATOM,
  node$centers.C_RESIDUES,
  node$centers.CONVEX.SHAPE_COEFFICIENT,
  node$centers.RADIUS_CYLINDER,
  node$centers.RADIUS_HULL,
  node$centers.SURFACE_HULL
  ),
  new_pocket)))
## PHARMACOPHORES SIMILARITY ##
v_ph_1 = unlist(alltree$`1`$`1`$pharmacophores_consensus_mean)
v_ph_2 = unlist(alltree$`2`$pharmacophores_consensus_mean)
length(v_ph_1)
alltree$Do(function(node) {
  v_ph_2 = unlist(node$pharmacophores_consensus_mean)
  if(length(v_ph_1) == length(v_ph_2)){
    t = table(v_ph_1, v_ph_2)
    print(node$path)
    #print(v_ph_2)
    print(t)
    if(ncol(t) != 1) {
      common_non_null = sum(t[2:nrow(t),2:ncol(t)])
      non_null_counts_Fa = sum(t[2:nrow(t),])
      non_null_counts_Fb = sum(t[,2:ncol(t)])
      node$val_similarity = common_non_null/min(non_null_counts_Fa,non_null_counts_Fb)
    }else {
      node$val_similarity = 0
    }
  } else {
    node$val_similarity = 0
  }
})
alltree$`9`$`10`$`9`$pharmacophores_consensus_mean
alltree$`9`$`10`$`9`$withinss
w_tree = alltree$Get("withinss", filterFun = isLeaf)
s_tree = alltree$Get("size", filterFun = isLeaf)
length(s_tree)
plot(s_tree, w_tree)

Sort(alltree, "dist", decreasing = FALSE)
print(alltree, "size", "withinss", "dist")

alltree$Get("dist", filterFun = function(node) node$level == 2)
tot = alltree$Get("totss", filterFun = isLeaf)
length(tot)
#### DATA VISUALIZATION ####
names_dt12 = c( "centers.p_polar_residues",
                "centers.p_Nlys_atom",
                "centers.p_aliphatic_residues",
                "centers.VOLUME_HULL",
                "centers.DIAMETER_HULL",
                "centers.p_Ooh_atom",
                "centers.hydrophobic_kyte",
                "centers.p_Ntrp_atom",
                "centers.p_Otyr_atom",
                "centers.C_RESIDUES",
                "centers.p_aromatic_residues",
                "centers.p_hydrophobic_residues")
#spiders of the clusters#
minmax_value = read.table(file = "../results/minmax_value_dt12.Rdata")
colnames(minmax_value) = names_dt12
colnames(new_pocket) = names_dt12
#get center values of first cluster for the minimum distance
min_dist = min(alltree$Get("dist", filterFun = function(node) node$level == 2))
cluster_center_min_dist = alltree$Get(function(node) c(
                                            centers.p_polar_residues = node$centers.p_polar_residues,
                                            centers.p_Nlys_atom = node$centers.p_Nlys_atom,
                                            centers.p_aliphatic_residues = node$centers.p_aliphatic_residues,
                                            centers.VOLUME_HULL =  node$centers.VOLUME_HULL,
                                            centers.DIAMETER_HULL =  node$centers.DIAMETER_HULL,
                                            centers.p_Ooh_atom = node$centers.p_Ooh_atom,
                                            centers.hydrophobic_kyte = node$centers.hydrophobic_kyte,
                                            centers.p_Ntrp_atom = node$centers.p_Ntrp_atom,
                                            centers.p_Otyr_atom = node$centers.p_Otyr_atom,
                                            centers.C_RESIDUES = node$centers.C_RESIDUES,
                                            centers.p_aromatic_residues = node$centers.p_aromatic_residues,
                                            centers.p_hydrophobic_residues = node$centers.p_hydrophobic_residues),
                                    filterFun = function(node) {
                                      if(node$level == 2 && node$dist == min_dist) {
                                        return(TRUE)
                                      } else {
                                        return(FALSE)
                                      }
                                      })
dt.visualize = NULL
dt.visualize = rbind(minmax_value, new_pocket)
dt.visualize = rbind(dt.visualize, t(cluster_center_min_dist))
#
#visualize firts 10 clusters
cluster_center_10_first = alltree$Get(function(node) c(
  centers.p_polar_residues = node$centers.p_polar_residues,
  centers.p_Nlys_atom = node$centers.p_Nlys_atom,
  centers.p_aliphatic_residues = node$centers.p_aliphatic_residues,
  centers.VOLUME_HULL =  node$centers.VOLUME_HULL,
  centers.DIAMETER_HULL =  node$centers.DIAMETER_HULL,
  centers.p_Ooh_atom = node$centers.p_Ooh_atom,
  centers.hydrophobic_kyte = node$centers.hydrophobic_kyte,
  centers.p_Ntrp_atom = node$centers.p_Ntrp_atom,
  centers.p_Otyr_atom = node$centers.p_Otyr_atom,
  centers.C_RESIDUES = node$centers.C_RESIDUES,
  centers.p_aromatic_residues = node$centers.p_aromatic_residues,
  centers.p_hydrophobic_residues = node$centers.p_hydrophobic_residues),
  filterFun = function(node) node$level == 2)
  
dt.visualize = NULL
dt.visualize =  rbind(minmax_value, new_pocket)
dt.visualize = rbind(dt.visualize, t(cluster_center_10_first))
# Color vector
colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )
radarchart(dt.visualize, axistype=1 , 
           #custom polygon
           pcol=colors_border ,plwd=4 , plty=1,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(-1,1,0.5), cglwd=0.8,
           #custom labels
           vlcex=0.8 
           )
# Add a legend
legend(x=0.6, y=1, legend = c("pocket_02", "cluster_4"), bty = "n", pch=20 , col=colors_in , cex=1.2, pt.cex=3)



#### VALIDATE TREE WITH TOUGH-M1 DATASET ####

dt_TOUGH_positive = read.table(file = "../data/dataset_TOUGH/TOUGH-M1_positive.list", sep=" ")
dt_TOUGH_negative = read.table(file = "../data/dataset_TOUGH/TOUGH-M1_negative.list", sep=" ")


dt_TOUGH_positive[,1] = substr(dt_TOUGH_positive[,1],1,4)
dt_TOUGH_positive[,2] = substr(dt_TOUGH_positive[,2],1,4)
head(dt_TOUGH_positive)

dt_TOUGH_negative[,1] = substr(dt_TOUGH_negative[,1],1,4)
dt_TOUGH_negative[,2] = substr(dt_TOUGH_negative[,2],1,4)
head(dt_TOUGH_negative)
dt_TOUGH_positive = dt_TOUGH_negative

vector_positive_c1 = NULL
for (i in 1:nrow(dt_TOUGH_positive)) {
  name_1 = 0
  name_2 = 0   
  for (j in 1:alltree$count) {
    if(length(grep(toupper(dt_TOUGH_positive[i,1]), unlist(alltree$children[[j]]$pockets_names))) != 0) {
      #print(toupper(dt_TOUGH_positive[i,1]))
      name_1 = j
    }
    if(length(grep(toupper(dt_TOUGH_positive[i,2]), unlist(alltree$children[[j]]$pockets_names))) != 0) {
      #print(toupper(dt_TOUGH_positive[i,2]))
      name_2 = j
    }
   
  }
  if(name_1 != 0 && name_2 != 0 && name_1 == name_2) {
    vector_positive_c1 = c(vector_positive_c1,1)
    print(i)
    print("ITS A MATCH")
  }
  if(name_1 != 0 && name_2 != 0 && name_1 != name_2) {
    vector_positive_c1 = c(vector_positive_c1,0)
    print(i)
    print("ITS NOT A MATCH")
  } 
}
length(vector_positive_c1)
length(which(vector_positive_c1 == 1))

#TOUGH negative : 17276 positifs / 153834 au total

print(grep(toupper(dt_TOUGH_positive[i,2]), unlist(alltree$children[[4]]$pockets_names)))
print(grep("2QEN", unlist(alltree$children[[4]]$pockets_names)))
name_test_4 =  unlist(alltree$children[[4]]$pockets_names)
name_test_4[3084]

print(grep(toupper(dt_TOUGH_positive[i,2]), unlist(alltree$`4`$pockets_names)))

t_0 = grep(toupper(dt_TOUGH_positive[i,1]), unlist(alltree$`4`$pockets_names), value = TRUE)
any(unlist(alltree$`4`$pockets_names) == toupper(dt_TOUGH_positive[i,2]))
#### VALIDATE TREE WITH TOUGH-C1 DATASET ####
names_TOUGH_C1_nucleotide = read.table(file = "../data/dataset_TOUGH-C1/osfstorage-archive/nucleotide.list", sep=" ")
names_TOUGH_C1_heme = read.table(file = "../data/dataset_TOUGH-C1/osfstorage-archive/heme.list", sep=" ")

n_nucleotide_cluster = NULL
for (i in 1:alltree$averageBranchingFactor) {
  n_nucleotide_cluster = c(n_nucleotide_cluster, length(intersect(toupper(substr(unlist(alltree$children[[i]]$pockets_names), 1, 4)),
                                                                  toupper(substr(levels(names_TOUGH_C1_nucleotide[1,1]),1,4)))))
}
print(n_nucleotide_cluster)

n_heme_cluster = NULL
for (i in 1:alltree$averageBranchingFactor) {
  n_heme_cluster = c(n_heme_cluster, length(intersect(toupper(substr(unlist(alltree$children[[i]]$pockets_names), 1, 9)),
                                                                  paste0(toupper(substr(levels(names_TOUGH_C1_heme[1,1]),1,4)), "_HEM_" ))))
  print(alltree$children[[i]]$path)
}
print(n_heme_cluster)


####HCLUST on Centroids ####
alltree$`4`$levelName
print(alltree$`4`,alltree$`4`$fields)

## DT 12 ##
dt.kmean_centers_10 = alltree$Get(function(node) c(
  centers.p_polar_residues = node$centers.p_polar_residues,
  centers.p_Nlys_atom = node$centers.p_Nlys_atom,
  centers.p_aliphatic_residues = node$centers.p_aliphatic_residues,
  centers.VOLUME_HULL =  node$centers.VOLUME_HULL,
  centers.DIAMETER_HULL =  node$centers.DIAMETER_HULL,
  centers.p_Ooh_atom = node$centers.p_Ooh_atom,
  centers.hydrophobic_kyte = node$centers.hydrophobic_kyte,
  centers.p_Ntrp_atom = node$centers.p_Ntrp_atom,
  centers.p_Otyr_atom = node$centers.p_Otyr_atom,
  centers.C_RESIDUES = node$centers.C_RESIDUES,
  centers.p_aromatic_residues = node$centers.p_aromatic_residues,
  centers.p_hydrophobic_residues = node$centers.p_hydrophobic_residues),
  filterFun = function(node) node$level == 2)
## DT 72 ##
dt.kmean_centers_10 = alltree$Get(function(node) c(
  centers.A = node$centers.A,
  centers.C = node$centers.C,
  centers.C_ATOM = node$centers.C_ATOM,
  centers.C_RESIDUES = node$centers.C_RESIDUES,
  centers.charge=node$centers.charge,
  centers.CONVEX.SHAPE_COEFFICIENT=node$centers.CONVEX.SHAPE_COEFFICIENT,
  centers.D=node$centers.D,
  centers.DIAMETER_HULL=node$centers.DIAMETER_HULL,
  centers.E=node$centers.E,
  centers.F=node$centers.F,
  centers.FACE=node$centers.FACE,
  centers.G=node$centers.G,
  centers.H=node$centers.H,
  centers.hydrophobic_kyte=node$centers.hydrophobic_kyte,
  centers.hydrophobicity=node$centers.hydrophobicity,
  centers.I=node$centers.I,
  centers.INERTIA_1=node$centers.INERTIA_1,
  centers.INERTIA_2=node$centers.INERTIA_2,
  centers.INERTIA_3=node$centers.INERTIA_3,
  centers.K=node$centers.K,
  centers.L=node$centers.L,
  centers.M=node$centers.M,
  centers.N=node$centers.N,
  centers.P=node$centers.P,
  centers.p_aliphatic_residues=node$centers.p_aliphatic_residues,
  centers.p_aromatic_residues=node$centers.p_aromatic_residues,
  centers.p_C_atom=node$centers.p_C_atom,
  centers.p_Car_atom=node$centers.p_Car_atom,
  centers.p_carbone_atom=node$centers.p_carbone_atom,
  centers.p_Carg_atom=node$centers.p_Carg_atom,
  centers.p_Ccoo_atom=node$centers.p_Ccoo_atom,
  centers.p_Cgln_atom=node$centers.p_Cgln_atom,
  centers.p_charged_residues=node$centers.p_charged_residues,
  centers.p_hyd_atom=node$centers.p_hyd_atom,
  centers.p_hydrophobic_atom=node$centers.p_hydrophobic_atom,
  centers.p_hydrophobic_residues=node$centers.p_hydrophobic_residues,
  centers.p_main_chain_atom=node$centers.p_main_chain_atom,
  centers.p_N_atom=node$centers.p_N_atom,
  centers.p_ND1_atom=node$centers.p_ND1_atom,
  centers.p_NE2_atom=node$centers.p_NE2_atom,
  centers.p_negative_residues=node$centers.p_negative_residues,
  centers.p_nitrogen_atom=node$centers.p_nitrogen_atom,
  centers.p_Nlys_atom=node$centers.p_Nlys_atom,
  centers.p_Ntrp_atom=node$centers.p_Ntrp_atom,
  centers.p_O_atom=node$centers.p_O_atom,
  centers.p_Ocoo_atom=node$centers.p_Ocoo_atom,
  centers.p_Ooh_atom=node$centers.p_Ooh_atom,
  centers.p_Otyr_atom=node$centers.p_Otyr_atom,
  centers.p_oxygen_atom=node$centers.p_oxygen_atom,
  centers.p_polar_residues=node$centers.p_polar_residues,
  centers.p_positive_residues=node$centers.p_positive_residues,
  centers.p_S_atom=node$centers.p_S_atom,
  centers.p_side_chain_atom=node$centers.p_side_chain_atom,
  centers.p_small_residues=node$centers.p_small_residues,
  centers.p_sulfur_atom=node$centers.p_sulfur_atom,
  centers.p_tiny_residues=node$centers.p_tiny_residues,
  centers.PCI=node$centers.PCI,
  centers.polarity=node$centers.polarity,
  centers.PSI=node$centers.PSI,
  centers.Q=node$centers.Q,
  centers.R=node$centers.R,
  centers.RADIUS_CYLINDER=node$centers.RADIUS_CYLINDER,
  centers.RADIUS_HULL=node$centers.RADIUS_HULL,
  centers.S=node$centers.S,
  centers.SMALLEST_SIZE=node$centers.SMALLEST_SIZE,
  centers.SURFACE_HULL=node$centers.SURFACE_HULL,
  centers.T=node$centers.T,
  centers.V=node$centers.V,
  centers.W=node$centers.W,
  centers.X._ATOM_CONVEXE=node$centers.X._ATOM_CONVEXE,
  centers.Y=node$centers.Y),
  filterFun = function(node) node$level == 2)

## DT PHYSICO ##
dt.kmean_centers_10 = alltree$Get(function(node) c(
  centers.p_aliphatic_residues = node$centers.p_aliphatic_residues,
  centers.p_ND1_atom = node$centers.p_ND1_atom,
  centers.p_NE2_atom = node$centers.p_NE2_atom,
  centers.p_Nlys_atom =  node$centers.p_Nlys_atom,
  centers.p_Ntrp_atom =  node$centers.p_Ntrp_atom,
  centers.p_Ooh_atom = node$centers.p_Ooh_atom,
  centers.p_Otyr_atom = node$centers.p_Otyr_atom),
  filterFun = function(node) node$level == 2) #&& strtoi(node$parent$name) == 7)
##
dt.kmean_centers_10_names = alltree$Get(function(node) node$levelName, filterFun = function(node) node$level == 2)# && strtoi(node$parent$name) == 7)

dt_centers.hclust = hclust(dist(t(dt.kmean_centers_10)), method = "ward.D2")
plot(dt_centers.hclust)

####HCLUST on closest cluster ####
library(dendextend)

hclust_best = hclust(dist(rbind(dt[unlist(alltree$`193`$pockets_names),sort(colnames(dt))],
                          new_pocket)),
                     method = "ward.D2")

which(hclust_best$labels == "")
hclust_best$labels[which(hclust_best$labels == "")] = "new_pocket"
tail(hclust_best$labels)

hclust_best = as.dendrogram(hclust_best)

labels_colors(hclust_best) = 1
labels_colors(hclust_best)["new_pocket"] = 2

hclust_best <- set(hclust_best, "labels_cex", 0.8)
plot(hclust_best)

#apply dist
start_time <- Sys.time()
res_d = apply(dt, 1, function(x) {dist(rbind(new_pocket,x))})
head(sort(res_d))
end_time <- Sys.time()

end_time - start_time

start_time <- Sys.time()
res_d_37 = apply(dt[unlist(alltree$`193`$pockets_names),], 1, function(x) {dist(rbind(new_pocket,x))})
head(sort(res_d),10)
end_time <- Sys.time()

end_time - start_time

sapply(strsplit(names(head(sort(res_d),10)), "_"), "[", 2)
head(sort(res_d),10)

