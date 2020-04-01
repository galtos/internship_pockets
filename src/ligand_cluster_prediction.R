#### AUTOMATED PROTOCOL FOR LIGANDS POCKETS DESCRIPTION ####
#!/usr/bin/env Rscript

#
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  args[1] = "../data/pockets_MD_NS1/Res_pocketConf0101-p0_atm/pocketConf101-0_atm.des"
  #stop("At least one argument must be supplied : path to pocket des file (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "../results/automated_test"
  #default kmeans file
  args[3] = "../results/dt_12clean_tree_mcqueen_seeds800.Rdata"
  #default scale center and sd
  args[4] = "../results/scaled:center_dt12clean.Rdata"
  args[5] = "../results/scaled:scale_dt12clean.Rdata"
  args[6] = "../results/minmax_value_dt12clean.Rdata"
  #default number of clusters to show
  args[7] = 10  
}

#### IMPORT LIBRARIES ####
print("Import libraries ...")
library(data.tree)
#library(treemap)
library(fmsb)
radarchart(langues.means)
library(dendextend)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####LOADING TREE AND INFORMATIONS FOR THE SCALE ####
print("Loading the clustered data ...")
load(file = args[3])
print("Loading the scaled information ...")
scaled_center_dt_t = read.table(file = args[4], col.names = F, row.names = 1)
scaled_scale_dt_t = read.table(file = args[5], col.names = F, row.names = 1)

scaled_center_dt = scaled_center_dt_t[,1]
names(scaled_center_dt) = rownames(scaled_center_dt_t)
scaled_scale_dt = scaled_scale_dt_t[,1]
names(scaled_scale_dt) = rownames(scaled_scale_dt_t)

####LOADING INFORMATION NEW POCKET####
##DES file
dt_new_pocket_des = read.table(args[1], row.names = 1)
## DATA 12 DESCRIPTORS##
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
## SCALE NEW POCKET DATA ##
new_pocket = scale(new_pocket, scaled_center_dt[sort(names(scaled_center_dt))], scaled_scale_dt[sort(names(scaled_scale_dt))])
## COMPUTE DISTANCE 12 DESCRIPTORS ##
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
## GET THE nth CLOSEST CLUSTERS ##
Sort(alltree, "dist", decreasing = FALSE)
as.integer(alltree$Get("name", filterFun = isLeaf)[1:args[7]])
print("The distance of your pocket to the centers of the closest clusters :")
print(alltree$Get("dist", filterFun = isLeaf)[1:args[7]])
best_cluster = as.integer(alltree$Get("name", filterFun = isLeaf)[1])
## DATA VISUALIZATION ##
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
minmax_value = read.table(file = args[6])
colnames(minmax_value) = names_dt12
#get center values of first cluster for the minimum distance
min_dist = min(alltree$Get("dist", filterFun = isLeaf))
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
###!!!!!!!!!!!!! ADD SAVING IMAGE IN OUTPUT FILE
legend(x=0.6, y=1, legend = c("new_pocket", paste0("cluster:",best_cluster)), bty = "n", pch=20 , col=colors_in , cex=1, pt.cex=3)
###!!!!!!!!!
####HCLUST on closest cluster ####
hclust_best = hclust(dist(rbind(dt[unlist(alltree$children[[best_cluster]]$pockets_names),sort(colnames(dt))],
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
res_d_37 = apply(dt[unlist(alltree$children[[best_cluster]]$pockets_names),], 1, function(x) {dist(rbind(new_pocket,x))})
head(sort(res_d),10)
end_time <- Sys.time()

end_time - start_time

sapply(strsplit(names(head(sort(res_d),10)), "_"), "[", 2)
head(sort(res_d),10)


