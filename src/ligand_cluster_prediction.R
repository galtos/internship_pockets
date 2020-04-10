#### AUTOMATED PROTOCOL FOR LIGANDS POCKETS DESCRIPTION ####
#"C:\Program Files\R\R-3.6.3\bin\Rscript.exe"
## Collect arguments
args <- commandArgs(TRUE)
## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}
## Help section
if("--help" %in% args) {
  cat("
      The R Script
      
      Arguments:
      [TO MODIFY]
      --des              - path to des file of the query pocket
      --o                - output path and name file (you should take off the extension)
      --name             - name of the pocket
      --n_show           - number of closest clusters to show in stdout with their distance
      [DO NOT MODIFY]
      --kmeantree        - [data.tree] the tree containing the informations
      --dt               - [csv] print this text
      --scale_sd         - [Rdata] data to scale the new pocket data
      --scale_center     - [Rdata] data to scale the new pocket data
      --scale_minmax     - [Rdata] data to scale the new pocket data
      
      Example:
      Rscript --des=\"../data/pocket.des\" --o=\"../results/prediction_\" --name=\"covid19_pocket\" \n\n")
  q(save="no")
}
 
## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
args = argsL

if(is.null(argsL$des)) {
  args$des = "../data/pockets_MD_NS1/Res_pocketConf0101-p0_atm/pocketConf101-0_atm.des"
}
print(argsL$des)
## Arg2 default
if(is.null(argsL$o)) {
  args$o = "../results/prediction_validation/hclust_centroids_"
}
## Arg3 default
if(is.null(argsL$name)) {
  args$name = "pocketConf0101"
}
## Arg4 default
if(is.null(argsL$n_show)) {
  args$n_show = 10
}
## Arg5 default
if(is.null(argsL$kmeantree)) {
  args$kmeantree =  "../results/K_Means_comparaison/dt_12clean_tree_classic_seeds800_dend.Rdata"
}
## Arg6 default
if(is.null(argsL$dt)) {
  args$dt = "../data/dt_12clean.csv"
}
## Arg7 default
if(is.null(argsL$scale_sd)) {
  args$scale_sd = "../results/scaled:scale_dt12clean.Rdata"
}
## Arg8 default
if(is.null(argsL$scale_center)) {
  args$scale_center = "../results/scaled:center_dt12clean.Rdata"
}
## Arg8 default
if(is.null(argsL$scale_minmax)) {
  args$scale_minmax = "../results/minmax_value_dt12clean.Rdata"
}
#### IMPORT LIBRARIES ####
print("Import libraries ...")
.libPaths(c( .libPaths(), "C:/Users/Guillaume/Documents/R/win-library/3.6"))
print(.libPaths())
library(data.tree)
#library(nnet)
#library(fmsb)
#radarchart(langues.means)
suppressPackageStartupMessages(library(dendextend))
library(colorspace)
library(plotfunctions)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####LOADING TREE AND INFORMATIONS FOR THE SCALE ####
print("Loading the clustered data ...")
load(file = args$kmeantree)
print("Loading the pocket data of the 12 descriptors ...")
dt = read.csv('../data/dt_12clean.csv',row.names = 1)
print("Loading the scaled information ...")
scaled_center_dt_t = read.table(file = args$scale_center, col.names = F, row.names = 1)
scaled_scale_dt_t = read.table(file = args$scale_sd, col.names = F, row.names = 1)

scaled_center_dt = scaled_center_dt_t[,1]
names(scaled_center_dt) = rownames(scaled_center_dt_t)
scaled_scale_dt = scaled_scale_dt_t[,1]
names(scaled_scale_dt) = rownames(scaled_scale_dt_t)

####LOADING INFORMATION NEW POCKET####
##DES file
dt_new_pocket_des = read.table(args$des, row.names = 1)
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
  node$centers.C_RESIDUES,
  node$centers.DIAMETER_HULL,
  node$centers.hydrophobic_kyte,
  node$centers.p_aliphatic_residues,
  node$centers.p_aromatic_residues,
  node$centers.p_hydrophobic_residues,
  node$centers.p_Nlys_atom,
  node$centers.p_Ntrp_atom,
  node$centers.p_Ooh_atom,
  node$centers.p_Otyr_atom,
  node$centers.p_polar_residues,
  node$centers.VOLUME_HULL
  ),
  new_pocket)))
## GET THE nth CLOSEST CLUSTERS ##
#as.integer(alltree$Get("name", filterFun = isLeaf)[1:args[7]])
print("The distance of your pocket to the centers of the closest clusters :")
print(sort(alltree$Get("dist", filterFun = isLeaf))[1:args$n_show])
best_cluster = as.integer(alltree$Get("name", filterFun = isLeaf)[1])
####PLOT DIST ON HCLUST CENTROIDS ####
print("Saving pdf file of the distance of your pocket to the centroids...")
#dendrogram
dt_centers.dend = alltree$cluster_dend[[1]]#as.dendrogram(dt_centers.hclust)
#change label name
pock_dist = alltree$Get("dist", filterFun = isLeaf)
#
centroid_green = which(labels(dt_centers.dend) == as.integer(names(which(pock_dist == min(pock_dist)))))
lab_cluster_dend = as.integer(labels(dt_centers.dend))
labels(dt_centers.dend) <- paste0(lab_cluster_dend, 
                                  paste0(";Dist:",
                                  round(pock_dist[as.character(lab_cluster_dend)], 2)))
labels(dt_centers.dend)[centroid_green] <- paste0("POCHE->;GRP",
                                                  paste0(lab_cluster_dend[centroid_green], 
                                                  paste0(";Dist:",
                                                         round(min(pock_dist), 2))))
#change label color
color.gradient <- function(x, colors=c("darkgreen","gold","white"), colsteps=100) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),mean(x), length.out=colsteps)) ] )
}
c_colors = color.gradient(pock_dist)
color.df<-data.frame(COLOR_VALUE=pock_dist, color.name=color.gradient(pock_dist))
labels_colors(dt_centers.dend) <- as.character(color.df[as.character(order.dendrogram(dt_centers.dend)), "color.name"])
# Open a PDF for plotting; units are inches by default
pdf(paste0(paste0(args$o,args$name),".pdf"), width=40, height=15)
# Do some plotting
par(cex=0.3, mar=c(5, 9, 7, 1))
plot(dt_centers.dend, type = "rectangle", ylab = "Height", cex.lab = 0.7, cex.main = 8,
     main = paste0(paste0("Distance de ligands ",args$name)," dans les groupes"))
legend("topright", title = "Distance de la poche aux centroides",
       legend = c(paste0("max:", round(max(pock_dist),2)),
                  rep("",3),
                  paste0("mean:", round(mean(pock_dist),2)),
                  rep("",4),
                  paste0("min:", round(min(pock_dist),2))), pt.cex = 15, cex = 8, bty = "n")
legend("topleft", title = "10 groupes les plus proches:",
       legend = paste("Groupe", 
                      paste(names(sort(pock_dist)[1:10]),
                            paste("| Dist", round(sort(pock_dist)[1:10],2)))),
       pt.cex = 15, cex = 8, bty = "n")
gradientLegend(valRange=c(-14,14), color=sort(c_colors), pos=c(680,40,730,63), coords=TRUE, 
               border.col=alpha('gray'), side=4)#pos.num = 3, length=.2, depth=.02)
garbage = dev.off()
####PLOT DIST ON HCLUST CLUSTERS ####
print("Saving pdf file of the distance of your pocket to the chosen groups...")
print("These groups should be chosen after watching the previous plot...")
cat("Enter the numbers for the groups you want to plot separated by comma (ex:1,9,480) : \n")
f <- file("stdin")
open(f)
groups_cluster <- readLines(f,1)
close(f)
groups_cluster <- strsplit(groups_cluster,',')
groups_cluster <- as.numeric(unlist(groups_cluster))
for (group in groups_cluster) {
  cluster_dist = alltree$Get("dist", filterFun = isLeaf)
  #compute dist
  pock_dist = apply(dt[unlist(alltree$children[[group]]$pockets_names),
                       sort(colnames(dt))],1, function(x) {dist(rbind(new_pocket,x))})
  #dendrogram
  dt_centers.dend = alltree$children[[group]]$cluster_dend[[1]]#as.dendrogram(dt_centers.hclust)
  #
  centroid_green = which(labels(dt_centers.dend) == names(which(pock_dist == min(pock_dist))))
  lab_cluster_dend = labels(dt_centers.dend)
  labels(dt_centers.dend) <- paste0(lab_cluster_dend, 
                                    paste0(";Dist:",
                                           round(pock_dist[lab_cluster_dend], 2)))
  labels(dt_centers.dend)[centroid_green] <- paste0("POCHE->;",
                                                    paste0(lab_cluster_dend[centroid_green], 
                                                           paste0(";Dist:",
                                                                  round(min(pock_dist), 2))))

  c_colors = color.gradient(pock_dist)
  color.df<-data.frame(COLOR_VALUE=pock_dist, color.name=color.gradient(pock_dist))
  labels_colors(dt_centers.dend) <- as.character(color.df[lab_cluster_dend, "color.name"])
  # Open a PDF for plotting; units are inches by default
  pdf(paste0(paste0(paste0(paste0(args$o,
                                  args$name),"_grp_"),group),".pdf"),width=40, height=15)
  
  # Do some plotting
  par(cex=1, mar=c(13, 9, 5, 5))
  plot(dt_centers.dend, type = "rectangle", ylab = "Height", cex.lab = 2, cex.main = 3,
       main = paste0(paste0("Distance de la poche ",args$name),
                     " dans le groupes ", group))
  legend("right", title = "Distance de la poche",
         legend = c(paste0("max:", round(max(pock_dist),2)),
                    rep("",1),
                    paste0("mean:", round(mean(pock_dist),2)),
                    rep("",1),
                    paste0("min:", round(min(pock_dist),2))), pt.cex = 4, cex = 2, bty = "n")
  legend("topleft", title = "10 poches les plus proches:",
         legend = paste("P:", 
                        paste(names(sort(pock_dist)[1:10]),
                              paste("|D", round(sort(pock_dist)[1:10],2)))),
         pt.cex = 4, cex = 2, bty = "n")
  
  gradientLegend(valRange=c(-14,14), color=color.gradient(sort(pock_dist)), 
                 border.col=alpha('gray'), side=4, pos.num = 3, inside = FALSE, length=.2, depth=.02)#pos=c(120,3.5,125,5.5),
  dev.off()
}
print("DONE")
# ## DATA VISUALIZATION ##
# names_dt12 = c( "centers.p_polar_residues",
#                 "centers.p_Nlys_atom",
#                 "centers.p_aliphatic_residues",
#                 "centers.VOLUME_HULL",
#                 "centers.DIAMETER_HULL",
#                 "centers.p_Ooh_atom",
#                 "centers.hydrophobic_kyte",
#                 "centers.p_Ntrp_atom",
#                 "centers.p_Otyr_atom",
#                 "centers.C_RESIDUES",
#                 "centers.p_aromatic_residues",
#                 "centers.p_hydrophobic_residues")
# minmax_value = read.table(file = args[6])
# colnames(minmax_value) = names_dt12
# #get center values of first cluster for the minimum distance
# min_dist = min(alltree$Get("dist", filterFun = isLeaf))
# cluster_center_min_dist = alltree$Get(function(node) c(
#   centers.p_polar_residues = node$centers.p_polar_residues,
#   centers.p_Nlys_atom = node$centers.p_Nlys_atom,
#   centers.p_aliphatic_residues = node$centers.p_aliphatic_residues,
#   centers.VOLUME_HULL =  node$centers.VOLUME_HULL,
#   centers.DIAMETER_HULL =  node$centers.DIAMETER_HULL,
#   centers.p_Ooh_atom = node$centers.p_Ooh_atom,
#   centers.hydrophobic_kyte = node$centers.hydrophobic_kyte,
#   centers.p_Ntrp_atom = node$centers.p_Ntrp_atom,
#   centers.p_Otyr_atom = node$centers.p_Otyr_atom,
#   centers.C_RESIDUES = node$centers.C_RESIDUES,
#   centers.p_aromatic_residues = node$centers.p_aromatic_residues,
#   centers.p_hydrophobic_residues = node$centers.p_hydrophobic_residues),
#   filterFun = function(node) {
#     if(node$level == 2 && node$dist == min_dist) {
#       return(TRUE)
#     } else {
#       return(FALSE)
#     }
#   })
# dt.visualize = NULL
# dt.visualize = rbind(minmax_value, new_pocket)
# dt.visualize = rbind(dt.visualize, t(cluster_center_min_dist))
# 
# colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
# colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )
# radarchart(dt.visualize, axistype=1 , 
#            #custom polygon
#            pcol=colors_border ,plwd=4 , plty=1,
#            #custom the grid
#            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(-1,1,0.5), cglwd=0.8,
#            #custom labels
#            vlcex=0.8 
# )
# # Add a legend
# ###!!!!!!!!!!!!! ADD SAVING IMAGE IN OUTPUT FILE
# legend(x=0.6, y=1, legend = c("new_pocket", paste0("cluster:",best_cluster)), bty = "n", pch=20 , col=colors_in , cex=1, pt.cex=3)
# ###!!!!!!!!!
# ####HCLUST on closest cluster ####
# hclust_best = hclust(dist(rbind(dt[unlist(alltree$children[[best_cluster]]$pockets_names),sort(colnames(dt))],
#                                 new_pocket)),
#                      method = "ward.D2")
# 
# which(hclust_best$labels == "")
# hclust_best$labels[which(hclust_best$labels == "")] = "new_pocket"
# tail(hclust_best$labels)
# 
# hclust_best = as.dendrogram(hclust_best)
# 
# labels_colors(hclust_best) = 1
# labels_colors(hclust_best)["new_pocket"] = 2
# 
# hclust_best <- set(hclust_best, "labels_cex", 0.8)
# plot(hclust_best)
# 
# #apply dist
# start_time <- Sys.time()
# res_d = apply(dt, 1, function(x) {dist(rbind(new_pocket,x))})
# head(sort(res_d))
# end_time <- Sys.time()
# 
# end_time - start_time
# 
# start_time <- Sys.time()
# res_d_37 = apply(dt[unlist(alltree$children[[best_cluster]]$pockets_names),], 1, function(x) {dist(rbind(new_pocket,x))})
# head(sort(res_d),10)
# end_time <- Sys.time()
# 
# end_time - start_time
# 
# sapply(strsplit(names(head(sort(res_d),10)), "_"), "[", 2)
# head(sort(res_d),10)