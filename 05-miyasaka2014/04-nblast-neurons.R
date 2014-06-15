# Setup
zm <- readRDS('neurons_canonical.rds')
library(nat.nblast)
if(!require('stringr')) install.packages("stringr")
if(!require('dendroextras')) install.packages("dendroextras")


# Convert neurons to dotprops objects
zmdps <- dotprops(zm, k=5, resample=1)
zm_class <- sub("-", "", str_match(names(zm), "[A-z]*-T?L?"))


# NBLAST neurons
zm_res <- nblast(zmdps, zmdps)
zm_res_norm <- scale(zm_res, center=FALSE, scale=diag(zm_res))
zm_res_mean <- (zm_res_norm + t(zm_res_norm)) / 2


# Discard self-matches
diag(zm_res_norm) <- NA
diag(zm_res_mean) <- NA


# Convert NBLAST results to lists and sort
zm_res_list <- apply(zm_res_norm, 2, list)
zm_res_list <- lapply(zm_res_list, function(x) sort(x[[1]], decreasing=TRUE))
zm_res_list <- zm_res_list[order(sapply(zm_res_list, function(x) x[[1]]), decreasing=TRUE)]

zm_res_mean_list <- apply(zm_res_mean, 2, list)
zm_res_mean_list <- lapply(zm_res_mean_list, function(x) sort(x[[1]], decreasing=TRUE))
zm_res_mean_list <- zm_res_mean_list[order(sapply(zm_res_mean_list, function(x) x[[1]]), decreasing=TRUE)]


# Convert scores to distances
zm_res_mean_dist <- 1 - zm_res_mean
diag(zm_res_mean_dist) <- 0
zm_res_mean_dist <- as.dist(zm_res_mean_dist)
zm_clust <- hclust(zm_res_mean_dist)
zm_dend <- as.dendrogram(zm_clust)

zm_cols <- c(maG="#3635fd", dG="#f08380", vaG="#df20be", lG="#6e9d3f", vpG="#6a6d26", mdGT="#87b3d8", mdGL="#dea76f", vmG="#454545")
par(mar=c(3,1,1,5))
zm_dend_col=set_leaf_colors(zm_dend, structure(zm_cols[zm_class],.Names=names(zm)), col_to_set='label')
zm_dend_col=color_clusters(zm_dend_col,k=4,col=rep('black',4),groupLabels=as.roman)
plot(zm_dend_col, horiz=TRUE)


# Helper display functions
show_hits <- function(neuron_name, num_hits=5) {
  plot3d(zm[neuron_name], lwd=2, col='black', soma=T)
  plot3d(zm[names(zm_res_list[[neuron_name]][1:num_hits])], soma=T)
}

show_hits_mean <- function(neuron_name, num_hits=5) {
  plot3d(zm[neuron_name], lwd=2, col='black', soma=T)
  plot3d(zm[names(zm_res_mean_list[[neuron_name]][1:num_hits])], soma=T)
}
