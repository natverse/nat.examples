## This script assumed that you have run the file "05-miyasaka2014/00-setup.R", and the subsequent R files

# Setup
zm <- readRDS('neurons_canonical.rds')

# Convert neurons to dotprops objects
zmdps <- dotprops(zm, k=5, resample=0.1)
saveRDS(zmdps,file='zmdps.rds')
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

# Perform the same NBLAST, a little more smoothly
result.zfishob = nblast_allbyall(zmdps, normalisation = "mean")
pdf("images/nat_zfishOB_nblast_dendrogram.pdf", width = 15, height = 10)
hclust=nhclust(scoremat = result.zfishob)
k = 7
dkcs=colour_clusters(hclust, k = k)
plot(dkcs)
dev.off()

# Plot and colour in 3D
nopen3d(userMatrix = structure(c(0.999464690685272, 0.0153909111395478, 
                                 -0.0288552921265364, 0, -0.0159958377480507, 0.999654471874237, 
                                 -0.0208465103060007, 0, 0.0285244900733232, 0.0212973672896624, 
                                 0.999365568161011, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.746215641498566, 
        windowRect = c(1440L, 45L, 2942L, 1063L))
plot3d(hclust, db = zm, soma = TRUE, k = k, lwd = 2)
rgl.snapshot(filename = "images/nat_zfishOB_nblast_neurons.png", fmt = "png")
