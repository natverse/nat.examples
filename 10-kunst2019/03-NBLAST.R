# set the scene
## This script assumed that you have run the file "10-kunst2019/00-setup.R", and the subsequent R files

# get data
fishmesh = read.hxsurf("data/zfish_simplified.surf")
load("zfish_neurons_right.rda")

# show the result of an NBLAST
zfish_neurons_right.dps = dotprops(zfish_neurons_right, stepsize = 1, OmitFailures = TRUE)
result.zafish = nblast_allbyall(zfish_neurons_right.dps, normalisation = "mean")
pdf("images/nat_zfish_nblast_dendrogram.pdf", width = 5, height = 3)
hclust=nhclust(scoremat = result.zafish)
k = 50
dkcs=colour_clusters(hclust, k = k) %>% set("branches_lwd", 1.5)
plot(dkcs, leaflab = "none")
dev.off()

# show neurons coloured by NBLAST
nopen3d(userMatrix = structure(c(-0.330627202987671, -0.943536818027496, 
                                 0.0206080116331577, 0, -0.940241873264313, 0.331199616193771, 
                                 0.0790761932730675, 0, -0.0814364328980446, 0.0067677479237318, 
                                 -0.996656060218811, 0, -38.2051658630371, 65.2883987426758, 8.2633752822876, 
                                 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, windowRect = c(1460L, 
                                                                                                 65L, 3039L, 1043L))

plot3d(fishmesh,alpha=0.1,col="lightgrey")
plot3d(hclust, db = zfish_neurons_right, soma = TRUE, k = k)
rgl.snapshot(filename = "images/nat_zfish_nblast_neurons.png", fmt = "png")

# Show the individual clusters
cols = rainbow(k)
for(i in 1:k){
  clear3d()
  plot3d(fishmesh,alpha=0.1,col="lightgrey")
  plot3d(hclust, db = zfish_neurons_right, soma = TRUE, groups = i, k = k, col = cols[i])
  rgl.snapshot(filename = paste0("images/nat_zfish_nblast_neurons_",i,".png"), fmt = "png")
}

## So part of the advtange here is that all neurons are mirrored to the right
### So that neurons that originated on either hemisphere, but are of the same cell type
### can be classed together. Let's look at a quick example of that
clear3d()
load("zfish_neurons_original.rda")
r = c("T_16119_13_3", "T_161019_8_1", "T_161116_7_4", "T_160502_4_2_corrected", 
      "T_150707_14_3_full", "T_150910_2_3c_full_tracing", "T_160406_4_1", 
      "T_160520_9_4", "T_150714_5_3b_full_tracing", "20170117_1013_BGUG_HuC_ltRFP_d5_F4", 
      "20170901_1013_BGUG_HuC_ltRFP_d6_F1_Neuron2", "20161005_1013_BGUG_HuC_ltRFP_d6_F11_Neuron2", 
      "20160920_1013_BGUG_HuC_ltRFP_d6_F11", "20170506_1013_BGUG_HuC_ltRFP_d5_F3_Neuron2", 
      "20170805_1013_BGUG_HuC_ltRFP_d7_F4", "T_180102_20_1", "T_161109_1_1", 
      "T_16119_7_1_1", "T_150714_5_3a_full_tracing", "20170608_1013_BGUG_HuC_ltRFP_d6_F1", 
      "T_16119_20_4_CHECK")
l = c("T_160524_13_1", "T_160705_17_2", "T_161031_12_4", "T_161019_8_4", 
      "T_170529_11_3", "T_161116_3_2", "T_160527_3_1", "T_16119_9_2", 
      "T_151130_4_2", "20170313_1013_BGUG_HuC_ltRFP_d5_F8", "20170228_1013_BGUG_HuC_ltRFP_d5_F1", 
      "T_160920_11_2", "20170818_1013_BGUG_HuC_ltRFP_d6_F3", "T_150714_3_3_full_tracing")
plot3d(fishmesh,alpha=0.1,col="lightgrey")
plot3d(zfish_neurons_original[l], col = "cyan", lwd = 2, soma = TRUE)
plot3d(zfish_neurons_original[r], col = "red", lwd = 2, soma = TRUE)
rgl.snapshot(filename = paste0("images/nat_zfish_superclusters_apart.png"), fmt = "png")
clear3d()
plot3d(fishmesh,alpha=0.1,col="lightgrey")
plot3d(zfish_neurons_right[l], col = "cyan", lwd = 2, soma = TRUE)
plot3d(zfish_neurons_right[r], col = "red", lwd = 2, soma = TRUE)
rgl.snapshot(filename = paste0("images/nat_zfish_superclusters_mirrored.png"), fmt = "png")

