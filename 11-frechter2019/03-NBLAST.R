# set the scene
## This script assumed that you have run the file "11-frechter2019/00-setup.R"

# Show the result of an NBLAST
opns = subset(lhns::most.lhins, grepl("Olf", modality))
opns.dps = dotprops(opns, stepsize = 1, OmitFailures = TRUE)
result = nblast_allbyall(opns.dps, normalisation = "mean")
pdf("images/nat_adultfly_nblast_dendrogram.pdf", width = 5, height = 3)
hclust=nhclust(scoremat = result)
k = 30
dkcs=colour_clusters(hclust, k = k) %>% set("branches_lwd", 2)
plot(dkcs, leaflab = "none")
dev.off()

# Show neurons coloured by NBLAST
nopen3d(userMatrix = structure(c(0.905018091201782, 0.231813445687294, 
                                 -0.356658041477203, 0, -0.0334783419966698, -0.797041535377502, 
                                 -0.602996289730072, 0, -0.424053758382797, 0.557662665843964, 
                                 -0.713576078414917, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, 
        windowRect = c(12L, 44L, 1398L, 904L))
plot3d(hclust, db = opns, soma = TRUE, k = k)
plot3d(FCWB,alpha=0.1,col="lightgrey")
rgl.snapshot(filename = "images/nat_adultfly_nblast_neurons.png", fmt = "png")

# Show the individual clusters
cols = rainbow(k)
for(i in 1:k){
  clear3d()
  plot3d(FCWB,alpha=0.1,col="lightgrey")
  plot3d(hclust, db = opns, soma = TRUE, groups = i, k = k, col = cols[i])
  rgl.snapshot(filename = paste0("images/nat_adultfly_nblast_neurons_",i,".png"), fmt = "png")
}
