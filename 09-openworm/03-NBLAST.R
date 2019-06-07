# set the scene
## This script assumed that you have run the file "08-openworm/00-setup.R"

# get data
celegans.neurons = read.neurons("data/generatedNeuroML/")

# Show the result of an NBLAST
celegans.neurons = celegans.neurons[!grepl("L$",names(celegans.neurons))]
pdf("images/nat_roundworm_nblast_dendrogram.pdf", width = 15, height = 10)
celegans.neurons.dps = dotprops(celegans.neurons, stepsize = 1, OmitFailures = TRUE)
result = nblast_allbyall(celegans.neurons.dps, normalisation = "mean")
hclust=nhclust(scoremat = result)
dkcs=colour_clusters(hclust, k= 10)
plot(dkcs)
dev.off()

# Show neurons coloured by NBLAST
clear3d()
plot3d(hclust, db = celegans.neurons, soma = TRUE, k = 10)
rgl.snapshot(filename = "images/nat_roundworm_nblast_neurons.png", fmt = "png")
