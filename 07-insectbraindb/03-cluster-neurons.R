# Use NBLAST to cluster neurons by structure

library(nat.nblast)

load("tedore_neurons.rda")

# subset neuronlist down to the Monarch neurons
dpn=subset(tedoren, Species=="Danaus plexippus")
dpn[,'side']=ifelse(grepl("(right|R[0-9]+)", dpn[,'Detail.Page']),"R","L")

# convert to dotprops representation for nblast
# resample every 5µm since these neurons are big
# note use of progress bar since this is a bit slow!
message("converting neurons to dotprops for nblast")
dpn2=dotprops(dpn, resample=5, .progress='text')
dpn2=dpn2/5

message("calculating all by all nblast scores for these neurons")
aba=nblast_allbyall(dpn2, .progress='text')

# plot our clustering
# note use of the Neuron name as a label 
# (fetched from the data.frame attached to dpn neuronlist)
# Looks like sensible relationships between neurons with similar names
hcdn=nhclust(scoremat=aba)
plot(hcdn, labels=with(dpn, Neuron ))
plot(hcdn, labels=with(dpn, paste(Neuron, side)))
plot(hcdn, labels=with(dpn, Detail.Page ))

# plot colouring each cluster
library(dendroextras)
hcdn2=color_clusters(hcdn, k=5)
labels(hcdn2)=with(dpn, paste(Neuron, side))
par(mar=c(5,8,5,2))
plot(hcdn2)

clear3d()
plot3d(hcdn, db=dpn, k=5, soma=5)

# TODO
# Delevop a mapping of left to right for this brain region
# likely based on this
# http://onlinelibrary.wiley.com/doi/10.1002/cne.23054/abstract
# to better compare neurons on left and right side of brain.

