# Use NBLAST to cluster neurons by structure

library(nat.nblast)

load("tedore_neurons.rda")

# subset neuronlist down to the Monarch neurons
dpn=subset(tedoren, Species=="Danaus plexippus")
dpn[,'side']=ifelse(grepl("(right|R[0-9]+)", dpn[,'Detail.Page']),"R","L")

# convert to dotprops representation for nblast
# resample every 5µm since these neurons are big, also use 5 nearest neighbours
# to compute tangent vectors
# note use of progress bar since this is a bit slow!
message("converting neurons to dotprops for nblast")
dpn2=dotprops(dpn, resample=5, .progress='text', k=5)
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
library(dendextend)
hcdn.dend=color_branches(hcdn, k=5, col=rainbow)
# nb it is necessary to permute the labels used for assignment so that they are
# in dendrogram order 
dendextend::labels(hcdn.dend)=dpn[labels(hcdn.dend),'Detail.Page']

par(mar=c(8,5,5,2))
plot(hcdn.dend)

clear3d()
plot3d(hcdn, db=dpn, k=5, soma=5)

# TODO
# Delevop a mapping of left to right for this brain region
# likely based on this
# http://onlinelibrary.wiley.com/doi/10.1002/cne.23054/abstract
# to better compare neurons on left and right side of brain.

# OK, let's try with monarch_cc.surf'
# Stanley Heinze says that CX, LX, and AOTU are based on the average template

# here we have a simple affine transform computerd by surface registration in
# Amira between the original surface and a mirror flipped version
m=structure(c(0.999891, 0.00954086, 0.0112588, 0, -0.00955665, 
            0.999953, 0.00134745, 0, -0.0112458, -0.00145435, 0.999935, 0, 
            3.36369, 0.210229, -0.0921294, 1), .Dim = c(4L, 4L))
# Now combine that with a mirror flip 
mm=m%*%scaleMatrix(-1,1,1)

# apply the mirroring to all neurons marked as on the right
# i.e. mapping everyone onto the left
dpnm=nat::xform(dpn, reg=mm, subset=dpn[,'side']=="R")
# and also to the dotprops versions of the neurons
dpn2m=nat::xform(dpn2, reg=mm, subset=dpn2[,'side']=="R")
message("calculating all by all nblast scores for all neurons mapped to left")
abam=nblast_allbyall(dpn2m, .progress='text')
hcdnm=nhclust(scoremat=abam)

plot(hcdnm, labels=with(dpn, paste(Neuron, side)))

# now, let's colour the clusters
hcdnm.dend=color_branches(hcdnm, k=5, col=rainbow)
dendextend::labels(hcdnm.dend)=dpn[labels(hcdnm.dend),'Detail.Page']
par(mar=c(8,5,5,2))
plot(hcdnm.dend)
clear3d()
plot3d(hcdnm, db=dpn, k=5, soma=8, lwd=2)
wire3d(monarch_cc.surf, col='grey', lwd=0.3)

