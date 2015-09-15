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

# and plotting mirrored neurons, rather than originals
clear3d()
plot3d(hcdnm, db=dpnm, k=5, soma=8, lwd=2)
wire3d(monarch_cc.surf, col='grey', lwd=0.3)


# Now let's compare the clustering with and without mirroring
tanglegram(hcdn.dend, hcdnm.dend, margin_inner=10)

# Looking closely, we can see some differences and similarities. First the 
# similarities
# 1. dp-TuLAL.* are all in the same position
# 2. TL2a and TL2b alteranate properly, left then right. Furtheremore TL4 now joins them
# 3. TU-pdl alternate left/right and are all together
# now the main thing that needs attention is groups 2 and 4.

clear3d()
plot3d(hcdnm, db=dpn, k=5, soma=8, lwd=2,groups=c(2,4))

# Group 4 now consists of a pair of CPU1b neurons. These clearly belong together 
# the only question is whether they should be tagged as left or right. 
# They are tagged as left, presumably since they have dendritic and axonal 
# processes on the right. However they have their somata on the left as far as 
# I can see, so it seems that their projections are all contralateral.

# Now finally let's review group 2
ids_to_review=subset(hcdnm, k=5, groups=2)
# recluster
hcdnmr=nhclust(ids_to_review, scoremat=abam)
plot(hcdnmr)
# let's cut at four clusters
hcdnmr.dend=dendextend::color_branches(hcdnmr, k=4, col=rainbow)
dendextend::labels(hcdnmr.dend)=dpn[labels(hcdnmr.dend),'Detail.Page']
par(mar=c(8,5,5,2))
plot(hcdnmr.dend)

# plot in 3d without mirroring 
clear3d()
plot3d(hcdnmr, db=dpn, k=4, soma=8, lwd=2)

# ... and with mirroring
clear3d()
plot3d(hcdnmr, db=dpnm, k=4, soma=8, lwd=2)

# Hmm, mabye let's try k=5

hcdnmr.dend=dendextend::color_branches(hcdnmr, k=5, col=rainbow)
dendextend::labels(hcdnmr.dend)=dpn[labels(hcdnmr.dend),'Detail.Page']
par(mar=c(8,5,5,2))
plot(hcdnmr.dend)

clear3d()
plot3d(hcdnmr, db=dpn, k=5, soma=8, lwd=2)
clear3d()
plot3d(hcdnmr, db=dpnm, k=5, soma=8, lwd=2)


# So I have a lot of sympathy for the assignments made here

in particular subgroup 1 and 3
clear3d()
plot3d(hcdnmr, db=dpnm, k=5, soma=8, lwd=2, groups = c(1,3))

# the thing that looks most discordant is that there are a CPU1a type 1/2
# seem to be spread out (into groups 1,3,6 when cut at h=0.5). 

# Let's review in more detail:

hcdnmr.dend=dendextend::color_branches(hcdnmr, h=0.5, col=rainbow)
dendextend::labels(hcdnmr.dend)=dpn[labels(hcdnmr.dend),'Detail.Page']
par(mar=c(8,5,5,2))
plot(hcdnmr.dend)

# So looking at those groups 1,3,6
clear3d()
wire3d(monarch_cc.surf, col='grey', lwd=0.3)
plot3d(hcdnmr, db=dpnm, h=0.5, soma=8, lwd=2, groups = c(1, 3, 6))
snapshot3d("6.group2-h=0.5-subgroups-1-3-6.png")

# one can see that magenta neurons CL1b-L7 and CPU1a type1-R7 really are very 
# similar in everything but the lamination of their dendrites in CBU vs CBL

# This actually includes the axon terminal arborisation zones, which look very
# similar across these 3 plotted subgroups  -  even for the CL1b-L7 and CPU1a type1-R7.

# Conclusion: NBLAST clustering does not reveal all the subtypes
# according to the same scheme as human anatomists because it uses different
# invariants (in this case position across the whole neuron favouring staves of
# the central complex rather than lamination between the ellipsoid body vs
# fan shaped body.

# However the finest subtypes in this final dendrogram to co-cluster, e.g.
# CPU1a L5/R5/R4
# this suggests to me that if we had rather more data, NBLAST would indeed
# reliably co-cluster the finest level of cell type in this region.
