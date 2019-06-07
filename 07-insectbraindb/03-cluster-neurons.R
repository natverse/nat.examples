## This script assumed that you have run the file "07-insectbraindb/00-setup.R"

## First, we call the read neurons function, with ids set to NULL
insect.neurons = insectbraindb_read_neurons(ids = NULL)

## Subset neuronlist down to the Monarch neurons
butterfly.neurons = subset(insect.neurons, common_name == "Monarch Butterfly")

## convert to dotprops representation for nblast
### resample every 5µm since these neurons are big, also use 5 nearest neighbours
### to compute tangent vectors
### note use of progress bar since this is a bit slow!
message("converting neurons to dotprops for nblast")
butterfly.neurons.dps=dotprops(butterfly.neurons, resample=5, .progress='text', k=5)
butterfly.neurons.dps=butterfly.neurons.dps/5

## Calculating all by all nblast scores for these neurons
aba=nblast_allbyall(butterfly.neurons.dps, .progress='text')

## Plot our clustering
#### note use of the Neuron name as a label 
#### (fetched from the data.frame attached to butterfly.neurons neuronlist)
#### Looks like sensible relationships between neurons with similar names
hcdn=nhclust(scoremat=aba)
plot(hcdn, labels=with(butterfly.neurons, short_name ))
plot(hcdn, labels=with(butterfly.neurons, paste(short_name, hemisphere)))

## Plot colouring each cluster
hcdn.dend=color_branches(hcdn, k=5, col=rainbow)
#### nb it is necessary to permute the labels used for assignment so that they are
#### in dendrogram order 
dendextend::labels(hcdn.dend)=butterfly.neurons[labels(hcdn.dend),'short_name']

## Plot
par(mar=c(8,5,5,2))
plot(hcdn.dend)

## Plot the neurons in 3D
nat::nopen3d(userMatrix = structure(c(0.998683869838715, -0.00999264605343342, 
                                      -0.050302866846323, 0, -0.0184538997709751, -0.985155284404755, 
                                      -0.170669943094254, 0, -0.0478506684303284, 0.171373650431633, 
                                      -0.984043300151825, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.600000023841858, 
             windowRect = c(1440L, 45L, 3209L, 1063L))
plot3d(hcdn, db=butterfly.neurons, k=5, lwd=2, soma=5)

## Delevop a mapping of left to right for this brain region,  based on (http://onlinelibrary.wiley.com/doi/10.1002/cne.23054/abstract)
### to better compare neurons on left and right hemisphere of brain, or here the central complex.
#### Stanley Heinze says that CX, LX, and AOTU are based on the average template
#### here we have a simple affine transform computerd by surface registration in
#### Amira between the original surface and a mirror flipped version
m=structure(c(0.999891, 0.00954086, 0.0112588, 0, -0.00955665, 
            0.999953, 0.00134745, 0, -0.0112458, -0.00145435, 0.999935, 0, 
            3.36369, 0.210229, -0.0921294, 1), .Dim = c(4L, 4L))
### Now combine that with a mirror flip 
mm=m%*%scaleMatrix(-1,1,1)

## Apply the mirroring to all neurons marked as on the right
#### i.e. mapping everyone onto the left
butterfly.neurons.m=nat::xform(butterfly.neurons, reg=mm, subset=butterfly.neurons[,'hemisphere']=="Right")
#### and also to the dotprops versions of the neurons
butterfly.neurons.m.dps=nat::xform(butterfly.neurons.dps, reg=mm, subset=butterfly.neurons.dps[,'hemisphere']=="Right")
### Calculating all by all nblast scores for all neurons mapped to left
abam=nblast_allbyall(butterfly.neurons.m.dps, .progress='text')
hcdnm=nhclust(scoremat=abam)

## Plot!
plot(hcdnm, labels=with(butterfly.neurons, paste(short_name, hemisphere)))

## Now, let's colour the clusters
hcdnm.dend=color_branches(hcdnm, k=5, col=rainbow)
dendextend::labels(hcdnm.dend)=butterfly.neurons[labels(hcdnm.dend),'short_name']
par(mar=c(8,5,5,2))
plot(hcdnm.dend)
clear3d()
plot3d(hcdnm, db=butterfly.neurons, k=5, soma=8, lwd=2)

## And add the butterfly central complex
butterfly.brain = insectbraindb_read_brain(species = "Danaus plexippus")
central.complex = c("Anterior lip",
                    "Central body lower division",
                    "Central body upper division",
                    "Nodulus",
                    "Protocerebral bridge",
                    "Posterior optic tubercle")
central.complex.obj = butterfly.brain$RegionList[butterfly.brain$neuropil_full_names%in%central.complex]
plot3d(butterfly.brain,paste(central.complex.obj,collapse="|"), alpha = 0.5)

# and plotting mirrored neurons, rather than originals
clear3d()
plot3d(butterfly.brain,paste(central.complex.obj,collapse="|"), alpha = 0.1)
plot3d(butterfly.brain, alpha = 0.1)
rgl.snapshot(file = paste0("images/monarchmirror/CC_neurons.png"), fmt = "png")
npop3d()
plot3d(hcdn, db=butterfly.neurons, k=5, soma=5, lwd=2)
rgl.snapshot(file = paste0("images/monarchmirror/CC_neurons_NBLAST.png"), fmt = "png")
npop3d()
plot3d(hcdnm, db=butterfly.neurons.m, k=5, soma=5, lwd=2)
rgl.snapshot(file = paste0("images/monarchmirror/CC_neurons_mirrored_NBLAST.png"), fmt = "png")

## Now let's compare the clustering with and without mirroring
pdf("images/monarchmirror/CC_neurons_mirrored_NBLAST.pdf", width = 15, height = 10)
dendlist(hcdn.dend, hcdnm.dend) %>%
  untangle(method = "step1side") %>% 
  tanglegram(
    highlight_distinct_edges = FALSE, # Turn-off dashed lines
    common_subtrees_color_lines = FALSE, # Turn-off line colors
    common_subtrees_color_branches = FALSE, # Color common branches 
    margin_inner=5
  )
dev.off()

## Looking closely, we can see some differences and similarities. First the similarities
### 1. dp-TuLAL.* are all in the same position
### 2. TL2a and TL2b alteranate properly, left then right. Furtheremore TL4 now joins them
### 3. TU-pdl alternate left/right and are all together

## Now the main thing that needs attention is groups 2 and 4.
clear3d()
plot3d(hcdnm, db=butterfly.neurons, k=5, soma=8, lwd=2,groups=c(2,4))

## Group 4 now consists of a pair of CPU1b neurons. These clearly belong together 
## the only question is whether they should be tagged as left or right. 
## They are tagged as left, presumably since they have dendritic and axonal 
## processes on the right. However they have their somata on the left as far as 
## I can see, so it seems that their projections are all contralateral.

## Now finally let's review group 2
ids_to_review=subset(hcdnm, k=5, groups=2)
### recluster
hcdnmr=nhclust(ids_to_review, scoremat=abam)
plot(hcdnmr)
#### let's cut at four clusters
hcdnmr.dend=dendextend::color_branches(hcdnmr, k=4, col=rainbow)
dendextend::labels(hcdnmr.dend)=butterfly.neurons[labels(hcdnmr.dend),'short_name']
par(mar=c(8,5,5,2))
plot(hcdnmr.dend)

## plot in 3d without mirroring 
clear3d()
plot3d(hcdnmr, db=butterfly.neurons, k=4, soma=8, lwd=2)

## ... and with mirroring
clear3d()
plot3d(hcdnmr, db=butterfly.neurons.m, k=4, soma=8, lwd=2)

## Hmm, mabye let's try k=5
hcdnmr.dend=dendextend::color_branches(hcdnmr, k=5, col=rainbow)
dendextend::labels(hcdnmr.dend)=butterfly.neurons[labels(hcdnmr.dend),'short_name']
par(mar=c(8,5,5,2))
plot(hcdnmr.dend)

## And in 3D
clear3d()
plot3d(hcdnmr, db=butterfly.neurons, k=5, soma=8, lwd=2)
clear3d()
plot3d(hcdnmr, db=butterfly.neurons.m, k=5, soma=8, lwd=2)

## So I have a lot of sympathy for the assignments made here
### in particular subgroup 1 and 3
clear3d()
plot3d(hcdnmr, db=butterfly.neurons.m, k=5, soma=8, lwd=2, groups = c(1,3))

## The thing that looks most discordant is that there are a CPU1a type 1/2
## seem to be spread out (into groups 1,3,6 when cut at h=0.5). 

## Let's review in more detail:
hcdnmr.dend=dendextend::color_branches(hcdnmr, h=0.5, col=rainbow)
dendextend::labels(hcdnmr.dend)=butterfly.neurons[labels(hcdnmr.dend),'short_name']
par(mar=c(8,5,5,2))
plot(hcdnmr.dend)

## So looking at those groups 1,3,6
clear3d()
plot3d(hcdnmr, db=butterfly.neurons.m, h=0.5, soma=8, lwd=2, groups = c(1, 3, 6))

## one can see that magenta neurons CL1b-L7 and CPU1a type1-R7 really are very 
## similar in everything but the lamination of their dendrites in CBU vs CBL

## This actually includes the axon terminal arborisation zones, which look very
## similar across these 3 plotted subgroups  -  even for the CL1b-L7 and CPU1a type1-R7.

## Conclusion: NBLAST clustering does not reveal all the subtypes
#### according to the same scheme as human anatomists because it uses different
#### invariants (in this case position across the whole neuron favouring staves of
#### the central complex rather than lamination between the ellipsoid body vs
#### fan shaped body.

## However the finest subtypes in this final dendrogram to co-cluster, e.g.
### CPU1a L5/R5/R4
### this suggests to me that if we had rather more data, NBLAST would indeed
### reliably co-cluster the finest level of cell type in this region.
