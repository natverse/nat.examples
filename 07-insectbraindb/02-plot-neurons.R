## This script assumes that you have run the files
# "07-insectbraindb/00-setup.R"
# "07-insectbraindb/01-download-neurons.R"

## And let's plot them
nat::nopen3d(userMatrix = structure(c(0.998683869838715, -0.00999264605343342, 
                                      -0.050302866846323, 0, -0.0184538997709751, -0.985155284404755, 
                                      -0.170669943094254, 0, -0.0478506684303284, 0.171373650431633, 
                                      -0.984043300151825, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.600000023841858, 
             windowRect = c(1440L, 45L, 3209L, 1063L))
plot3d(butterfly.neurons, lwd = 2, soma = 5)

## Cool! But maybe we also want to see it's template brain? 
## Let's check if they have it
available.brains = insectbraindb_species_info()
available.brains

## Great, they do, let's get it
butterfly.brain = insectbraindb_read_brain(species = "Danaus plexippus")

## And plot as a wireframe
wire3d(butterfly.brain, col='gold')
rgl.snapshot(file = paste0("images/MonarchButterly_ReppertLab_.png"), fmt = "png")

## Oop, that's a lot of neuropils. 
## Let's go for only a subset. What's available?
butterfly.brain$RegionList
butterfly.brain$neuropil_full_names

## There lateral horn (LH) and the antennal lobe (AL) are my favourites.
## Let's plot those
clear3d()
plot3d(butterfly.neurons, lwd = 2, soma = 5)
plot3d(subset(butterfly.brain, "LH|AL"), alpha = 0.3)

## Sometimes it's useful to make simplified versions of neurons / meshes

# the neurons seem to be rather densely sampled
steplengths=unlist(sapply(butterfly.neurons,seglengths, sumsegment=F))
# (for details, ?seglengths)
hist(steplengths, br=400, xlim=c(0,1), col='red')

# Let's resample to 1µm which is probably fine for visualisation
butterfly.neurons.rs1=resample(butterfly.neurons, stepsize = 1)
nclear3d()
plot3d(butterfly.neurons.rs1, lwd = 2, soma = 5)

# Plotting all neuropils in one go with alpha is pretty slow.
# Instead we can combine them into a single mesh
wholebrain.hires=as.mesh3d(butterfly.brain)
# subsample to 10% of number of nodes
wholebrain=Rvcg::vcgQEdecim(wholebrain.hires, percent = 0.1)
# add vertex normals for smoother looking plot
wholebrain=addNormals(wholebrain)
# now plot
shade3d(wholebrain, col='gold', alpha=.1)

