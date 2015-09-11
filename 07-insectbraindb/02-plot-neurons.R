library(nat)
# load in cached version of processed neurons from disk
load("tedore_neurons.rda")

# what kind of neurons do we have?
head(tedoren)

# the with.neuronlist method gives you convenient access to the attached metadata frame
# here we summarise the "Species" for each neuron
with(tedoren, summary(Species))

# plot neurons from specific structural cluster
nopen3d()
plot3d(subset(tedoren, Species=="Danaus plexippus"), lwd=2, soma=5)

## add surface object for monarch central complex

# newest versions of rgl package can read wavefront obj files.
if(exists("readOBJ", 'package:rgl')){
  if(!file.exists('monarch_brain_optimized_hh1SWlD.obj'))
    download.file("http://www.tedore.net/media/brain-models/5/monarch_brain_optimized_hh1SWlD.obj", 
                  destfile = 'monarch_brain_optimized_hh1SWlD.obj')
  monarch_cc.surf=rgl::readOBJ("monarch_brain_optimized_hh1SWlD.obj")
} else {
  # get a cached copy from flybrain website
  download.file("http://flybrain.mrc-lmb.cam.ac.uk/si/nblast/nat.examples/07-insectbraindb/monarch_cc.rds", "monarch_cc.rds")
  monarch_cc.surf=readRDS("monarch_cc.rds")
}
wire3d(monarch_cc.surf, col='grey', lwd=0.5)
