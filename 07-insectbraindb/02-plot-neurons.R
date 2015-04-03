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
plot3d(subset(tedoren, Species=="Danaus plexippus"))
