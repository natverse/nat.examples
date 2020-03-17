# set the scene
## This script assumed that you have run the file "14-mouselight/00-setup.R"

## First we can quickly just plot the outer mesh for the brain
outline = mouselight_read_brain(type = "outline")
plot3d(outline, col = "pink", alpha = 0.3)

## This is cool, but maybe what we really want are its sub-divisions.
mousebrain = mouselight_read_brain(type = "brain_areas")
clear3d()
plot3d(mousebrain)
### This takes a long time the first time you call this function per session
### .obj files are saved in a temporary folders
### You could save this brain locally, to use quickly in future.
#### save(mousebrain,"YOUR_PATH/mouselight_brain.rda")

## What brain regions are on offer?
print(mousebrain$neuropil_full_names)

## Or if we want more information, we can get it like this:
mbr = mouselight_brain_region_info()
View(mouselight_brain_region_info)

## Perhaps we want to plot just the amygdala?
### To do this we can do
amygdala.codes = mousebrain$RegionList[grepl("amygdala",mousebrain$neuropil_full_names,
ignore.case = TRUE)]
plot3d(outline, col = "pink", alpha = 0.1)
plot3d(subset(mousebrain, amygdala.codes), alpha = 0.5)

## So now we want some neurons
## What neurons data is available?
ndf=mouselight_neuron_info()
#### How many tracings per neurons?
table(table(ndf$neuron.id))
#### This is because many 'neurons' have a separate axon and dendrite skeleton

## We can download all of these neurons, and their meta-data
### Typically two tracings, and axon and a dendrite, per neuron
mlns = mouselight_read_neurons(ndf$tracing.id, meta = TRUE)

## Let's read in all the amygdalal neurons
### Since each
in.amyg = mouselight_nodes_in_region(mlns, brain.areas = amygdala.codes, labels = NULL)
amyg.ids = names(in.amyg)[in.amyg>0]
amyg.neurons = mlns[amyg.ids]

## And plot!
plot3d(amyg.neurons)
### Quite wide-ranging!


