# This script assumed that you have run the file "15-hemibrain/00-setup.R"

## What do we want to do with all this hemibrain data?
### Perhaps we could look at some popular neuron types and see if they connect to each other?
### Mushroom body output neurons (MBONs) are a well-studied neuron class.
### They read-out stored associative memories, 'located' in the mushroom body of the flies.
### These mushroom bodies are split into compartments (e.g. a, alpha or y, gamma)
### That have different learning properties.
### Let's find some in our HemiBrain data:
mbon.info = neuprint_search(".*MBON.*")
mbon.info[!is.na(mbon.info$type),"type"] = mbon.info[!is.na(mbon.info$type),"name"]
#### Add MB compartment information from the name
name.split = strsplit(mbon.info$name,split="\\(")
comps = sapply(name.split, function(x) tryCatch(gsub("\\)|_|bi.*|R.*|L.*|\\?.*|P.*","",x[[2]]), error = function(e) ""))
mbon.info$compartment = comps
mbon.info = subset(mbon.info, compartment!="")

## Let's see what meta-info we get, along with stuff with 'MBON' in the name
View(mbon.info)
### Looks like we get their unique ID numbers (bodyIds), cell names,
### types, their size (voxels) and synapse numbers (pre=output, post=input).

## Let's quickly read one neuron, and have a look at it!
mbon = neuprint_read_neuron(mbon.info$bodyid[1])
### visualise:
nopen3d(userMatrix = structure(c(0.947127282619476, 0.222770735621452, 
                                 -0.230919197201729, 0, 0.247805655002594, -0.0506941750645638, 
                                 0.967482566833496, 0, 0.203820571303368, -0.973551869392395, 
                                 -0.103217624127865, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.58467960357666, 
        windowRect = c(1460L, 65L, 2877L, 996L)) # set view
plot3d(mbon, WithConnectors = TRUE, col = lacroix[["brown"]], lwd = 2)
rgl.snapshot(filename = "images/hemibrain_an_mbon.png", fmt ="png")
### output synapses in red, input ones in blue
### You can see the auto-segmentation goes a bit crazy where the soma is

## Now that we have all of the bodyIds, we can read these neurons from neuPrint:
mbons = neuprint_read_neurons(mbon.info$bodyid)
clear3d()
plot3d(mbons, col = sample(lacroix,length(mbons)), lwd = 2)
### Why is this function so slow?
### Because it fetches fragmented neuron skeletons assigned to the same neuron
### Stitched them into one neuron, grabs all the neuron's synapses
### and tries to work out where the soma is. If you want to move faster
### You can grab bits separately and quickly:
?neuprint_read_neuron_simple
?neuprint_get_synapses
?neuprint_assign_connectors
?neuprint_locate_soma

## Do we have to know what's in the name to find interesting neurons? No!
### We can also read neurons from a region of interest (ROI) of the brain.
### The regions of interest available are:
rois = sort(neuprint_ROIs())
rois

## We don't we try fetching all of the neurons in a compartment of the mushroombody lobes?
### The compartment a'3 is interesting, because it is involved in aversive odour learning
### and novetly detection (Aso et al. 2014; Hattori et al. 2016)
ap3.info = neuprint_find_neurons(
  input_ROIs = "a'3(R)",
  output_ROIs = NULL,
  all_segments = FALSE # if true, fragments smaller than 'neurons' are returned as well
)

## So how many neurons is that?
print(nrow(ap3.info))

## A lot! But what about that pesky load of Kenyon cells
### The most populous neurons in the brain
ap3.info = subset(ap3.info, !is.na(bodyname) & neuronStatus == "Traced")
ap3.info = ap3.info[!grepl("KC",ap3.info$bodyname),]

## So how many neurons is that?
print(nrow(ap3.info))

## Getting all that data took a while, so let's save it for later
save(mbons,"hemibrain_mbons.rda")

