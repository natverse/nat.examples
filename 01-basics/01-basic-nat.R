### This script assumes you have run "/00-basics/00-setup.R"

# Let's have a go at plotting some basic nat data types. 
## Here, we will look at olfactory projection neurons in the vinegar fly, Drosophila melanogaster
## These data originate from the FlyCircuit.tw project, and we are accessing them using
## the lhns r package, developed for Frechter et al. 2019, eLife.
## Images will be saved in the images/ subfolder of this directory

# Get uniglomerular ACh neurons
AL_mALT_PN1 = subset(lhns::most.lhins, anatomy.group == "AL-mALT-PN1")
AL_mALT_PN1[,"data.type"] = "light"
da1s = subset(AL_mALT_PN1, glomerulus == "DA1")
n = da1s[[2]]
## The DA1 PN is one of the most studied neurons in the fly brain.
## It responds to fly sex pheromones.

# Image of light-level neuron
points = nat::xyzmatrix(n)
epoints = nat::xyzmatrix(n)[nat::endpoints(n),]
bpoints = nat::xyzmatrix(n)[nat::branchpoints(n),]
rpoints = nat::xyzmatrix(n)[nat::rootpoints(n),]
nopen3d(userMatrix = structure(c(0.99775356054306, -0.0018308530561626, 
                                 -0.0669662356376648, 0, -0.0177687853574753, -0.97105461359024, 
                                 -0.238195270299911, 0, -0.0645917728543282, 0.238850027322769, 
                                 -0.968905925750732, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.55, 
        windowRect = c(20L, 65L, 1402L, 887L))
plot3d(n, lwd = 3, col = "green")
spheres3d(epoints, col = "blue", radius = 0.5)
spheres3d(bpoints, col = "red", radius = 0.75)
spheres3d(rpoints, col = "darkgreen", radius = 2)
rgl.snapshot(filename = "images/basic_nat/LM_DA1_PN.png", fmt = "png")

# Great stuff, but maybe we also want to see a neuron from a different data source
# What about an EM neuron? Reconstructed by Zheng et al. 2018, Cell, in a CATMAID instance?
## Connect to the public FAFB instance (Zheng et al. 2018) hosted publicly by Virtual Fly Brain
adult.conn = catmaid_login(server="https://catmaid-fafb.virtualflybrain.org/")
## Get the neuron in question
da1.em = fetchn_fafb("name:PN glomerulus DA1", mirror = FALSE, reference = FCWB)[1:5] # Bridge to the same space as our light-level neuron!!
n.em = da1.em[[1]]
n.em = unspike(n.em, threshold = 100) # correct data abberation fromn mis-registrations

# Image of light-level neuron
points = nat::xyzmatrix(n.em)
epoints = nat::xyzmatrix(n)[nat::endpoints(n),]
bpoints = nat::xyzmatrix(n)[nat::branchpoints(n),]
rpoints = nat::xyzmatrix(n)[nat::rootpoints(n),]
nopen3d(userMatrix = structure(c(0.99775356054306, -0.0018308530561626, 
                                 -0.0669662356376648, 0, -0.0177687853574753, -0.97105461359024, 
                                 -0.238195270299911, 0, -0.0645917728543282, 0.238850027322769, 
                                 -0.968905925750732, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.55, 
        windowRect = c(20L, 65L, 1402L, 887L))
plot3d(n.em, lwd = 3, col = "green", WithConnectors = TRUE)
rgl.snapshot(filename = "images/EM_DA1_PN.png", fmt = "png")

# Hmm, cool. But let's also get an example of a mammalia neuron. Also an olfactory projection; a mitral cell.
# Let's get our  mammmalian neuron by searching the NeuroMorpho database.
mitral.df = neuromorpho_search(search_terms = c("brain_region:main olfactory bulb"))
mitral.cells = neuromorpho_read_neurons(neuron_name = mitral.df$neuron_name, batch.size = 5, nat = TRUE, progress = TRUE, OmitFailures = TRUE)
mitral.cells = mitral.cells[mitral.cells[,"species"]=="mouse"]
mitral.cells = mitral.cells[grepl("mitral",mitral.cells[,"cell_type"])]
mt  = mitral.cells["86520"][[1]]

# Image of cortical neuron
points = nat::xyzmatrix(mt)
epoints = nat::xyzmatrix(mt)[nat::endpoints(mt),]
bpoints = nat::xyzmatrix(mt)[nat::branchpoints(mt),]
rpoints = nat::xyzmatrix(mt)[nat::rootpoints(mt),]
nopen3d(userMatrix = structure(c(0.512131154537201, -0.81005847454071, 
                                 0.285529524087906, 0, 0.840945184230804, 0.540534257888794, 0.0251813232898712, 
                                 0, -0.174736797809601, 0.22721853852272, 0.958039045333862, 0, 
                                 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.556837499141693, windowRect = c(1440L, 
                                                                                                          45L, 2875L, 959L))
plot3d(mt, lwd = 3, col = "green")
spheres3d(epoints, col = "blue", radius = 1.5)
spheres3d(bpoints, col = "red", radius = 3)
spheres3d(rpoints, col = "darkgreen", radius = 20)
rgl.snapshot(filename = "images/basic_nat/nat_basic_neuron_object_MC.png", fmt = "png")

# We have plotted data from different sources. Now, it might be nice to see how different data types are plotted.

# plot neuron object
nopen3d(userMatrix = structure(c(0.99775356054306, -0.0018308530561626, 
                                 -0.0669662356376648, 0, -0.0177687853574753, -0.97105461359024, 
                                 -0.238195270299911, 0, -0.0645917728543282, 0.238850027322769, 
                                 -0.968905925750732, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.55, 
        windowRect = c(20L, 65L, 1402L, 887L)) # set window view
plot3d(n, soma = TRUE, lwd = 4)
plot3d(FCWB, alpha = 0.1, col = "lightgrey")
plot3d(subset(FCWBNP.surf,"LH_R"), alpha = 0.2, col = "grey")
rgl.snapshot(filename = "images/nat_basic_da1_neuron.png", fmt = "png")

# plot neuronlist
clear3d()
plot3d(AL_mALT_PN1, soma = TRUE)
plot3d(FCWB, alpha = 0.1, col = "lightgrey")
plot3d(subset(FCWBNP.surf,"LH_R"), alpha = 0.2, col = "grey", lwd = 2)
rgl.snapshot(filename = "images/nat_basic_upns_neuronlist.png", fmt = "png")

# Cool. But what if we only want a bit of a neuron? We can prune it!!
## Cut neuron
clear3d()
n.lh = prune_in_volume(n, brain = FCWBNP.surf, neuropil = "LH.R")
plot3d(n.lh, soma = FALSE)
plot3d(as.neuronlist(n), soma = TRUE, col = "darkgrey", lwd  =4)
plot3d(FCWB, alpha = 0.1, col = "lightgrey")
plot3d(subset(FCWBNP.surf,"LH_R"), alpha = 0.2, col = "grey")
rgl.snapshot(filename = "/GD/LMBD/Papers/2014bridgingtemplates/fig/EMLM/illustrator/images/basic_nat/nat_basic_da1_subset.png", fmt = "png")


# We might also want to use a popular analysis, looking at Strahler order (https://en.wikipedia.org/wiki/Strahler_number)
## Decompose into Strahler orders
clear3d()
so=strahler_order(n)
orders=1:max(so$points)
cols = rainbow(max(orders))
for (i in orders) {
  som = ifelse(i == max(orders),TRUE,FALSE)
  plot3d(as.neuronlist(subset(n, so$points==i)), col= cols[i], soma = som, lwd  = 4)
}
plot3d(as.neuronlist(n), col = "lightgrey", soma = FALSE, lwd = 4)
rgl.snapshot(filename = "images/basic_nat/nat_basic_da1_strahler.png", fmt = "png")

# Or we can just take the longest path through a neuron
## So let's extract just the spine
clear3d()
sp=spine(n)
plot3d(sp, col = "purple", soma = FALSE, lwd = 4)
plot3d(as.neuronlist(n), col = "lightgrey", soma = TRUE, lwd = 4)
rgl.snapshot(filename = "images/nat_basic_da1_spine.png", fmt = "png")

# Ooh, an let's plot it just as a graph structure. Note, this is 2D plot!
pdf("images/nat_basic_da1_graph.pdf", width = 10, height = 10)
sg=segmentgraph(n)
p <- plot(sg, layout = layout_with_fr, edge.arrow.size=.1, vertex.color="green", vertex.size=5, 
     
     vertex.frame.color="gray", vertex.label.color="black", 
     
     vertex.label = "", vertex.label.cex=0, vertex.label.dist=5, edge.curved=0.2) 
p
dev.off()
p

