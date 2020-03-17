### This script assumes you have run "/00-basics/00-setup.R"

# Let's have a go at plotting a neuron from EM data!!
## Reconstructed by Dolan et al. 2018, Neuron and Zheng et al. 2018, Cell, in a CATMAID instance?
## Connect to the public FAFB instance (Zheng et al. 2018) hosted publicly by Virtual Fly Brain
adult.conn = catmaid_login(server="https://catmaid-fafb.virtualflybrain.org/")

## Get the neuron in question
nopen3d(userMatrix = structure(c(0.99775356054306, -0.0018308530561626, 
                                 -0.0669662356376648, 0, -0.0177687853574753, -0.97105461359024, 
                                 -0.238195270299911, 0, -0.0645917728543282, 0.238850027322769, 
                                 -0.968905925750732, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.75, 
        windowRect = c(20L, 65L, 1402L, 887L))
n.em = da1.em[[1]]
n.em = unspike(n.em, threshold = 5)
n.em.lh = prune_in_volume(n.em, brain = FCWBNP.surf, neuropil = "LH.R") # Cut to justone neuropil sub-volume!
plot3d(subset(FCWBNP.surf,"LH_R"), alpha = 0.2, col = "grey") # Plot that sub-volume!
plot3d(as.neuronlist(n.lh), soma = FALSE, col = "chartreuse4", lwd  = 4)
rgl.snapshot(filename = "images/nat_basic_da1_LH_axon.png", fmt = "png")

# we can also put it in the context of the whole neuron, from which this branch was cut
plot3d(n.em.lh, soma = FALSE, WithConnectors = TRUE, col = "black", lwd = 4)
rgl.snapshot(filename = "/GD/LMBD/Papers/2014bridgingtemplates/fig/EMLM/illustrator/images/basic_nat/nat_basic_da1_LH_axon_EM.png", fmt = "png")

# And then with the whole template brain plotted as well!
clear3d()
plot3d(subset(FCWBNP.surf,"LH_R"), alpha = 0.2, col = "green")
plot3d(n.em.lh, soma = FALSE, WithConnectors = TRUE, col = "black", lwd = 4)
rgl.snapshot(filename = "images/basic_nat/nat_basic_LH_axon_EM.png", fmt = "png")

# Good stuff.
