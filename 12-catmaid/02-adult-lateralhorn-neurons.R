### This script assumes you have run "/12-catmaid/00-setup.R"

# Connect to the public FAFB instance (Zheng et al. 2018) hosted publicly by Virtual Fly Brain
adult.conn = catmaid_login(server="https://catmaid-fafb.virtualflybrain.org/")

# PD2a1/b1 neurons - for MBON connectivity (Dolan et al. 2018)
em.pd2 = fetchn_fafb("name:^PD2a1/b1", mirror = F, reference = JFRC2)
em.pd2[,"type"] = "ON"

# Assign cell info
names(em.pd2) = sub(" .*", "", em.pd2[,"name"])
em.pd2[,"cell"] = names(em.pd2)
lhn_types=process_lhn_name(em.pd2[,'cell'])
em.pd2[,]=cbind(em.pd2[,], lhn_types[,c("pnt",'anatomy.group',"cell.type")])
em.pd2[,"compartment"] = "whole"
em.pd2[,"glom"] = "PD2a1/b1"
em.pd2.flow = catnat::flow.centrality(em.pd2, polypre= FALSE, mode = "centrifugal")
em.pd2.flow.dendrites = dendritic_cable(em.pd2.flow)
em.pns = readRDS("/GD/LMBD/Papers/2014bridgingtemplates/fig/EMLM/data/neurons/EM_PNs.rds")
em.pns[,"glomerulus"] = em.pns[,"glom"]
em.pns[,"data.type"] = "EM"
em.pns = subset(em.pns,tract == "mALT"&glom!="mPN")

# Show the lateral horn neuron PD2a1 split into axon and dendrite, and save a .png of it
nopen3d(userMatrix = structure(c(0.994485974311829, 0.0881130397319794, 
                                 0.056865967810154, 0, 0.094024084508419, -0.509014368057251, 
                                 -0.855607032775879, 0, -0.0464447215199471, 0.856236219406128, 
                                 -0.514492392539978, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.530321657657623, 
        windowRect = c(1440L, 45L, 3215L, 1052L))
neurons = flow.centrality(em.pd2["1299700"], mode = "centrifugal", polypre = FALSE)
seesplit3d(neurons,lwd=4,soma=TRUE, col = c("blue", "orange", "purple", "chartreuse",
                                            "grey", "pink"),WithConnectors = FALSE)
rgl.snapshot(filename = paste0("images/PD2a1_Split.png"), fmt = "png")

# Now with synapses
clear3d()
plot3d(neurons, col = "darkgrey", lwd = 4 , soma = TRUE)
rgl::spheres3d(nat::xyzmatrix(get.synapses(neurons,target = "POST")), col = 'cyan', radius = 0.5)
rgl::spheres3d(nat::xyzmatrix(get.synapses(neurons,target = "PRE")), col = 'red', radius = 0.75)
rgl.snapshot(filename = paste0("images/PD2a1_Synapses.png"), fmt = "png")

# Show PD2a1 split by microtubules
clear3d()
mt = prune_microtubules(neurons[1][[1]], microtubules = TRUE)
twigs = prune_microtubules(neurons[1][[1]], microtubules = FALSE)
rgl::plot3d(mt, col = "burlywood4", WithNodes = FALSE, soma = TRUE, 
            WithConnectors = FALSE, lwd = 4)
rgl::plot3d(twigs, col = "burlywood1", WithNodes = FALSE, 
            WithConnectors = FALSE, lwd = 4)
rgl.snapshot(filename = "images/PD2a1_microtubules.png", fmt = "png")

# extract just the spine
clear3d()
sp=spine(neurons[[1]])
plot3d(sp, col = "purple", soma = FALSE, lwd = 4)
plot3d(neurons, col = "lightgrey", soma = TRUE, lwd = 4)
rgl.snapshot(filename = "images/PD2a1_spine.png", fmt = "png")

# Plot the branch points and end nodes
clear3d()
points = nat::xyzmatrix(neurons)
epoints = nat::xyzmatrix(neurons)[nat::endpoints(neurons[[1]]),]
bpoints = nat::xyzmatrix(neurons)[nat::branchpoints(neurons[[1]]),]
rpoints = nat::xyzmatrix(neurons)[nat::rootpoints(neurons[[1]]),]
plot3d(neurons, lwd = 3, col = "green")
spheres3d(epoints, col = "blue", radius = 0.2)
spheres3d(bpoints, col = "red", radius = 0.2)
spheres3d(rpoints, col = "purple", radius = 2)
rgl.snapshot(filename = "images/PD2a1_neuron_object.png", fmt = "png")
