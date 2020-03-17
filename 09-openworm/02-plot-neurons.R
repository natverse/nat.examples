# set the scene
## This script assumed that you have run the file "08-openworm/00-setup.R"

## get data
celegans.neurons = read.neurons("data/generatedNeuroML/")

## plot
celegans.neurons = celegans.neurons[!grepl("L$",names(celegans.neurons))]
points  = xyzmatrix(celegans.neurons)
wormesh = as.mesh3d(alphashape3d::ashape3d(points, alpha = 30))
celegans.neurons = celegans.neurons[!grepl("L$",names(celegans.neurons))]
nopen3d(userMatrix = structure(c(0.755564987659454, 0.106279626488686, 
                                 -0.646394729614258, 0, -0.629285573959351, 0.391895711421967, 
                                 -0.671131134033203, 0, 0.181991666555405, 0.91385018825531, 0.362982928752899, 
                                 0, -4.64235009001107, -33.9248777045344, -2.67488689509321e-06, 
                                 1), .Dim = c(4L, 4L)), zoom = 0.431433469057083, windowRect = c(1460L, 
                                                                                                 65L, 3254L, 1045L))
plot3d(wormesh,alpha=0.1,col="lightgrey", add = TRUE)
plot3d(celegans.neurons, lwd = 2, soma = TRUE)
rgl.snapshot(filename ="images/nat_roundworm_neurons.png" ,fmt = "png")
