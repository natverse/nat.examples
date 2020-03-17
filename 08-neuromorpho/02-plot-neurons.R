# set the scene
## This script assumed that you have run the file "08-neuromorpho/00-setup.R" and "01-download/00-setup.R"

## load data
load("Jacobs_principal_neurons.rda")
load("Jacobs_interneurons.rda")

## In order to see everything in thei relative sizes, we need to plot the bounding box of the biggest brain
bbx = nat::boundingbox(xyzmatrix(principals))
bbxxm = as.mesh3d(bbx)

## plot neurons
nopen3d(userMatrix = structure(c(0.999965906143188, 0.00817600265145302, 
                                 0.00111126154661179, 0, 0.00816878490149975, -0.999946534633636, 
                                 0.00633470714092255, 0, 0.00116293877363205, -0.00632557645440102, 
                                 -0.999979138374329, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 1, 
        windowRect = c(1440L, 45L, 2754L, 879L))
plot3d(bbxxm, alpha = 0, add = TRUE)
plot3d(principals, soma = TRUE)
rgl.snapshot(filename ="images/nat_neuromorpho_jacobs_principalneurons.png" ,fmt = "png")
clear3d()
plot3d(bbxxm, alpha = 0, add  = TRUE)
plot3d(interneurons, soma = TRUE)
rgl.snapshot(filename ="images/nat_neuromorpho_jacobs_interneurons.png" ,fmt = "png")

## re-plot the neurons by species
species = unique(principals[,"species"])
for(s in species){
  clear3d()
  plot3d(bbxxm, alpha = 0, add = TRUE)
  plot3d(subset(principals, species == s), soma = TRUE)
  rgl.snapshot(filename = paste0("images/nat_neuromorpho_jacobs_principalneurons_",s,".png") ,fmt = "png")
}



