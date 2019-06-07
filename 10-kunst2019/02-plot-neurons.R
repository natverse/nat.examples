# set the scene
## This script assumed that you have run the file "10-kunst2019/00-setup.R", and the subsequent R files

# get data
fishmesh = read.hxsurf("data/zfish_simplified.surf")
load("zfish_neurons_right.rda")
load("zfish_neurons_original.rda")

# plot original
nopen3d(userMatrix = structure(c(-0.330627202987671, -0.943536818027496, 
                                 0.0206080116331577, 0, -0.940241873264313, 0.331199616193771, 
                                 0.0790761932730675, 0, -0.0814364328980446, 0.0067677479237318, 
                                 -0.996656060218811, 0, -38.2051658630371, 65.2883987426758, 8.2633752822876, 
                                 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, windowRect = c(1460L, 
                                                                                                 65L, 3039L, 1043L))
plot3d(fishmesh,alpha=0.1,col="lightgrey")
plot3d(zfish_neurons_original, lwd = 2, soma = TRUE)
rgl.snapshot(filename ="images/nat_zfish_neurons.png" ,fmt = "png")

# plot right
nopen3d(userMatrix = structure(c(-0.330627202987671, -0.943536818027496, 
                                 0.0206080116331577, 0, -0.940241873264313, 0.331199616193771, 
                                 0.0790761932730675, 0, -0.0814364328980446, 0.0067677479237318, 
                                 -0.996656060218811, 0, -38.2051658630371, 65.2883987426758, 8.2633752822876, 
                                 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, windowRect = c(1460L, 
                                                                                                 65L, 3039L, 1043L))
plot3d(fishmesh,alpha=0.1,col="lightgrey")
plot3d(zfish_neurons_right, lwd = 2, soma = TRUE)
rgl.snapshot(filename ="images/nat_zfish_neurons_right.png" ,fmt = "png")
