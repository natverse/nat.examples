## This script assumed that you have run the file "05-miyasaka2014/00-setup.R", and the subsequent R files

zm <- readRDS('neurons.rds')

# Mirror all neurons past middle of x range
mid_x <- 248.857/2
root_xs<-sapply(zm,function(x) x$d$X[x$StartPoint])
neurons_to_mirror <- names(root_xs[root_xs > mid_x])

# Have a look at what we want to mirror
nopen3d(userMatrix = structure(c(0.999464690685272, 0.0153909111395478, 
                                 -0.0288552921265364, 0, -0.0159958377480507, 0.999654471874237, 
                                 -0.0208465103060007, 0, 0.0285244900733232, 0.0212973672896624, 
                                 0.999365568161011, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.746215641498566, 
        windowRect = c(1440L, 45L, 2942L, 1063L))
plot3d(zm[neurons_to_mirror], gene=='lhx2a',  lwd=2, soma=T, col = "cyan")
plot3d(zm, gene=='lhx2a',  lwd=2, soma=T, col = "red")
rgl.snapshot(filename ="images/nat_zfishOB_neurons_leftright.png" ,fmt = "png")
clear3d()
zm.canonical[neurons_to_mirror] <- mirror(zm[neurons_to_mirror], warpfile="mirroring_registration.list/", mirrorAxisSize=248.857)
plot3d(zm.canonical[neurons_to_mirror], gene=='lhx2a',  lwd=2, soma=T, col = "cyan")
plot3d(zm.canonical, gene=='lhx2a',  lwd=2, soma=T, col = "red")
rgl.snapshot(filename ="images/nat_zfishOB_neurons_mirror.png" ,fmt = "png")

# Save results
saveRDS(zm.canonical, file="neurons_canonical.rds")


