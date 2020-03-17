## This script assumed that you have run the file "05-miyasaka2014/00-setup.R", and the subsequent R files

# Miyasaka et al present analyse the projection patterns in two domains
# vTel and dTel

# find bounding boxes
mask.box=read.im3d("mask-box.nrrd",ReadByteAsRaw=TRUE)
mask.box.coords=imexpand.grid(mask.box)[mask.box>0,]
mask.box.bbox=apply(mask.box.coords,2,range)
mask.box.vTel=read.im3d("mask-box-vTel.nrrd",ReadByteAsRaw=TRUE)
mask.box.vTel.coords=imexpand.grid(mask.box.vTel)[mask.box.vTel>0,]
mask.box.vTel.bbox=apply(mask.box.vTel.coords,2,range)

# convert bounding box to a selection function
bbox2sel3d<-function(b) {
    rfun<-function(x, y=NULL, z=NULL) with(xyz.coords(x,y,z), b[1,1]<=x & x<=b[2,1] & b[1,2]<=y & y<=b[2,2] & b[1,3]<=z & z<=b[2,3])
}

# subset every neuron in zmdps neuronlist down to respective bounding box
zmdps=readRDS('zmdps.rds')
zm.mask=nlapply(zmdps, function(x) subset(x, bbox2sel3d(mask.box.bbox)))
zm.mask.vTel=nlapply(zmdps, function(x) subset(x, bbox2sel3d(mask.box.vTel.bbox)))

# Plot in 3D
nopen3d(userMatrix = structure(c(0.999464690685272, 0.0153909111395478, 
                                 -0.0288552921265364, 0, -0.0159958377480507, 0.999654471874237, 
                                 -0.0208465103060007, 0, 0.0285244900733232, 0.0212973672896624, 
                                 0.999365568161011, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.746215641498566, 
        windowRect = c(1440L, 45L, 2942L, 1063L))
plot3d(zm.mask, lwd  =2, soma = FALSE)
plot3d(zm, soma = TRUE, col = "lightgrey", lwd = 1)
rgl.snapshot(filename = "images/nat_zfishOB_nblast_neurons_dTel.png", fmt = "png")
clear3d()
plot3d(zm.mask.vTel, lwd  = 2, soma = FALSE)
plot3d(zm, soma = TRUE, col = "lightgrey", lwd = 1)
rgl.snapshot(filename = "images/nat_zfishOB_nblast_neurons_vTel.png", fmt = "png")

# Save
saveRDS(zm.mask,file='zmdps.mask.rds')
saveRDS(zm.mask.vTel,file='zm.mask.vTel.rds')
