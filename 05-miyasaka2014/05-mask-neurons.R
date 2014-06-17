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
    rfun<-function(x, y=NULL, z=NULL) with(xyz.coords(x,y,z), b[1,1]<=x & x<=b[2,1] & b[1,2]<=y & y<=b[2,1] & b[1,3]<=z & z<=b[2,3])
}


zmdps=readRDS('zmdps.rds')
zm.mask=nlapply(zmdps, function(x) subset(x, bbox2sel3d(mask.box.bbox)))
zm.mask.vTel=nlapply(zmdps, function(x) subset(x, bbox2sel3d(mask.box.vTel.bbox)))
saveRDS(zm.mask,file='zmdps.mask.rds')
saveRDS(zm.mask.vTel,file='zm.mask.vTel.rds')
