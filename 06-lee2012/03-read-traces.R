message("Reading neurons ...")
load("flycircuit.md.rda")
kcraw=read.neurons('tracings',patt="lineset\\.am\\.gz$",
  neuronnames=function(f) sub("_lineset.*","",basename(f)),
  df=flycircuit.md)

# flip neurons (Z and X) if required
####
needflipdf=subset(kcraw, reverseandflip, rval='data.frame')
# nmapply is a specialised mapply for neuronlists
flipz=nmapply(mirror,kcraw[needflipdf$gene_name],mirrorAxis='Z',mirrorAxisSize=needflipdf$Z2)
kcraw[needflipdf$gene_name]=nmapply(mirror,flipz,mirrorAxis='X',mirrorAxisSize=needflipdf$X2)

## alternatively
# for(n in needflipdf$gene_name){
#   zflip=mirror(kcraw[[n]],mirrorAxis='Z',mirrorAxisSize= needflipdf[n,'Z2'])
#   kcraw[[n]]=mirror(zflip,mirrorAxis='X',mirrorAxisSize= needflipdf[n,'X2'])
# }

# scale neurons to original physical coords
# i.e how they would have appeared in the confocal microscope stack
####
df=attr(kcraw, 'df')
kcscaled=kcraw
for(n in names(kcscaled)){
  kcscaled[[n]]=with(df[n,],kcscaled[n]*c(VoxelSizeX/dx, VoxelSizeY/dy, VoxelSizeZ/dz))
}

####
# => these are now ready for registrations to be applied
# however those registrations are unfortunately not yet available for download.

# Therefore as a simple expedient, centre them according to their centroid
kccent=nlapply(kcscaled, function(n) n-colMeans(xyzmatrix(n) ))

# some of them are from the fly's right, let's flip
# would have identified neurons to flip interactively by making a selection box
# in rgl window.
# plot3d(kccent)
# right=find.neuron(db=kccent)
right=c("ChaMARCM-F000586_seg002", "FruMARCM-F000085_seg001", "FruMARCM-F000188_seg001", 
"FruMARCM-F000706_seg001", "FruMARCM-F001115_seg002", "FruMARCM-F001494_seg002", 
"FruMARCM-M001051_seg002", "FruMARCM-M001339_seg001", "GadMARCM-F000050_seg001", 
"GadMARCM-F000071_seg001", "GadMARCM-F000423_seg001", "GadMARCM-F000442_seg002", 
"GadMARCM-F000476_seg001")
left=setdiff(names(kccent),right)

# flip the ones on fly's right and recentre everybody
kccent[right]=kccent[right]*c(-1,1,1)
kccent=nlapply(kccent, function(n) n-colMeans(xyzmatrix(n) ))

# ok, so now these neurons are coarsely aligned
plot3d(kccent)