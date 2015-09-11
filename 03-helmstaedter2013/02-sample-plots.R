load("sk.uniq.rda")
load("skeleton_metadata.rda")

# All neurons, coloured by raw type id
plot3d(sk.uniq,col=factor(typeid))
# add (uncoloured) somata
spheres3d(t(t(skeleton_metadata$kn.e2006.ALLSKELETONS.FINAL2012.allSomata[,1:3])*c(16.5,16.5,25)),rad=2000)

# All neurons, coloured by raw type id
plot3d(sk.uniq,col=factor(typeid))

# All neurons, coloured by sorted type id
plot3d(sk.uniq,col=factor(stypeid))

# just bipolar cells, types 58-71
plot3d(sk.uniq,class=='bipolar',col=factor(stypeid))

# just bipolar cells, types 58-71
plot3d(sk.uniq,class=='bipolar',col=factor(stypeid))

# Glial and starburst amacrines - potential coordinate systems
# TODO Check "glial" identity
clear3d();plot3d(sk.uniq,class=='glial',col='grey')


clear3d()
plot3d(sk.uniq,class=='bipolar',col='grey')
plot3d(sk.uniq,ntype=='starburst amacrine',col=stypeid)
plot3d(sk.uniq,class=='ganglion',col='green')
bb=boundingbox(sk.uniq)
plot3d(bb)
