library(nat)
# load in cached version of processed neurons from disk
load("sumbuln.rda")

# plot neurons from specific structural cluster
plot3d(subset(sumbuln, cluster=="Z"))
# two clusters
plot3d(subset(sumbuln, cluster%in%c("A","Z")), col=cluster)
# add decorations
axes3d()
title3d(xlab='X', ylab='Y', zlab='IPL depth')

# plot everybody except GFP neurons
# start to see some clear stratification
plot3d(subset(sumbuln, geneticLine!="GFP-YFP"), col=rainbow(nlevels(geneticLine))[geneticLine])

# the with.neuronlist method gives you convenient access to the attached metadata frame
# here we summarise the "cluster" letter for each neuron
with(sumbuln, summary(cluster))
# ... and the genetic line
with(sumbuln, summary(geneticLine))
# here we plot the cell body positions as spheres
# coloured by the geneticLine used
with(sumbuln,
  spheres3d(X,Y,Z, radius=4, col=rainbow(nlevels(geneticLine))[geneticLine])
)