# set the scene
## This script assumed that you have run the file "03-helmstaedter2013/00-setup.R" and "03-helmstaedter2013/01-download.R"

# load the skeleton object we saved earlier
load("sk.uniq.rda")

# turn skel objects into a SWC type structure
skel_to_neuron <- function(x){
  y = as.data.frame(x$nodes)
  colnames(y) = c("X", "Y", "Z", "W")
  y$Parent = -1
  y[x$edges[,2],"Parent"] = x$edges[,1]
  nat::as.neuron(y)
}
reroot_neuron <- function (x) {
  x = nat::as.neuron(nat::as.ngraph(x), origin = which(x$d$Label == 
                                                         4))
  x$d$Label = 0
  x$d$Label[x$StartPoint] = 1
  x
}

# make sure all its elements are neuron objects
skn = nlapply(sk.uniq, skel_to_neuron, OmitFailures = TRUE, .progress='text')

# save neurons in native R format
save(skn, file='skn.rda')

# save a zip archive of SWC format neurons for all reconstructions
write.neurons(skn, dir='skuniq.swc.zip', files=names(skn), format='swc')


neuron(d = s$nodes, )
