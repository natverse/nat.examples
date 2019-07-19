# Get our neurons of interest
load("Jacobs_principal_neurons.rda")
laod("Jacobs_interneurons.rda")

# Do a sholl type analysis
sholl_analysis <- function(neuron, start = colMeans(xyzmatrix(neuron)), starting.radius = radius.step, ending.radius = NULL, radius.step = ending.radius/100){
  unit.vector <- function(x) {x / sqrt(sum(x^2))}
  dend = neuron$d
  dend$dists = nabor::knn(data = matrix(start,ncol=3), query = nat::xyzmatrix(neuron),k=1)$nn.dists
  if(is.null(ending.radius)){
    ending.radius = max(dend$dists)
  }
  radii = seq(from = starting.radius, to = ending.radius, by = radius.step)
  sholl = data.frame(radii = radii, intersections = 0)
  for(n in 1:length(radii)){
    r = radii[n]
    segments = neuron$SegList
    for(segment in segments){
      p = dend[segment,]
      dists = (nabor::knn(data = matrix(start,ncol=3), query = nat::xyzmatrix(p),k=1)$nn.dists - r) >= 0
      sholl[n,]$intersections = sholl[n,]$intersections + lengths(regmatches(paste(dists,collapse=""), gregexpr("TRUEFALSE|FALSETRUE", paste(dists,collapse=""))))
    }
  }
  sholl
}
