# Here, we want to test out some data that the Akerman group in Oxford has derived from confocal
# imaging of Layer 5 pyramidal cells, and processed through simple neurite 
# tracer and other data extraction tools in FIJI. (skeletonisation tool in Amira
# 6.4 is probably much more time efficient)
# Data courtesy of A. Vourvoukelis, A. von Klemperer and C.J. Akerman.

# First of all, we need to load the packages we will be using
# Make sure the location of the present R file is your working directory
# See getwd()
setwd(here::here("00-setup"))

source("02-packages.R")

# Now let's read in that .swc file
neuron = read.neuron("data/test/axon_traces_27_06_19_s4_c7_Pom.traces")
stepsize = 0.1
neuron = nat::resample(neuron, stepsize = stepsize)

# And let's also read the bouton file
boutons = read.csv("data/test/27_06_19_s4_c7_Pom_boutons.csv")

# Let's plot our neuron. 
# We could set soma to TRUE, or a valu like 100, and this would plot a sphere at the root
# of this skeleton. However, I do not know if the root of the skeleton is actually the
# soma? In AMIRA there is an easy way to root the skeleton at a point. It's good if 
# you can do that for the soma location. Otherwise, you can use correct_soma function in
# the nat package to update the soma/root of a skeleton interactively, in R.
plot3d(as.neuronlist(neuron), col = "black", soma = FALSE, lwd = 4)

# And now we can plot the boutons
points3d(boutons[,c("X","Y","Slice")], col = "red")

# So the boutons seem (systematically?) displace from the skeleton? 
# Seems like an error.
# Let's snap them back to the skeleton
points = boutons[,c("X","Y","Slice")]
near = nabor::knn(query = points, data = xyzmatrix(neuron),
                  k = 1)$nn.idx
points.new = xyzmatrix(neuron)[near,]
boutons[,c("X","Y","Z")] = points.new
boutons$treenode = near

# Plot the new boutons
points3d(boutons[,c("X","Y","Z")], col = "cyan", size = 5)

# Or, we can move the skeleton towards the boutons
trans = computeTransform(x =points, y  = points.new, type = "affine")
## And then snap the boutons in
near = nabor::knn(query = points, data = xyzmatrix(neuron),
                  k = 1)$nn.idx
points.new.2 = xyzmatrix(neuron)[near,]
boutons[,c("X","Y","Z")] = points.new.2
boutons$treenode = near

# Plot the new boutons
xyzmatrix(neuron) <- applyTransform(xyzmatrix(neuron), trans)
plot3d(as.neuronlist(neuron), col = "green", soma = FALSE, lwd = 4)
points3d(boutons[,c("X","Y","Z")], col = "magenta", size = 5)

# So now we can quikly work out the boutons per cable length
info = summary(neuron)
print(info)
boutons.per.cable = nrow(boutons)/info$cable.length

# So 0.05 boutons per micron of cable. We could also see branchpoints
# per cable
bps.per.cable = info$branchpoints/info$cable.length

# That's 0.006, so about 10 x more boutons per cable than branch points.
# Could branchpoints be another interesting mesure? Or leaf nodes? I.e. endpoints?

# What does the break down look like on a per branch basis?
boutons$seg = sapply(boutons$treenode, function(x)
  which(sapply(neuron$SegList, function(y) x%in%y))[1])
t = table(boutons$seg)
print(t)

# So it's clear some branches have more boutons than others
# But is this helpful? This depends on how the skeletonisation
# procedure has been done. A better measure, may be Strahler order
# NOTE: for this to work properly, the root of the neuron MUST be the soma.
neuron = assign_strahler(neuron)

# What has this done? We can take a quick look
clear3d()
plot3d(as.neuronlist(neuron), col = "black", soma = FALSE, lwd = 1)
points3d(xyzmatrix(neuron), col = neuron$d$strahler_order)

# Here, the branches have been coloured by Strahler order
# What about the boutons
boutons$strahler = neuron$d$strahler_order[boutons$treenode]
t2 = table(boutons$strahler)
print(t2)

# So it is clear that most of the boutons are on Strahler order 1 branches,
# i.e. the leaf branches of the neuron

# To get the number of branchpoints per cable for these three orders
# We need to get the cable length. This could be doen by just counting the
# Number of points, if we resample the neuron to one point per micron
# Luckily, we did that at the start, so let's just calculate it
t3 = table(neuron$d$strahler_order)
# calculate boutons per micron
t2[setdiff(names(t3),names(t2))] = 0
print(t2/t3)

# So the boutons per micron is different by Strahler order

# We can also run a sholl analysis, another popular measure for mammalian neurons
# Becuase people haven't thought of a better analysis for mammalian neurons somehow...

# First, let us write a fucntion 
# perform a sholl type analysis
sholl_analysis <- function(neuron, start = colMeans(xyzmatrix(neuron)), 
                           starting.radius = radius.step, ending.radius = 1000, 
                           radius.step = ending.radius/100){
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

# Run this function on our neuron
sholla = sholl_analysis(neuron)

# plot the result
ggplot(sholla, aes(x=radii, y=intersections, color="darkgrey")) +
  geom_line()+
  theme_minimal()

# We could also plot the distance to boutons on this graph I guess...
root = neuron$d[rootpoints(neuron),c("X","Y","Z")]
near2 = nabor::knn(data = root, query = xyzmatrix(boutons),
                   k = 1)$nn.dist
bdist = data.frame(intersections = 0, dist = near2)

# Now plot again!
ggplot(sholla, aes(x=radii, y=intersections, color="darkgrey")) +
  geom_line()+
  geom_density(data = bdist, aes(x=dist, y=intersections, color="red"))+
  theme_minimal()

# Him, this visualisation is a little hard to deal with, could also do a density plot
ggplot(data = bdist, aes(x=dist, color="red", fill = "red")) + 
  geom_density()+
  theme_minimal()

# Okay, cool. So now we need to work out which layers these boutons are in
## for this, Valex has a high-res image (from which he got the neuron data) and a low-res one (2D, for the layers)
### first name is the lower layer from border, second is the upper layer from border (greater Y value)
surface = read.csv("data/test/layer_lines_csv/surface_MAX_27_06_19_s4_c7_Pom.csv")[,c("X","Y")]
surface$lboundary = "l1_surface"
l23 = read.csv("data/test/layer_lines_csv/L23_L1_MAX_27_06_19_s4_c7_Pom.csv")[,c("X","Y")]
l23$lboundary = "l23_l1"
l4 = read.csv("data/test/layer_lines_csv/L4_L23_MAX_27_06_19_s4_c7_Pom.csv")[,c("X","Y")]
l4$lboundary = "l4_23"
l5 = read.csv("data/test/layer_lines_csv/L5_L4_MAX_27_06_19_s4_c7_Pom.csv")[,c("X","Y")]
l5$lboundary = "l5_l4"
l6 = read.csv("data/test/layer_lines_csv/L6_L5_MAX_27_06_19_s4_c7_Pom.csv")[,c("X","Y")]
l6$lboundary = "l6_l5"
wm = read.csv("data/test/layer_lines_csv/wm_L6_MAX_27_06_19_s4_c7_Pom.csv")[,c("X","Y")]
wm$lboundary = "wm_l6"
## wm is white matter

# All together now
layers = rbind(surface,
               l23,
               l4,
               l5,
               l6,
               wm)
layers$lower = gsub("_.*","",layers$lboundary)
layers$upper = gsub(".*_","",layers$lboundary)

# So these coordinates are 2D, and in pixels. Let's convert to mictons
layers[,c("X","Y")] = layers[,c("X","Y")]/2.4089

# And lets add a Z dimension
layers$Z = colMeans(xyzmatrix(neuron))["Z"]

# Okay great, so we have some layers. How does it look?
clear3d()
plot3d(as.neuronlist(neuron), col = "black", lwd = 3)
points3d(boutons[,c("X","Y","Z")], col = "darkred", size = 5)
layer.cols = rainbow(length(unique(layers$lower)))
names(layer.cols) = unique(layers$lower)
points3d(xyzmatrix(layers), col = layer.cols[layers$lower])

# So now we need to assign the boutons to layers
## to assign the boutons, we need to find the closest point, and see if it is above or below the layer line
near = nabor::knn(query = xyzmatrix(boutons), data = xyzmatrix(layers), k = 1)$nn.idx
boutons$nearest.lboundary = layers[near,"lboundary"]
boutons$nearest.lboundary.x = layers[near,"X"]
boutons$layer = ifelse(boutons$nearest.lboundary.x>boutons$X, 
                       gsub("_.*","",boutons$nearest.lboundary),
                       gsub(".*_","",boutons$nearest.lboundary))

# Now plot again, colouring the 
clear3d()
plot3d(as.neuronlist(neuron), col = "black", lwd = 3)
spheres3d(boutons[,c("X","Y","Z")], col = layer.cols[boutons$layer], radius = 10)
points3d(xyzmatrix(layers), col = layer.cols[layers$lower])

# Great!! Now we can take an rgl.snapshot to save it!!
open3d(userMatrix = structure(c(0.286752492189407, 0.947574555873871, 
                                0.140977531671524, 0, 0.956106007099152, -0.292328834533691, 
                                0.0201254710555077, 0, 0.0602807253599167, 0.129018723964691, 
                                -0.989808678627014, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.814054071903229, 
       windowRect = c(-41L, 108L, 1383L, 955L))
plot3d(as.neuronlist(neuron), col = "black", lwd = 3)
spheres3d(boutons[,c("X","Y","Z")], col = layer.cols[boutons$layer], radius = 10)
points3d(xyzmatrix(layers), col = layer.cols[layers$lower])
rgl.snapshot(filename = "images/test_neuron_layers_bouton.png", fmt = "png")

# So we want to look at the the bouton density per layer, which means we need to know which treenodes are in which layer
## let's do the same again then
near2 = nabor::knn(query = xyzmatrix(neuron), data = xyzmatrix(layers), k = 1)$nn.idx
neuron$d$nearest.lboundary = layers[near2,"lboundary"]
neuron$d$nearest.lboundary.x = layers[near2,"X"]
neuron$d$layer = ifelse(neuron$d$nearest.lboundary.x>neuron$d$X, 
                       gsub("_.*","",neuron$d$nearest.lboundary),
                       gsub(".*_","",neuron$d$nearest.lboundary))

# Calculate
node.layers = table(neuron$d$layer)
boutons.layers = table(boutons$layer)[names(node.layers)]
bouton.density = boutons.layers / node.layers

# Make data frame for plotting
df = rbind(node.layers,
           boutons.layers,
           bouton.density)
df = melt(df)
colnames(df) = c("type", "layer", "value")

# Plot, and save plots!
pdf("images/nodes_per_layer_bar_chart.pdf", height = 10, width = 5)
ggplot(subset(df, grepl("node.layers",type)), aes(x=type, y=value, fill=layer)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = layer.cols)+
  theme_minimal()
dev.off()

pdf("images/boutons_per_layer_bar_chart.pdf", height = 10, width = 5)
ggplot(subset(df, grepl("boutons.layers",type)), aes(x=type, y=value, fill=layer)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = layer.cols)+
  theme_minimal()
dev.off()

pdf("images/bouton_density_per_layer_bar_chart.pdf", height = 10, width = 5)
ggplot(subset(df, grepl("bouton.density",type)), aes(x=type, y=value, fill=layer)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = layer.cols)+
  theme_minimal()
dev.off()


# So we also have a midline for the pyramidal cell that we can read and look at
axis = read.csv("data/test/Axis_line/Axis_line_27_06_19_s4_c7.csv")[,c("X","Y")]
axis[,c("X","Y")] = axis[,c("X","Y")]/2.4089
axis$Z = colMeans(xyzmatrix(neuron))["Z"]

# Euclidean bouton distances to the midline
bouton.points = boutons[,c("X","Y","Slice")] # let's use the original bouton locations for this
boutons$dist.to.midline = unlist(nabor::knn(query = bouton.points, data = axis, k = 1)$nn.dists)

# Plot!
pdf("images/boutons_euclidean_distance_to_midline.pdf", height = 10, width = 5) # 50 um bins
ggplot(boutons, aes(x=dist.to.midline, color = layer, fill = layer)) + 
  geom_histogram(binwidth=50) +
  scale_fill_manual(values = layer.cols)+
  scale_color_manual(values = layer.cols)+
  theme_minimal()
dev.off()

# Now lets say we want to calculate the geodesic 
## It is VERY VERY important that the function resample has been run on the neuron
## before we start, i.e. at the top of this file. Stepsize = 1 means our each point
## on our neuron is about 1 um from the next point in the tree.
## to calculate geodesic distance, we are just going to count the number of points
## that separate boutons.

### And let's work this out per layer
graph = as.ngraph(neuron) # we can now use this representation with functions from the package iGraph.
## A directed graph!

# A single value per layer
boutons$interbouton.distance  = NA
interbouton = data.frame()
for(l in unique(boutons$layer)){
  # Just boutons on one layers
  bouton.nodes = unique(subset(boutons, layer == l)$treenode)
  # So we want to take the mean distance for the bouton behind, and in front, of each bouton of interest
  ## the distance to the next bouton, away from the root.
  dists1 = distances(graph, v = bouton.nodes, to = bouton.nodes, mode = c("out"), weights = NULL)
  diag(dists1) = NA
  dist.out = apply(dists1, 1, function(x) min(x, na.rm = TRUE)) # node that cannot be reached have an infinite value
  dist.out = dist.out[!is.infinite(dist.out)]
  ## And towards it
  dists2 = distances(graph, v = bouton.nodes, to = bouton.nodes, mode = c("in"), weights = NULL)
  diag(dists2) = NA
  dist.in = apply(dists2, 1, function(x) min(x, na.rm = TRUE)) # node that cannot be reached have an infinite value
  dist.in = dist.in[!is.infinite(dist.in)]
  interbouton = rbind(interbouton, data.frame(
    layer = l,
    mean = mean( c(dist.out,dist.in) ) / stepsize,
    sd = sd( c(dist.out,dist.in) ) / stepsize,
    number = length(bouton.nodes)
  ))
  ## Calculate on a single neuron basis
  dist.both = intersect(names(dist.in), names(dist.out))
  dists = c(dist.in[dist.both] + dist.out[dist.both])/2
  missing1 = names(dist.out)[!names(dist.out)%in%dist.both]
  missing2 = names(dist.in)[!names(dist.in)%in%dist.both]
  dists[missing1] = dist.out[missing1]
  dists[missing2] = dist.in[missing2]
  dists = dists / stepsize
  match.to.bouton = match(names(dists), boutons$treenode)
  boutons[match.to.bouton,]$interbouton.distance = dists
}
View(interbouton)

## Plot jitter plot!
pdf("images/boutons_euclidean_distance_to_midline.pdf", height = 10, width = 5)
ggplot(boutons, aes(x=layer, y=interbouton.distance, color = layer)) + 
  geom_jitter(position=position_jitter(0.2)) + 
  stat_summary(fun.data="mean_sdl", mult=1, 
               geom="crossbar", width=0.5, col = "grey70")+ 
  stat_summary(fun.data=mean_sdl, mult=1, 
               geom="pointrange", color="red")+
  scale_color_manual(values = layer.cols)+
  theme_minimal()
dev.off()

## NOTE: In order to learn more about plotting information from data frames with ggplot2 have a look here:
### ggplot2: http://www.sthda.com/english/wiki/ggplot2-essentials
### ggpubr: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/ # support for popular statistical comparisons !



