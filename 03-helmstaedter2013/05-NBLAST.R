# set the scene
## This script assumed that you have run the file "03-helmstaedter2013/00-setup.R" and "03-helmstaedter2013/01-download.R"

# load the neuronlist class skeleton data that we have made
load('skn.rda')
load("skeleton_metadata.rda")
soma.points = t(t(skeleton_metadata$kn.e2006.ALLSKELETONS.FINAL2012.allSomata[,1:3])*c(16.5,16.5,25))

# Center all the neurons at 0,0,0 before clustering
zero_neuron <- function(neuron, align.pc = TRUE){
  center = colMeans(xyzmatrix(neuron))
  xyzmatrix(neuron) = xyzmatrix(neuron) - matrix(center,ncol=3, nrow = nrow(xyzmatrix(neuron)), byrow = TRUE)
  if(align.pc){
    xyzmatrix(neuron) = Morpho::pcAlign(xyzmatrix(neuron)) 
  }
  neuron
}
skn.soma = subset(skn, cellid.soma%in%1:nrow(soma.points))
skn.summary = summary(skn.soma)
keep = rownames(skn.summary)[skn.summary$cable.length>500000]
skn.soma = skn.soma[keep]
skn_zeroed = nlapply(skn.soma, zero_neuron, align.pc = FALSE, OmitFailures = TRUE)

## Let's have a look at the effect of zeroing a neuron
### Two similar neurons, far away
nopen3d(userMatrix = structure(c(-0.164959013462067, -0.981998920440674, 
                                 -0.0920134484767914, 0, 0.80710905790329, -0.0807799398899078, 
                                 -0.584850072860718, 0, 0.566889405250549, -0.170741215348244, 
                                 0.80590558052063, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 1, 
        windowRect = c(1460L, 65L, 3007L, 1043L))
examples = c("sk2551", "sk2953")
examples.soma = soma.points[skn[examples,"cellid.soma"],]
examples.neurons = skn[examples]
spheres3d(examples.soma, radius = 5000, col = c("red", "cyan"))
plot3d(examples.neurons, lwd = 3, col = c("red", "cyan"))
rgl.snapshot(filename = paste0("images/nat_mouse_retina_similar_neurons.png"), fmt = "png")
### Two similar neurons, put together
nopen3d()
for(i in 1:length(examples.neurons)){
  center = colMeans(xyzmatrix(examples.neurons[i]))
  xyzmatrix(examples.neurons[i]) = xyzmatrix(examples.neurons[i]) - matrix(center,ncol=3, nrow = nrow(xyzmatrix(examples.neurons[i])), byrow = TRUE)
  examples.soma[i,] = examples.soma[i,] - center
}
spheres3d(examples.soma, radius = 5000, col = c("red", "cyan"))
plot3d(examples.neurons, lwd = 3, col = c("red", "cyan"))
rgl.snapshot(filename = paste0("images/nat_mouse_retina_similar_neurons_zeroed.png"), fmt = "png")

# get surrounding mesh
points = unique(xyzmatrix(subset(skn, class == "glial")))
points = points[sample(1:nrow(points),10000),]
retina.chunk = alphashape3d::ashape3d(points, alpha = 10000, pert = TRUE)
retina.chunk = as.mesh3d(retina.chunk)

# show neurons
nopen3d(userMatrix = structure(c(-0.0486575998365879, 0.996040105819702, 
                                              -0.0744036510586739, 0, -0.0785038694739342, 0.0704479366540909, 
                                              0.994421482086182, 0, 0.995725631713867, 0.0542270317673683, 
                                              0.0747653767466545, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.75, 
        windowRect = c(20L, 65L, 1346L, 902L))
plot3d(skn, col = rainbow(length(skn))[sample(length(skn))])
spheres3d(soma.points, radius = 2000, col = rainbow(nrow(soma.points))[sample(nrow(soma.points))])
plot3d(retina.chunk, alpha = 0.1, col = "lightgrey", add = TRUE)
rgl.snapshot(filename ="images/mouse_inner_plexiform_neurons_uncoloured.png" ,fmt = "png")

# and now, coloured by cell type
clear3d()
plot3d(subset(skn, class == "amacrine"), col = "red")
spheres3d(soma.points[subset(skn, class == "amacrine")[,"cellid.soma"],],radius = 2000, col = "darkred")
plot3d(subset(skn, class == "bipolar"), col = "green")
spheres3d(soma.points[subset(skn, class == "bipolar")[,"cellid.soma"],],radius = 2000, col = "darkgreen")
plot3d(subset(skn, class == "ganglion"), col = "cyan")
spheres3d(soma.points[subset(skn, class == "ganglion")[,"cellid.soma"],],radius = 2000, col = "blue")
plot3d(subset(skn, class == "glial"), col = "magenta")
spheres3d(soma.points[subset(skn, class == "glial"&cellid.soma<1161)[,"cellid.soma"],],radius = 2000, col = "purple")
plot3d(retina.chunk, alpha = 0.1, col = "lightgrey", add = TRUE)
rgl.snapshot(filename ="images/mouse_inner_plexiform_neurons.png" ,fmt = "png")

# show the result of an in situ NBLAST, for neurons with a soma and over a certain cable length
k = 24
skn.dps = dotprops(skn.soma, stepsize = 1000, OmitFailures = TRUE)/1000
result.retina = nblast_allbyall(skn.dps, normalisation = "mean")
pdf("images/nat_mouse_retina_nblast_dendrogram.pdf", width = 5, height = 3)
hclust=nhclust(scoremat = result.retina)
dkcs=colour_clusters(hclust, k = k) %>% set("branches_lwd", 2)
plot(dkcs, leaflab = "none")
dev.off()

# Show neurons coloured by NBLAST
clear3d()
plot3d(hclust, db = skn, soma = TRUE, k = k)
for(i in 1:k){
  neuronids = subset(hclust, groups = i, k = k)
  somaids = skn[neuronids,"cellid.soma"]
  somaids = somaids[somaids<=1160]
  spheres3d(soma.points[somaids,],radius = 2000, col = cols[i])
}
plot3d(retina.chunk, alpha = 0.1, col = "lightgrey", add = TRUE)
rgl.snapshot(filename = "images/nat_mouse_retina_nblast_neurons.png", fmt = "png")

# Show the result of an NBLAST, now for zeroed neurons
skn.zeroed.dps = dotprops(skn_zeroed, stepsize = 1000, OmitFailures = TRUE)/1000
result_zeroed = nblast_allbyall(skn.zeroed.dps, normalisation = "mean")
pdf("images/nat_mouse_retina_nblast_zeroed_dendrogram.pdf", width = 5, height = 3)
hclust=nhclust(scoremat = result_zeroed)
dkcs=colour_clusters(hclust, k = k) %>% set("branches_lwd", 2)
plot(dkcs, leaflab = "none")
dev.off()

# Show neurons coloured by NBLAST
clear3d()
plot3d(hclust, db = skn, soma = TRUE, k = k)
cols = rainbow(k)
for(i in 1:k){
  neuronids = subset(hclust, groups = i, k = k)
  somaids = skn[neuronids,"cellid.soma"]
  somaids = somaids[somaids<=1160]
  spheres3d(soma.points[somaids,],radius = 2000, col = cols[i])
}
plot3d(retina.chunk, alpha = 0.1, col = "lightgrey", add = TRUE)
rgl.snapshot(filename = "images/nat_mouse_retina_zeroed_nblast_neurons.png", fmt = "png")

# Show neurons coloured by NBLAST
nopen3d(userMatrix = structure(c(-0.0486575998365879, 0.996040105819702, 
                                 -0.0744036510586739, 0, -0.0785038694739342, 0.0704479366540909, 
                                 0.994421482086182, 0, 0.995725631713867, 0.0542270317673683, 
                                 0.0747653767466545, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.75, 
        windowRect = c(20L, 65L, 1346L, 902L))
cols = rainbow(k)
for(i in 1:k){
  clear3d()
  plot3d(retina.chunk, alpha = 0.1, col = "lightgrey", add = TRUE)
  plot3d(hclust, db = skn, soma = TRUE, groups = i, k = k, col = cols[i])
  neuronids = subset(hclust, groups = i, k = k)
  somaids = skn[neuronids,"cellid.soma"]
  somaids = somaids[somaids<=1160]
  spheres3d(soma.points[somaids,],radius = 2000, col = cols[i])
  nam = names(sort(table(skn[neuronids,c("class")]), decreasing = TRUE))[1]
  rgl.snapshot(filename = paste0("images/nat_mouse_retina_nblast_zeroed_neurons_",i,"_", nam,".png"), fmt = "png")
}

