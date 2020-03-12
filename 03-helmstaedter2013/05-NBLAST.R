# set the scene
## This script assumed that you have run the file "03-helmstaedter2013/00-setup.R" and "03-helmstaedter2013/01-download.R"

# load the neuronlist class skeleton data that we have made
load('skn.rda')
skn = subset(skn, stypeid != 0)
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
nopen3d(userMatrix = structure(c(0.0361445173621178, -0.98055773973465, 
                                 -0.192872509360313, 0, 0.683494925498962, -0.116543605923653, 
                                 0.72059154510498, 0, -0.72905969619751, -0.157872825860977, 0.665993869304657, 
                                 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 1, windowRect = c(-71L, 
                                                                                             184L, 1369L, 1040L))
examples = c("sk2551", "sk2953")
examples.soma = soma.points[skn[examples,"cellid.soma"],]
examples.neurons = skn[examples]
spheres3d(examples.soma, radius = 5000, col = c("red", "cyan"))
plot3d(examples.neurons, lwd = 3, col = c("red", "cyan"))
rgl.snapshot(filename = paste0("images/nat_mouse_retina_similar_neurons.png"), fmt = "png")
### Two similar neurons, put together
nopen3d(userMatrix = structure(c(-0.133114844560623, -0.990953505039215, 
                                 0.0170866698026657, 0, 0.427550107240677, -0.0729692280292511, 
                                 -0.901041984558105, 0, 0.894137382507324, -0.112636663019657, 
                                 0.433395445346832, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 1, 
        windowRect = c(84L, 76L, 1436L, 913L))
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

# types
amacrine = subset(skn, class == "amacrine")
bipolar = subset(skn, class == "bipolar")
glia = subset(skn, class == "glial"&cellid.soma<1161)
ganglion = subset(skn, class == "ganglion")

# and now, coloured by cell type
clear3d()
plot3d(amacrine, col = "red")
spheres3d(soma.points[amacrine[,"cellid.soma"],],radius = 2000, col = "darkred")
plot3d(bipolar, col = "cyan")
spheres3d(soma.points[bipolar[,"cellid.soma"],],radius = 2000, col = "blue")
plot3d(ganglion, col = "green")
spheres3d(soma.points[ganglion[,"cellid.soma"],],radius = 2000, col = "darkgreen")
plot3d(glia, col = "magenta")
spheres3d(soma.points[glia[,"cellid.soma"],],radius = 2000, col = "purple")
plot3d(retina.chunk, alpha = 0.1, col = "lightgrey", add = TRUE)
rgl.snapshot(filename ="images/mouse_inner_plexiform_neurons.png" ,fmt = "png")

# and now, colour by subtype
clear3d()
plot3d(retina.chunk, alpha = 0.1, col = "lightgrey", add = TRUE)
cols = rainbow(length(unique(amacrine[,"stypeid"])))
names(cols) = unique(amacrine[,"stypeid"])
plot3d(amacrine, col = cols[as.character(amacrine[,"stypeid"])])
spheres3d(soma.points[amacrine[,"cellid.soma"],],radius = 2000, col = cols[as.character(amacrine[,"stypeid"])])
rgl.snapshot(filename ="images/mouse_inner_plexiform_amacrine_types.png" ,fmt = "png")

clear3d()
plot3d(retina.chunk, alpha = 0.1, col = "lightgrey", add = TRUE)
cols = rainbow(length(unique(bipolar[,"stypeid"])))
names(cols) = unique(bipolar[,"stypeid"])
plot3d(bipolar, col = cols[as.character(bipolar[,"stypeid"])])
spheres3d(soma.points[bipolar[,"cellid.soma"],],radius = 2000, col = cols[as.character(bipolar[,"stypeid"])])
rgl.snapshot(filename ="images/mouse_inner_plexiform_bipolar_types.png" ,fmt = "png")

clear3d()
plot3d(retina.chunk, alpha = 0.1, col = "lightgrey", add = TRUE)
cols = rainbow(length(unique(ganglion[,"stypeid"])))
names(cols) = unique(ganglion[,"stypeid"])
plot3d(ganglion, col = cols[as.character(ganglion[,"stypeid"])])
spheres3d(soma.points[ganglion[,"cellid.soma"],],radius = 2000, col = cols[as.character(ganglion[,"stypeid"])])
rgl.snapshot(filename ="images/mouse_inner_plexiform_ganglion_types.png" ,fmt = "png")

clear3d()
plot3d(retina.chunk, alpha = 0.1, col = "lightgrey", add = TRUE)
cols = rainbow(length(unique(glia[,"stypeid"])))
names(cols) = unique(glia[,"stypeid"])
plot3d(glia, col = cols[as.character(glia[,"stypeid"])])
spheres3d(soma.points[glia[,"cellid.soma"],],radius = 2000, col = cols[as.character(glia[,"stypeid"])])
rgl.snapshot(filename ="images/mouse_inner_plexiform_glia_types.png" ,fmt = "png")

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
cols = rainbow(k)
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
result_zeroed = nblast_allbyall(skn.zeroed.dps, normalisation = "raw")
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

# Show tSNE plot
nb <- result_zeroed
nb_scaled <- scale(nb)
nb = as.data.frame(nb)

## Curating the database for analysis with both t-SNE and PCA
labels <- skn_zeroed[,"stypeid"]

## for plotting
colors = rainbow(length(unique(labels)))
names(colors) = unique(labels)

# Skeleton.types
shapes = skn_zeroed[rownames(nb),"class"]
shapes = as.character(shapes)

## Executing the algorithm on our data
tsne <- Rtsne(nb, dims = 2, perplexity=30, verbose=TRUE, max_iter = 20000, 
              pca  = TRUE, initial_dims = 50, pca_center = TRUE, 
              pca_scale = TRUE, check_duplicates = FALSE)

# Ggplot
pdf("images/nat_retinal_NBLAST_tSNE.pdf", width = 10, height = 5)
tsne.df <- data.frame(tSNE1 = tsne$Y[,1], tSNE2 = tsne$Y[,2], col = colors[labels], cell.type = labels, shape = shapes)
tsne.df <- subset(tsne.df, shape != "")
hull_cyl <- tsne.df %>%
  group_by(cell.type) %>%
  slice(chull(tSNE1, tSNE2))
ggplot(tsne.df) + 
  geom_point(aes(x=tSNE1, y=tSNE2, color= as.character(cell.type), shape = as.character(shape))) + 
  scale_shape_manual(values = c(16, 17, 15, 1)) +
  theme(legend.position="none") + 
  geom_polygon(data = hull_cyl, aes(x = tSNE1, y = tSNE2, fill = as.character(cell.type)),alpha = 0.1)+
  theme_minimal() +
  guides(fill=FALSE, color = FALSE, shape = FALSE) + 
  ggtitle("")
dev.off()
# circles, amacrine
# triangles, bipolar
# square, ganglion
# open cirlces, glia

# So it seems like our NBLAST clusters have some relationship to, 
# but do not exactly fit the manual cell type classifications
# We can scan through the NBLAST clusters, and assess if they really belong
nopen3d(userMatrix = structure(c(-0.0486575998365879, 0.996040105819702, 
                                 -0.0744036510586739, 0, -0.0785038694739342, 0.0704479366540909, 
                                 0.994421482086182, 0, 0.995725631713867, 0.0542270317673683, 
                                 0.0747653767466545, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.75, 
        windowRect = c(20L, 65L, 1346L, 902L))
plot3d(retina.chunk, alpha = 0.1, col = "lightgrey", add = TRUE)
group.1 = subset(hclust, soma = TRUE, groups = 1, k = k)
plot3d(skn[group.1], col = "darkgrey")
somaids = skn[group.1,"cellid.soma"]
somaids = somaids[somaids<=1160]
spheres3d(soma.points[somaids,],radius = 2000, col = "darkgrey")
rgl.snapshot(filename = paste0("images/findneuron_group1_.png"), fmt = "png")
nlscan(skn[group.1], col = "red", lwd = 3)

# You can also find a specific neuron, using find.neuron, try it out
f = find.neuron(db = skn[group.1])
npop3d()
if (!exists("f")){
 f = "sk1324" 
}
plot3d(skn[f], col = "cyan", lwd = 4)
spheres3d(soma.points[skn[f,"cellid.soma"],],radius = 3000, col = "cyan")
rgl.snapshot(filename = paste0("images/findneuron_group1_found.png"), fmt = "png")

# The soma are detached in this scheme. We can also scan throguh with them though.
for(g1 in group.1){
  clear3d()
  plot3d(retina.chunk, alpha = 0.1, col = "lightgrey", add = TRUE)
  plot3d(skn[group.1], col = "darkgrey")
  plot3d(skn[g1], col = "red", lwd = 4)
  somaids = skn[group.1,"cellid.soma"]
  somaids = somaids[somaids<=1160]
  spheres3d(soma.points[somaids,],radius = 2000, col = "darkgrey")
  spheres3d(soma.points[skn[g1,"cellid.soma"],],radius = 3000, col = "red")
  nam = names(sort(table(skn[g1,c("class")]), decreasing = TRUE))[1]
  rgl.snapshot(filename = paste0("images/nlscan_group1_",nam,"_",g1,".png"), fmt = "png")
}








