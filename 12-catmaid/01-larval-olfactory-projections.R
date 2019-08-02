### This script assumes you have run "/12-catmaid/00-setup.R"

# Connect to the larval L1 CATMAID instance (Ohyama et al. 2015) hosted publicly by Virtual Fly Brain
larva.conn = catmaid_login(server="https://l1em.catmaid.virtualflybrain.org/")

# Let's fetch some uniglomerular olfactory projection neurons from the larva, reconstructed by Berck et al. 2016
upns = assignside(read.neurons.catmaid("annotation:^ORN PNs$", .progress='text'))/1000
upns = flow.centrality(upns, polypre = TRUE, mode = "centrifugal")
r.upns = c(subset(upns, grepl("right|Right|_r|R$", name))) # add  the bilateral neuron
l.upns = c(subset(upns, grepl("left|Left|_l|L$", name)))

# Also get some neuroanatomical landmarks and other data
load("data/almesh.rda")
load("data/lhmesh.rda")
load("data/mbmesh.rda")
load("data/camesh.rda")
l1mesh=read.hxsurf("data/l1mesh.microns.surf")
l1 = xyzmatrix(l1mesh)
points = cbind(l1mesh$Vertices[,-4],1)
cortexmesh = tmesh3d(t(points),t(l1mesh$Regions$Inside))

# Let's have a look at what one larval PN looks like
nopen3d(userMatrix = structure(c(0.999479055404663, 0.0313623659312725, 
                                 -0.00758622027933598, 0, 0.0314145050942898, -0.999482750892639, 
                                 0.00685763731598854, 0, -0.00736722396686673, -0.00709245353937149, 
                                 -0.999947726726532, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.440431654453278, 
        windowRect = c(1483L, 49L, 2784L, 1064L))
dm1l = subset(upns, name == "42b PN right" )
syn = get.synapses(dm1l, "POST")[,c(1:3)]
plot3d(cortexmesh, alpha = 0, add = TRUE)
plot3d(almesh, alpha = 0.1, add = TRUE, col = "red")
plot3d(lhmesh, alpha = 0.1, add = TRUE, col = "green")
plot3d(camesh, alpha = 0.1, add = TRUE, col = "purple")
plot3d(mbmesh, alpha = 0.1, add = TRUE, col = "magenta")
plot3d(dm1l, WithConnectors = TRUE, col = "black", lwd = 2, soma = TRUE)
rgl.snapshot(filename = "images/l1_glomerulus_42b_PN.png", fmt = "png")
clear3d()
plot3d(cortexmesh, alpha = 0, add = TRUE)
plot3d(almesh, alpha = 0.1, add = TRUE, col = "red")
plot3d(lhmesh, alpha = 0.1, add = TRUE, col = "green")
plot3d(camesh, alpha = 0.1, add = TRUE, col = "purple")
plot3d(mbmesh, alpha = 0.1, add = TRUE, col = "magenta")
points3d(syn, col = "cyan", lwd = 2)
rgl.snapshot(filename = "images/l1_glomerulus_42b_PN_post.png", fmt = "png")
clear3d()
plot3d(cortexmesh, alpha = 0, add = TRUE)
plot3d(almesh, alpha = 0.1, add = TRUE, col = "red")
syn = syn[pointsinside(syn, almesh),]
points3d(syn, col = "cyan", lwd = 2)
rgl.snapshot(filename = "images/l1_glomerulus_42b_PN_ALpost.png", fmt = "png")
a = as.mesh3d(ashape3d(x = syn, pert = TRUE, alpha = 3))
plot3d(a, col = "green", add = TRUE)
rgl.snapshot(filename = "images/l1_glomerulus_42b_PN_shape.png", fmt = "png")
clear3d()

# Calculate volumes
upns.resampled = resample(upns, 0.1)
volumes = post = pre = cable.length = c()
axon.post = axon.pre = axon.cable.length = c()
ca.post = ca.pre = lh.post = lh.pre = c()
glomeruli = list()
nopen3d(userMatrix = structure(c(0.967234313488007, -0.0323957242071629, 
                                 -0.251809895038605, 0, -0.0299824699759483, -0.999460399150848, 
                                 0.0134155247360468, 0, -0.252108603715897, -0.00542604364454746, 
                                 -0.967683672904968, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.45, 
        windowRect = c(20L, 65L, 1047L, 837L))
plot3d(cortexmesh, alpha = 0.05, col = "lightgrey", add = TRUE)
plot3d(almesh, alpha = 0.1, col = "darkgrey", add = TRUE)
cols = rainbow(length(upns))[sample(length(upns))]
for(i in 1:length(upns)){
  synapses.all = as.data.frame(get.synapses(upns[[i]], "BOTH"))
  nodes.all = xyzmatrix(upns.resampled[[i]])
  # AL
  synapses = synapses.all[pointsinside(synapses.all, almesh),]
  post = c(post, nrow(subset(synapses, prepost==1)))
  pre = c(pre, nrow(subset(synapses, prepost==0)))
  a = ashape3d(x = unique(as.matrix(synapses[,c(1:3)])), pert = TRUE, alpha = 5)
  v = volume_ashape3d(a)
  volumes = c(volumes, v)
  nodes = nrow(nodes.all[pointsinside(nodes.all, almesh),])
  cable.length = c(cable.length, nodes)
  glomeruli[[i]] = as.mesh3d(a)
  plot3d(glomeruli[[i]], col = cols[i], add = TRUE, alpha = 0.3)
  # Axon
  synapses.upper = subset(synapses.all, Y<60)
  axon.post = c(axon.post, nrow(subset(synapses.upper, prepost==1)))
  axon.pre = c(axon.pre, nrow(subset(synapses.upper, prepost==0)))
  nodes = subset(as.data.frame(nodes.all), Y<60)
  axon.cable.length = c(axon.cable.length, nrow(nodes))
  # Calyx
  ca.synapses = synapses.upper[pointsinside(synapses.upper, camesh),]
  ca.post = c(ca.post, nrow(subset(ca.synapses, prepost==1)))
  ca.pre = c(ca.pre, nrow(subset(ca.synapses, prepost==0)))
  # LH
  ln.synapses = synapses.upper[!pointsinside(synapses.upper, camesh),]
  lh.post = c(lh.post, nrow(subset(ln.synapses, prepost==1)))
  lh.pre = c(lh.pre, nrow(subset(ln.synapses, prepost==0)))
}
rgl.snapshot(filename = "images/l1_glomerular_volumes.png", fmt = "png")

# Put data together
df = upns[,]
df$volume = volumes
df$n.post = post
df$n.pre = pre
df$cable.length = cable.length/10
df$axon.n.post = axon.post
df$axon.n.pre = axon.pre
df$axon.cable.length = axon.cable.length/10
df$ca.n.post = ca.post
df$ca.n.pre = ca.pre
df$lh.n.post = lh.post
df$lh.n.pre = lh.pre
df$glom = gsub(" .*","",upns[,"name"])

# Plot
pdf("images/l1_PN_glomeruli_volumes.pdf", width = 10, height = 7)
ggplot(df, aes(y=volume, x=n.post)) + 
  geom_point(aes(size=n.post,color=glom))+
  geom_smooth(method=lm, col = "black")+geom_text(label= df$glom, size=3)+
  theme_minimal()
dev.off()

# Plot
pdf("images/l1_PN_input_vs_output.pdf", width = 10, height = 7)
ggplot(df, aes(y=axon.n.pre, x=n.pre)) + 
  geom_point(aes(size=axon.n.pre,color=glom))+
  geom_smooth(method=lm, col = "black")+geom_text(label= df$glom, size=3)+
  theme_minimal()
dev.off()

# Plot
pdf("images/l1_PN_glom_vs_outputcable.pdf", width = 10, height = 7)
ggplot(df, aes(y=axon.cable.length, x=volume)) + 
  geom_point(aes(size=axon.cable.length,color=glom))+
  geom_smooth(method=lm, col = "black")+geom_text(label= df$glom, size=3)+
  theme_minimal()
dev.off()

# Plot
pdf("images/l1_PN_glom_vs_outputcable.pdf", width = 10, height = 7)
ggplot(df, aes(y=ca.n.pre, x=volume)) + 
  geom_point(aes(size=ca.n.pre,color=glom))+
  geom_smooth(method=lm, col = "black")+geom_text(label= df$glom, size=3)+
  theme_minimal()
dev.off()

# Plot Ca and LH together
m = melt(df[,c("glom", "volume","ca.n.pre", "lh.n.pre")], id.vars =c("glom","volume"))
pdf("images/l1_PN_CA_LH_correlation_AL_volume.pdf", width = 5, height = 5)
ggplot(m, aes(y=value, x=volume, color = variable)) + 
  geom_point()+
  geom_smooth(method=lm) +
  scale_color_manual(values = c(ca.n.pre = "purple", lh.n.pre = "chartreuse3")) +
  theme_minimal()+
  theme(legend.position = "none")+
  stat_cor(aes(color = variable), method = "pearson", label.x = 3, label.y = c(70, 65))  # Add correlation coefficient
dev.off()
