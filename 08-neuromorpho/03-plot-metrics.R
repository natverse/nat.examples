# set the scene
## This script assumed that you have run the file "08-neuromorpho/00-setup.R" and "01-download/00-setup.R"

## Get our neurons of interest
load("Jacobs_principal_neurons.rda")
load("Jacobs_interneurons.rda")

## plot differences in cable length between species
principals.summary = summary(principals)
principals.summary$species  = principals[, "species"]
principals.summary$bbx = sapply(principals[rownames(principals.summary)], function(x) prod(abs(boundingbox(x)[1,]-boundingbox(x)[2,])))
a = aggregate(principals.summary$branchpoints/principals.summary$cable.length, list(principals.summary$species), mean)
species.order = a[order(a$x, decreasing = TRUE),]$Group.1
principals.summary$species = factor(principals.summary$species, levels = species.order)

## Now we can explore things, lke the density of branchpoints
pdf("images/nat_neuromorpho_principals_branchpoints_over_cable.pdf", width = 10, height = 7)
principals.summary = subset(principals.summary, branchpoints/cable.length < 0.02)
ggplot(principals.summary, aes(x=species, y=branchpoints/cable.length, color=species)) +
  geom_jitter()+
  geom_boxplot(position=position_dodge(0.8), col = "darkgrey", alpha = 0)+
  theme_minimal() +
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1),
        legend.position="none")
dev.off()

## Or the volume of the bounding box around each neuron
pdf("images/nat_neuromorpho_principals_volume_over_cable.pdf", width = 10, height = 7)
principals.summary = subset(principals.summary, branchpoints/cable.length < 0.02)
ggplot(principals.summary, aes(x=species, y=bbx, color=species)) +
  geom_jitter()+
  geom_boxplot(position=position_dodge(0.8), col = "darkgrey", alpha = 0)+
  theme_minimal() +
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1),
        legend.position="none")
dev.off()

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

# Calculate sholl analysis
principals.sholl = nlapply(principals, sholl_analysis)
interneurons.sholl = nlapply(interneurons, sholl_analysis)
load("interneurons_sholl.rda")
load("principals_sholl.rda")
df = data.frame()
for(i in 1:length(principals.sholl)){
  ps = principals.sholl[[i]]
  ps = cbind(ps,principals.sholl[i,])
  df = rbind(df, ps)
}
a = aggregate(df$intersection, list(df$species, df$radii), mean)
colnames(a) = c("species","radii","intersections")

# Plot
a$species = factor(a$species, levels = species.order)
pdf("images/nat_neuromorpho_principals_sholl.pdf", width = 10, height = 7)
ggplot(a, aes(x=radii, y=intersections, color=species)) +
  geom_line()+
  theme_minimal()
dev.off()

# Plot interneurons
df = data.frame()
for(i in 1:length(interneurons.sholl)){
  ps = principals.sholl[[i]]
  ps = cbind(ps,interneurons.sholl[i,])
  df = rbind(df, ps)
}
a = aggregate(df$intersection, list(df$species, df$radii), mean)
colnames(a) = c("species","radii","intersections")
a$species = factor(a$species, levels = species.order)
pdf("images/nat_neuromorpho_interneurons_sholl.pdf", width = 10, height = 7)
ggplot(a, aes(x=radii, y=intersections, color=species)) +
  geom_line()+
  theme_minimal()
dev.off()

