### This script assumes you have run "/00-basics/00-setup.R"

# Let's see if we can cut through the main tract of some neurons, and then see how they are poistioned in that tract
## For extyra fun, we can compate a dataset derived from light microscopy, and one from electron micrographs.
## We can then see the benefits from a lot of data from a single brain, rather than sparse registered data from many.

# Connect to the public FAFB instance (Zheng et al. 2018) hosted publicly by Virtual Fly Brain
adult.conn = catmaid_login(server="https://catmaid-fafb.virtualflybrain.org/")

# Get all EM D. melanogaster uniglomerular olfactory projeciton neurons
em.pns = read.neurons.catmaid("annotation:^PN glomerulus")

# Ah! But they are in a non-standard EM space?
## We need to put them into a light-level template space
em.pns.fcwb = xform_brain(em.pns, sample = FAFB14, reference = FCWB)
### Phew! That's a big bottleneck we just jumped around.

# Let's update the meta-data, to record that they are EM reconstructions
em.pns.fcwb[,"data.type"] = "EM"

# And we can get our light-level data from the lhns package.
## These are mostly FlyCircuit neurons
AL_mALT_PN1 = subset(lhns::most.lhins, anatomy.group == "AL-mALT-PN1")
AL_mALT_PN1[,"data.type"] = "light"
### These are already in our template space, FCWB

# All together now
all.pns = c(AL_mALT_PN1,em.pns.fcwb)

# Plot axon tract locations on a plane, light versus EM neurons, to reveal EM organisation
## Let's do a bit of cutting to get a better look in at this tract
spines = nlapply(all.pns, simplify_neuron, n = 0)
cropped = nlapply(spines, subset, Y<150 & Y > 75 & Z>55 & Z<75, OmitFailures = TRUE)
tract_vector <- function(n, xval=250, thresh=100) {
  p=xyzmatrix(n)
  near=abs(p[,'X']-xval)<thresh
  pc=prcomp(p[near,])
  c(colMeans(p[near,]), pc$rotation[,'PC1'])
}
tvs=sapply(cropped, tract_vector)
mean_cent=rowMeans(tvs[1:3,])
mean_vec=rowMeans(tvs[4:6,])
make_centred_plane <- function(point, normal, scale = 1) {
  # find the two orthogonal vectors to this one
  uv = Morpho::tangentPlane(normal)
  uv = sapply(uv, "*", scale, simplify = F)
  qmesh3d(
    cbind(
      point + uv$z,
      point + uv$y,
      point - uv$z,
      point - uv$y
    ),
    homogeneous = F,
    indices = 1:4
  )
}

# Plot
nopen3d(userMatrix = structure(c(0.958560526371002, 0.0707437470555305, 
                                 -0.275965690612793, 0, -0.173726290464401, -0.622576534748077, 
                                 -0.763031601905823, 0, -0.225789546966553, 0.779354691505432, 
                                 -0.584487438201904, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.550000011920929, 
        windowRect = c(20L, 65L, 1402L, 887L))
plane=make_centred_plane(mean_cent, mean_vec, scale=15)
plot3d(AL_mALT_PN1, col = "lightgrey", soma = TRUE)
plot3d(em.pns.fcwb, col='black', soma = TRUE)
shade3d(plane, col='red')
rgl.snapshot(filename = "images/basic_nat/nat_basic_mALT_tract_intersection.png", fmt = "png")

# Find and center intersections
plc=plane_coefficients(mean_cent, mean_vec)
intersections=t(sapply(cropped, intersect_plane, plane = plc, closestpoint = mean_cent))
meta = spines[sapply(intersections, length)>0,]
intersections = do.call(rbind, intersections)
intersections.cent=scale(intersections, center = T, scale = F)
d=sqrt(rowSums(intersections.cent^2))

# Finally we can construct a 2D coordinate system on the plane and project the intersection positions onto that
# find pair of orthogonal tangent vectors of plane
uv = Morpho::tangentPlane(mean_vec)

# centre our points on mean
sp=scale(intersections, scale = F, center=T)
DotProduct=function(a,b) sum(a*b)

# find coordinates w.r.t. to our two basis vectors
xy=data.frame(u=apply(sp, 1, DotProduct, uv[[1]]),
              v=apply(sp, 1, DotProduct, uv[[2]]))

# maintain consist x-y aspect ratio 
plot(xy, pch=19, asp=1, main = "Position of Axon intersection in plane", 
     xlab="u /µm", ylab="v /µm")

# now ggplot2 it
xy = cbind(meta, xy)
hull_cyl <- xy %>%
  group_by(data.type) %>%
  slice(chull(u, v))

pdf("images/basic_nat/nat_basic_PN_mALT_intersections.pdf", width = 10, height = 5)
ggplot(xy, aes(x= u, y= v, color = data.type, fill = data.type)) +
  geom_point() +
  facet_grid(data.type~.) +
  geom_polygon(data = hull_cyl,alpha = 0.1)+
  theme_minimal() + 
  theme(legend.position="none", text = element_text(size=20))
dev.off()

pdf("images/basic_nat/nat_basic_PN_mALT_intersections_legend.pdf", width = 10, height = 5)
ggplot(xy, aes(x= u, y= v, color = data.type, fill = data.type)) +
  geom_point() +
  facet_grid(data.type~.) +
  geom_polygon(data = hull_cyl,alpha = 0.1)+
  theme_minimal()
dev.off()
