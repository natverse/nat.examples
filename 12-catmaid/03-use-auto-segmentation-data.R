### This script assumes you have run "/12-catmaid/00-setup.R"

# The aim here is to use auto-segmentation to turn a manually traced, EM skeleton into a volumetric neuron!

# Connect to the public FAFB instance (Zheng et al. 2018) hosted publicly by Virtual Fly Brain
adult.conn = catmaid_login(server="https://catmaid-fafb.virtualflybrain.org/")

# PD2a1/b1 neurons - for MBON connectivity (Dolan et al. 2018)
em.pd2 = fetchn_fafb("name:^PD2a1/b1", mirror = F, reference = JFRC2)
em.pd2[,"type"] = "ON"

# Get its connectors in the LH
xyzmatrix(lhregion) <- xyzmatrix(Morpho::applyTransform(x=xyzmatrix(lhregion*1.25),trafo = getTrafo4x4(lhbigger)))
upstream.connectors = do.call(rbind,lapply(post,function(x) cbind(subset(x$connectors,prepost==1))))
i = pointsinside(upstream.connectors,lhregion,rval = "logical")
lh.upstream.connectors = upstream.connectors[i,]

# Search for potential connected segments
skids = "1299700"
pd2a1.mapping = map_fafbsegs_to_neuron(someneuronlist = pd2a1, node.match = 3, return.unmatched = FALSE, volume = volume)
pd2a1.mapping = subset(pd2a1.mapping,ngl_id!=0)

# Get meshes
## snapshots for all for PD2a1
ids = unique(pd2a1.mapping$ngl_id)
segs=find_merged_segments(ids)
nopen3d(userMatrix = structure(c(0.536078035831451, 0.423577815294266, 
                                 -0.730207145214081, 0, -0.355168461799622, -0.671545386314392, 
                                 -0.650294661521912, 0, -0.765817105770111, 0.60795521736145, 
                                 -0.20955890417099, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.310068160295486, 
        windowRect = c(20L, 65L, 1098L, 843L))
plot3d(FAFB14.surf, alpha = 0)
colours = rainbow(length(segs))
colours = sample(colours,length(segs))
for(i in 1:length(segs)){
  message(i)
  bmm = tryCatch(read_brainmaps_meshes(segs[i], volume = volume), error = function(e) NULL)
  plot3d(bmm, alpha = 1, add = TRUE, col = colours[i])
}
rgl.snapshot(filename = paste0( image.folder ,'images/visualise_segments_PD2a1_4.png'),fmt="png")
plot3d(lhns.done.flow["1299700"], lwd = 1, col = "darkgrey", soma = T)
rgl.snapshot(filename = paste0( image.folder ,'images/visualise_segments_PD2a1_PD2a1_4.png'),fmt="png")
