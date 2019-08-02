# set the scene
## This script assumed that you have run the file "11-frechter2019/00-setup.R"

# Now we are going to try to break the LH down into different voxels, based in its inputs and its outputs. 
# To do this I am going to create 25 supervoxels (I have tried a higher number, and the end result is similar) 
# based on predicted connectiivty groups. This involves breaking the rows (i.e. LH inputs) of the possible synapse matrix
# into 25 clusters, and using the nodes of the neurons in each cluster to generate a kernel density estimate in 3D space. 
# Each point in each input and LHN dendrite will then be given a score that represents it inclusion within each
# of the 25 3D kernel density plots.
# 
# (The result might be quite different if we use the input axons to define the voxels, rather than the outputs' dendrites,
# compared with the opposite scenario. We will look at both options, and compare with an NBlast clustering of severed dendrites
# and axons within the LH volume as well, for good measure)
#   
# We will also break down the output zones of the LH into 25 different voxels. This time, the LHON 
# axons are morphologically clustered via NBlast and the kernel density estimates are built on these 25 point clouds.

################################
### Get neuron data together ###
################################

# Get LH input neurons
inputs = lhns::most.lhins
inputs = subset(inputs,!anatomy.group%in%c("Centrifugal","notLHproper","MBON-a'2-a'3")) # Remove non-PNs
upns = subset(inputs,type=="uPN" & !is.na(glomerulus))
pns = c(pn.axons.light, upns)

# Choose which are LHONs
lhns.chosen = subset(lhns::most.lhns,cell.type!="notLHproper")
lhons = subset(lhns.chosen,type=="ON")
lhlns = subset(lhns.chosen,type=="LN")
lhn.dendrites = dendritic_cable(lhns.chosen)
lhn.axons = axonic_cable(lhns.chosen)

# Get LHN dendrites, near these PN axons
lhn.dendrites.lh = prune_in_volume(lhn.dendrites, neuropil = "LH_R",OmitFailures=TRUE)

# Neurons belonging to large cell types and dye fills
big.cts = names(table(most.lhns[,"cell.type"]))[table(most.lhns[,"cell.type"])>2]
bigs = subset(most.lhns,cell.type%in%big.cts&coreLH=="TRUE")
bigs = c(bigs,dye.fills[!names(dye.fills)%in%names(bigs)])
bigs[is.na(bigs[,"type"]),"type"] = "LN"

# Get LH volumes
LHR = ashape3d(xyzmatrix(subset(FCWBNP.surf, "LH_R")),alpha=30)
LHL = ashape3d(xyzmatrix(subset(FCWBNP.surf, "LH_L")),alpha=30)

# Get PN termini inside of this shape
points = xyzmatrix(pns)[inashape3d(points=xyzmatrix(pns),as3d=LHR),]
pns.termini = nlapply(pns, nat::prune, target = points, keep = 'near', maxdist = 0,OmitFailures=T)
lhr = as.mesh3d(LHR)


#############################
### Generate Super Voxels ###
#############################


# Number of super voxels that we want to create
clust.no = 25

# NBLAST the input axons
neuropil = LHR
points = xyzmatrix(inputs)
points = points[inashape3d(points=points,as3d=LHR),]
inputs.termini = nlapply(inputs, nat::prune, target = points, keep = 'near', maxdist = 0,OmitFailures=T)
inputs.termini = nlapply(inputs.termini, resample, stepsize = 1, OmitFailures = TRUE)
inputs.termini.dps = dotprops(inputs.termini, OmitFailures = TRUE)
lh.axon.morph = nblast_allbyall(inputs.termini.dps, normalisation = "mean", distance = F)
clust.nat = nhclust(scoremat = lh.axon.morph)
fsc=dendroextras::color_clusters(clust.nat, k=clust.no)
grp.input.nat = cutree(clust.nat,k=clust.no)
plot(fsc,leaflab="none",main="NBlast Clustering of Input Axons")

# Break axons into anatomical clusters, but remove high Strahler order branches
inputs.strahler = assign_strahler(inputs[sapply(inputs, function(x) x$nTrees==1)], OmitFailures = TRUE)
inputs.fragments = neuronlist()
for(ii in 1:length(inputs.strahler)){
  id = names(inputs.strahler[ii])
  df = inputs.strahler[ii][,]
  i = inputs.strahler[id][[1]]
  x = inputs.strahler[id][[1]]
  relevant.points = subset(x$d, PointNo%in%i$d$PointNo)
  i$d$strahler_order = relevant.points[match(i$d$PointNo,relevant.points$PointNo),]$strahler_order
  s = max(i$d$strahler_order)
  if(s>1){
    v = rownames(subset(i$d,strahler_order%in%c(s)))
    i = nat::prune_vertices(i, verticestoprune = v, invert = FALSE)
    relevant.points = subset(x$d, PointNo%in%i$d$PointNo)
    i$d$strahler_order = relevant.points[match(i$d$PointNo,relevant.points$PointNo),]$strahler_order
    for(t in 1:i$nTrees){
      if(length(unlist(i$SubTrees[t]))>1){
        subt = subtree(i, subtree = t)
        subt = as.neuronlist(subt)
        attr(subt,"df") = df
        names(subt) = paste0(id,"_",t)
        inputs.fragments = c(inputs.fragments,subt)  
      }
    }
  }else{
    i = as.neuronlist(i)
    attr(i,"df") = df
    names(i) = paste0(id,"_",t)
    inputs.fragments = c(inputs.fragments,as.neuronlist(i))
  }
}
inputs.fragments = inputs.fragments[summary(inputs.fragments)$cable.length>12] # Get rid of the short ones
inputs.fragments = nlapply(inputs.fragments, nat::prune, target = points, keep = 'near', maxdist = 0,OmitFailures=T)
inputs.fragments.dps = dotprops(inputs.fragments, OmitFailures = TRUE)

# NBLAST these axon fragments
inputs.fragments.morph = nblast_allbyall(inputs.fragments.dps, normalisation = "mean", distance = F)
clust.nat.frag = nhclust(scoremat = inputs.fragments.morph, maxneurons = 5000)
fsc2=dendroextras::color_clusters(clust.nat.frag, k=clust.no)
grp.input.nat.frag = cutree(clust.nat.frag,k=clust.no)
plot(fsc2,leaflab="none",main="NBlast Clustering of Fragmented Input Axons")

# Voxelise brain by PN domains
ims1 = list()
for(n in sort(unique(grp.input.nat.frag))){
  l = n
  points = xyzmatrix(inputs.fragments[names(grp.input.nat.frag)[grp.input.nat.frag==n]])
  H.pi <- Hpi(points,binned=TRUE)
  fhat <- kde(points, H=H.pi, xmin = boundingbox(FCWBNP.surf)[1,], xmax = boundingbox(FCWBNP.surf)[2,], w = rep(1/nrow(points),nrow(points))) # Area under curve = 1
  i = im3d(fhat$estimate,BoundingBox=boundingbox(FCWBNP.surf))
  attr(i,"contours") = fhat$cont
  attr(i,"fhat") = fhat
  ims1[[n]] <- i
}
names(ims1) = 1:length(ims1)
lhn.dendrites.marked.1 = nlapply(lhn.dendrites,assign_supervoxel.neuron,ims=ims1)

# Now we can look at output regions, elsewhere in the brain
lhon.axons.dps = dotprops(lhn.axons, OmitFailures = TRUE)
lhon.axon.morph = nblast_allbyall(lhon.axons.dps, normalisation = "mean", distance = F)
clust.nat.2 = nhclust(scoremat = lhon.axon.morph, method = "ward.D2")
fsc=dendroextras::color_clusters(clust.nat.2, k=clust.no)
plot(fsc,leaflab = "none")
grp.nat <- cutree(clust.nat.2, k = clust.no)

# Calculate voxel scores
ims2 = list()
for(n in sort(unique(grp.nat))){
  l = n
  points = catnat:::axonic.points.neuronlist(lhn.axons[names(lhn.axons)%in%names(grp.nat)[grp.nat==n]])
  # call the plug-in bandwidth estimator
  H.pi <- Hpi(points,binned=TRUE) ## b is a matrix of x,y,z points
  fhat <- kde(points, H=H.pi, xmax = boundingbox(FCWB)[2,],xmin = boundingbox(FCWB)[1,])
  i = im3d(fhat$estimate,BoundingBox=boundingbox(FCWB))
  attr(i,"contours") = fhat$cont
  attr(i,"fhat") = fhat
  ims2[[n]] <- i
}
names(ims2) = 1:length(ims2)
lhn.axons.marked.2 = nlapply(lhn.axons,assign_supervoxel.neuron,ims=ims2)

#######################################
### Calculate Voxel Inclusion Score ###
#######################################

# Big cell types
bigs.marked = subset(lhn.dendrites.marked.1, cell.type%in%bigs[,"cell.type"])
bigs.lns.marked = nlapply(subset(most.lhns, cell.type%in%subset(bigs,type=="LN")[,"cell.type"]),assign_supervoxel.neuron,ims=ims1)
bigs.marked = c(bigs.marked,bigs.lns.marked)

# uPNs
pns.termini = subset(pns.termini,!Glomerulus%in%c("VL2pP"))
pns.marked = nlapply(pns.termini,assign_supervoxel.neuron,ims=ims1)
p.df = attr(pns.marked,"df")
p.df[is.na(p.df$Glomerulus),]$Glomerulus = "V"
attr(pns.marked,"df") = p.df

# Calculate for PN axons
pns.overlap.matrix = matrix(0,ncol=length(unique(pns.marked[,"Glomerulus"])),nrow=clust.no)
colnames(pns.overlap.matrix) = unique(pns.marked[,"Glomerulus"])
cts = unique(pns.marked[,"Glomerulus"])
cts = cts[!is.na(cts)]
for (ct in cts){ # Average cell types
  pns.marked.ct = subset(pns.marked,Glomerulus==ct)
  pn.scores = rep(0,clust.no)
  for(id in names(pns.marked.ct)){
    neuron = pns.marked[id][[1]]
    # Calculate scores
    pn.scores = pn.scores + apply(neuron$d[,c(8:(7+clust.no))],2,sum)
  }
  pn.scores = pn.scores/length(pns.marked.ct)
  pns.overlap.matrix[,ct] = pn.scores
}
rownames(pns.overlap.matrix) = paste0(1:nrow(pns.overlap.matrix),"i")# Columns are letters for the LH regions
colnames(pns.overlap.matrix) = cts
pns.overlap.matrix.norm = apply(pns.overlap.matrix,2,function(x)x/max(x)) # Normalise axonic output
pns.rowclus = color_clusters(hclust(dist(pns.overlap.matrix,method = "euclidean"), method = "ward.D2" ),k=7)
pns.rowclus = dendextend::set(as.dendrogram(pns.rowclus),"branches_lwd",3)
pns.colclus = hclust(dist(t(pns.overlap.matrix),method = "euclidean"), method = "ward.D2" )
pns.colclus = dendextend::set(as.dendrogram(pns.colclus),"branches_lwd",3)

# What about other sorts of PNs? Like those strange multi-modal ones?
inputs.termini.marked = nlapply(inputs.termini,assign_supervoxel.neuron,ims=ims1)
inputs.termini.marked =  subset(inputs.termini.marked,!is.na(anatomy.group))
inputs.termini.overlap.matrix = matrix(0,ncol=length(unique(inputs.termini.marked[,"anatomy.group"])),nrow=clust.no)
colnames(inputs.termini.overlap.matrix) = unique(inputs.termini.marked[,"anatomy.group"])
cts = unique(inputs.termini.marked[,"anatomy.group"])
for (ct in cts){ # Average cell types
  inputs.termini.marked.ct = subset(inputs.termini.marked,anatomy.group==ct)
  pn.scores = rep(0,clust.no)
  for(id in names(inputs.termini.marked.ct)){
    neuron =inputs.termini.marked[id][[1]]
    # Calculate scores
    pn.scores = pn.scores + apply(neuron$d[,c(8:(7+clust.no))],2,sum)
  }
  pn.scores = pn.scores#/length(inputs.termini.marked.ct)
  inputs.termini.overlap.matrix[,ct] = pn.scores
}
rownames(inputs.termini.overlap.matrix) = paste0(1:nrow(inputs.termini.overlap.matrix),"i")# Columns are letters for the LH regions
colnames(inputs.termini.overlap.matrix) = cts
# We can also do this by modality
modalities =  as.character(sapply(colnames(inputs.termini.overlap.matrix),function(x) paste0(subset(inputs,anatomy.group==x)[,"modality"][1],"-",subset(inputs,anatomy.group==x)[,"neurotransmitter"][1])))
names(modalities) = colnames(inputs.termini.overlap.matrix)
modalities[colnames(inputs.termini.overlap.matrix)%in%c("AL-mALT-PN1","AL-lALT-PN1")] = "Uniglomerular.Olfactory-ACh"
modalities[grepl("^Olfactory-ACh$",modalities)] = "Multiglomerular.Olfactory-ACh"
modalities[grepl("^Unknown",modalities)] = "Expansive"
modalities[grepl("mlALT",names(modalities))] = "Inhibitory.Olfactory-GABA"
modalities = gsub("\\-.*","",modalities)
modalities = modalities[modalities!="NA"]
colnames(inputs.termini.overlap.matrix) = modalities
mod.inputs.termini.overlap.matrix = t(apply(t(inputs.termini.overlap.matrix), 2, function(x) tapply(x, colnames(inputs.termini.overlap.matrix), sum)))
mod.inputs.termini.overlap.matrix =mod.inputs.termini.overlap.matrix[,!grepl("Expansive|Centro|Centri|NA",colnames(mod.inputs.termini.overlap.matrix))]
mod.inputs.termini.overlap.matrix.norm = apply(mod.inputs.termini.overlap.matrix,2,function(x)x/max(x)) # Normalise axonic output
inputs.termini.rowclus=hclust(dist(mod.inputs.termini.overlap.matrix,method = "euclidean"), method = "ward.D2" )
inputs.termini.colclus =hclust(dist(t(mod.inputs.termini.overlap.matrix),method = "euclidean"), method = "ward.D2" )
inputs.termini.colclus = dendextend::set(as.dendrogram(inputs.termini.colclus),"branches_lwd",3)

# Calculate for LH dendrites
bigs.overlap.matrix = matrix(0,ncol=length(unique(bigs.marked[,"cell.type"])),nrow=clust.no)
colnames(bigs.overlap.matrix) = unique(bigs.marked[,"cell.type"])
cts = unique(bigs.marked[,"cell.type"])
for (ct in cts){ # Average cell types
  bigs.marked.ct = subset(bigs.marked,cell.type==ct)
  bigs.scores = rep(0,clust.no)
  for(id in names(bigs.marked.ct)){
    neuron =bigs.marked[id][[1]]
    # Calculate scores
    bigs.scores = bigs.scores + apply(neuron$d[,c(8:(7+clust.no))],2,sum)
  }
  bigs.scores = bigs.scores/length(bigs.marked.ct)
  bigs.overlap.matrix[,ct] = bigs.scores
}
colnames(bigs.overlap.matrix) = cts
rownames(bigs.overlap.matrix) = paste0(1:nrow(bigs.overlap.matrix),"i")# Columns are letters for the LH regions
# Different normalisations
bigs.overlap.matrix.scaled = scale(bigs.overlap.matrix) # Normalise axonic output
bigs.overlap.matrix.scaled.2 = t(scale(t(bigs.overlap.matrix)))  # Normalise dendritic input
bigs.overlap.matrix.scaled[is.na(bigs.overlap.matrix.scaled)] = 0
bigs.overlap.matrix.norm = apply(bigs.overlap.matrix,2,function(x)x/max(x)) # Normalise axonic output
bigs.overlap.matrix.norm.2 = t(apply(bigs.overlap.matrix,1,function(x)x/max(x))) # Normalise dendritic input
bigs.overlap.matrix.norm[is.na(bigs.overlap.matrix.norm)] = 0
# Colors
bigs.cols = ifelse(colnames(bigs.overlap.matrix.norm)%in%subset(bigs,type=="LN")[,"cell.type"],"green","blue")
# Cluster columns
bigs.colclus =hclust(dist(t(bigs.overlap.matrix.norm),method = "euclidean"), method = "ward.D2" )
bigs.colclus = dendextend::set(as.dendrogram(bigs.colclus),"branches_lwd",3)
overlap.matrix = matrix(0,ncol=clust.no,nrow=clust.no)
LHpointssummed = colSums(do.call(rbind,lapply(lhn.dendrites.marked.1,function(x) colSums(x$d[,c(8:(7+clust.no))]))))
LHnormaliser = matrix(LHpointssummed,ncol=clust.no,nrow=clust.no, byrow = TRUE)

# Cell types
cts = unique(lhn.axons.marked.2[,"cell.type"])[!is.na(unique(lhn.axons.marked.2[,"cell.type"]))]
for (ct in cts){ # Average cell types
  lhons.ct = subset(lhn.axons.marked.2,cell.type==ct)
  overlap.matrix.ct = matrix(0,ncol=clust.no,nrow=clust.no)
  for(id in names(lhons.ct)){
    neuron =lhn.axons.marked.2[id][[1]]
    neuron2 = lhn.dendrites.marked.1[id][[1]]
    # Calculate scores
    matrix1 = matrix(apply(neuron$d[,c(8:(7+clust.no))],2,sum),ncol=clust.no,nrow=clust.no,byrow = T)
    #matrix1 = matrix1 / LHnormaliser
    matrix2 = matrix(apply(neuron2$d[,c(8:(7+clust.no))],2,sum),ncol=clust.no,nrow=clust.no,byrow = F)
    overlap.matrix.neuron = matrix1 * matrix2
    overlap.matrix.ct = overlap.matrix.ct + overlap.matrix.neuron
  }
  overlap.matrix = overlap.matrix + (overlap.matrix.ct/length(lhons.ct))
}
colnames(overlap.matrix) = paste0(1:ncol(overlap.matrix),"o")# Columns are letters for the output regions
rownames(overlap.matrix) = paste0(1:nrow(overlap.matrix),"i")# Columns are letters for the LH regions
grp.overlap.input = cutree(pns.rowclus, k = 7, order_clusters_as_data = FALSE)
rownames(overlap.matrix) = paste0(grp.overlap.input[rownames(overlap.matrix)],"i")
overlap.matrix = apply(overlap.matrix, 2, function(x) tapply(x, rownames(overlap.matrix), sum))
overlap.matrix.norm = apply(overlap.matrix,2,function(x)x/max(x)) # Normalise axonic output
overlap.matrix.norm[is.na(overlap.matrix.norm)] = 0
colclus = color_clusters(hclust(dist(t(overlap.matrix),method = "euclidean"), method = "ward.D2" ),k=8)
colclus = dendextend::set(as.dendrogram(colclus),"branches_lwd",3)

# set colours
valences = unlist(sapply(colnames(pns.overlap.matrix),function(x) subset(lhns::pn.info,Glomerulus==x)[,"Valence(simplified)"][1] ))
valences[c("VM7","VA1lm")] = "Food"
valences[c("VA1lm")] = "Pheromone"
valences[c("V")] = "Aversive"
valences[grepl("VP",names(valences))] = "Humidity"
col.valences = c()
for(v in valences){
  if(grepl("Aversive",v)){
    i = "darkred"
  }else if(grepl("Food|Attractive",v)){
    i = "darkgreen"
  }else if(grepl("Egg",v)){
    i = "darkorange"
  }else if(grepl("Humidity",v)){
    i = "blue"
  }else if(grepl("Phero",v)){
    i = "purple"
  }else{
    i = "darkgrey"
  }
  col.valences = c(col.valences,i)
}
colnames(inputs.termini.overlap.matrix) = colnames(mod.inputs.termini.overlap.matrix.norm) = c("Gust", "iPN", "Mech", "mPN", 
                                                                                        "Mod", "Olf+Gust", "Thermo", "uPN")
mod.cols = c("green","cyan","blue","orange","magenta","brown","red","yellow")


# Plot a final heatmap!
Heatmap(mod.inputs.termini.overlap.matrix.norm, name = "",
        cluster_rows = pns.rowclus, cluster_columns = inputs.termini.colclus, 
        column_title ="",row_title ="",
        row_names_gp = gpar(fontsize = 0, fontface = "plain",col = "chartreuse2"),column_names_gp = gpar(fontsize = 20, fontface = "plain",col=mod.cols))+
  Heatmap(pns.overlap.matrix.norm, name = "",
          cluster_rows = pns.rowclus, cluster_columns = pns.colclus, 
          column_title ="",row_title ="",
          row_names_gp = gpar(fontsize = 0, fontface = "plain",col = "chartreuse2"),column_names_gp = gpar(fontsize = 20, fontface = "plain",col=col.valences)) +
  Heatmap(bigs.overlap.matrix.norm, name = "",
          cluster_rows = pns.rowclus, cluster_columns = bigs.colclus, 
          column_title ="",row_title ="",
          row_names_gp = gpar(fontsize = 20, fontface = "plain",col = "chartreuse2"),column_names_gp = gpar(fontsize = 20, fontface = "plain",col=bigs.cols))

##############################
### Visualise Syper Voxels ###
##############################

# Visualise input voxels
open3d(userMatrix = structure(c(0.955159366130829, -0.285889625549316, 
                                 -0.0770552754402161, 0, -0.296017527580261, -0.927844762802124, 
                                 -0.226886421442032, 0, -0.00663087517023087, 0.23952242732048, 
                                 -0.970868229866028, 0, -0.810615502960502, -1.20646136620167, 
                                 -7.49578407532869e-06, 1), .Dim = c(4L, 4L)), zoom = 0.676839649677277, 
        windowRect = c(1460L, 65L, 2585L, 1058L))
for(n in sort(unique(grp.overlap.input))){
  clusters = as.numeric(gsub("i","",names(grp.overlap.input)[grp.overlap.input==n]))
  plot3d(lhr, col = "grey", alpha = 0.1, add = TRUE)
  for(c in clusters){
    fhat = attr(ims1[[c]], "fhat")
    col = rainbow(clust.no)[c]
    plot(fhat,cont=c(50,95),colors=c(col,col),drawpoints=F
         ,xlab="", ylab="", zlab="",size=2, ptcol="white", add=TRUE, box=FALSE, axes=FALSE)
  }
  rgl.snapshot(paste0("images/voxelcluster_","input_",n,".png"),fmt="png")
  clear3d()
}
open3d(userMatrix=structure(c(0.955159366130829, -0.285889625549316, -0.0770552754402161, 
                              0, -0.296017527580261, -0.927844762802124, -0.226886421442032, 
                              0, -0.00663087517023087, 0.23952242732048, -0.970868229866028, 
                              0, 0, 0, 0, 1), .Dim = c(4L, 4L)),zoom=0.746215641498566, windowRect = c(1440L, 45L, 2867L, 1040L))
text.location = boundingbox(lhr)[1,]+c(80,0,0)
for(c in 1:length(ims1)){
  plot3d(lhr, col = "grey", alpha = 0.1, add = TRUE)
  fhat = attr(ims1[[c]], "fhat")
  col = "green"
  plot(fhat,cont=c(50,95),colors=c(col,col),drawpoints=F
       ,xlab="", ylab="", zlab="",size=2, ptcol="white", add=TRUE, box=FALSE, axes=FALSE)
  text3d(text.location,texts = paste0(c,"i"),col=col,font=2, cex = 1)
  rgl.snapshot(paste0("images/supervoxels_","input_",c,".png"),fmt="png")
  clear3d()
}
# and a dorsal view
open3d(userMatrix = structure(c(0.954851269721985, 0.040564451366663, 
                                0.294301688671112, 0, 0.288105726242065, 0.115264102816582, -0.950636386871338, 
                                0, -0.0724844262003899, 0.992506504058838, 0.098373107612133, 
                                0, 1.55491255847045, -1.43691751944516, -1.79519264520422e-06, 
                                1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, windowRect = c(1500L, 
                                                                                                105L, 2846L, 1030L))
for(c in 1:length(ims1)){
  plot3d(lhr, col = "grey", alpha = 0.1, add = TRUE)
  fhat = attr(ims1[[c]], "fhat")
  col = "green"
  plot(fhat,cont=c(50,95),colors=c(col,col),drawpoints=F
       ,xlab="", ylab="", zlab="",size=2, ptcol="white", add=TRUE, box=FALSE, axes=FALSE)
  rgl.snapshot(paste0("images/supervoxels_dorsalview_","input_",c,".png"),fmt="png")
  clear3d()
}


# Visualise output voxels
grp.overlap.output = cutree(colclus, k = 8,order_clusters_as_data = FALSE)
open3d(userMatrix = structure(c(0.998224198818207, 0.0104750925675035, 
                                -0.0586415156722069, 0, 0.00851582363247871, -0.999400317668915, 
                                -0.0335615389049053, 0, -0.0589578822255135, 0.0330026298761368, 
                                -0.997714698314667, 0, 2.96389417437521, -12.2722030161793, 2.34399720966394e-06, 
                                1), .Dim = c(4L, 4L)), zoom = 0.436296999454498, windowRect = c(1526L, 
                                                                                                44L, 2811L, 1058L))
for(n in sort(unique(grp.overlap.output))){
  clusters = as.numeric(gsub("o","",names(grp.overlap.output)[grp.overlap.output==n]))
  plot3d(FCWB, col = "grey", alpha = 0.1, add = TRUE)
  for(c in clusters){
    fhat = attr(ims2[[c]], "fhat")
    col = rainbow(clust.no)[c]
    plot(fhat,cont=c(50,95),colors=c(col,col),drawpoints=F
         ,xlab="", ylab="", zlab="",size=2, ptcol="white", add=TRUE, box=FALSE, axes=FALSE)
  }
  rgl.snapshot(paste0("images/voxelcluster_","output_",n,".png"),fmt="png")
  clear3d()
}
open3d(userMatrix = structure(c(0.998224198818207, 0.0104750925675035, 
                                -0.0586415156722069, 0, 0.00851582363247871, -0.999400317668915, 
                                -0.0335615389049053, 0, -0.0589578822255135, 0.0330026298761368, 
                                -0.997714698314667, 0, 2.96389417437521, -12.2722030161793, 2.34399720966394e-06, 
                                1), .Dim = c(4L, 4L)), zoom = 0.436296999454498, windowRect = c(1526L, 
                                                                                                44L, 2811L, 1058L))
plot3d(FCWB, col = "grey", alpha = 0.1, add = TRUE)
text.location = boundingbox(FCWB)[1,]+c(100,0,0)
for(c in 1:length(ims2)){
  plot3d(FCWB, col = "grey", alpha = 0.1, add = TRUE)
  fhat = attr(ims2[[c]], "fhat")
  col = "magenta"
  plot(fhat,cont=c(50,95),colors=c(col,col),drawpoints=F
       ,xlab="", ylab="", zlab="",size=2, ptcol="white", add=TRUE, box=FALSE, axes=FALSE)
  rgl.snapshot(paste0("/images/supervoxels_","output_",c,".png"),fmt="png")
  clear3d()
}










