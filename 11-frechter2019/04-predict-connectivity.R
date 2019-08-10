# set the scene
## This script assumed that you have run the file "11-frechter2019/00-setup.R"

# Now we are going to try to break the LH down into different voxels, based in its inputs and its outputs. 
# To do this I am going to create 25 supervoxels (I have tried a higher number, and the end result is similar) 
# based on predicted connectiivty groups. This involves breaking the rows (i.e. LH inputs) of the possible synapse matrix
# into 25 clusters, and using the nodes of the neurons in each cluster to generate a kernel density estimate in 3D space. 
# Each point in each input and LHN dendrite will then be given a score that represents it inclusion within each
# of the 25 3D kernel density plots.
#   
# We will also break down the output zones of the LH into 25 different voxels. This time, the LHON 
# axons are morphologically clustered via NBlast and the kernel density estimates are built on these 25 point clouds.

################################
### Get neuron data together ###
################################

# Get LH input neurons
inputs = lhns::most.lhins
inputs = subset(inputs,!anatomy.group%in%c("Centrifugal","notLHproper","MBON-a'2-a'3")) # Remove non-PNs
inputs = c(lhns::pn.axons.light, inputs)

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

# Let's just see what happens if we cluster these arbours as they are
### NBLAST the input axons
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
inputs.fragments = remove_highest_strahler(someneuronlist = inputs, neuropil = "LH_R") # 7174 fragments
inputs.fragments = inputs.fragments[summary(inputs.fragments)$cable.length>10] # Get rid of the short ones
inputs.fragments.dps = dotprops(inputs.fragments, OmitFailures = TRUE, resample = 0.1) # 2652 fragments

# NBLAST these axon fragments
inputs.fragments.morph = nblast_allbyall(inputs.fragments.dps, normalisation = "mean", distance = F)
clust.nat.frag = nhclust(scoremat = inputs.fragments.morph, maxneurons = 5000)
fsc2=dendroextras::color_clusters(clust.nat.frag, k=clust.no)
grp.input.nat.frag = cutree(clust.nat.frag,k=clust.no)
plot(fsc2,leaflab="none",main="NBlast Clustering of Fragmented Input Axons")

# Voxelise brain by PN domains
ims1 = calculate_voxels(inputs.fragments, grp.input.nat.frag, lhr)
lhn.dendrites.marked.1 = nlapply(lhn.dendrites,assign_supervoxel.neuron,ims=ims1)

# To make sense of what we just did, let's have a look at this DA2 PN as an example
nopen3d(userMatrix = structure(c(0.992447674274445, 0.0267316102981567, 
                                 0.119722224771976, 0, 0.0416060760617256, -0.991469442844391, 
                                 -0.123521246016026, 0, 0.115398973226547, 0.12756934762001, -0.985093414783478, 
                                 0, 2.20140728769467, 1.7749549958703, 26.571202402014, 1), .Dim = c(4L, 
                                                                                                     4L)), zoom = 0.699999988079071, windowRect = c(1440L, 45L, 2726L, 
                                                                                                                                                    961L))
fc.id = "Cha-F-100085"
example.cols = c(lacroix["red"],lacroix["green"],lacroix["blue"])
plot3d(inputs[fc.id], lwd = 3, col = "grey30", soma= TRUE) # show chosen neuron
plot3d(lhns.chosen["5HT1A-F-300013"], lwd = 3, col = "grey70", soma= TRUE) # show output neuron
plot3d(subset(FCWBNP.surf,"LH_R"), col = 'grey', alpha = 0.1)
rgl.snapshot(file = "images/LH_DA2_PN.png", fmt = "png")
clear3d() # now show the cut up branches
example.branches = inputs.fragments[grepl(fc.id, names(inputs.fragments))]
plot3d(example.branches, lwd = 3, soma = FALSE, col = example.cols)
plot3d(inputs[fc.id], lwd = 3, col = "grey30", soma= TRUE)
plot3d(subset(FCWBNP.surf,"LH_R"), col = 'grey', alpha = 0.1)
rgl.snapshot(file = "images/Remove_Highest_Strahler_LH_DA2_PN.png", fmt = "png")
example.input.groups = grp.input.nat.frag[names(example.branches)] # and now the corresponding NBLAST clusters
example.input.groups = example.input.groups[!is.na(example.input.groups)]
plot3d(inputs.fragments[names(grp.input.nat.frag)[grp.input.nat.frag==example.input.groups[1]]], col = lacroix["red"])
plot3d(inputs.fragments[names(grp.input.nat.frag)[grp.input.nat.frag==example.input.groups[2]]], col = lacroix["green"])
plot3d(inputs.fragments[names(grp.input.nat.frag)[grp.input.nat.frag==example.input.groups[3]]], col = lacroix["blue"])
rgl.snapshot(file = "images/Remove_Highest_Strahler_LH_DA2_PN_NBLAST_groups.png", fmt = "png")
clear3d()
plot3d(example.branches, lwd = 3, soma = FALSE)
plot3d(inputs[fc.id], lwd = 3, col = "grey30", soma= TRUE)
plot3d(subset(FCWBNP.surf,"LH_R"), col = 'grey', alpha = 0.1)
for(c in seq_along(example.input.groups)){
  g = example.input.groups[c]
  fhat = attr(ims1[[g]], "fhat")
  col = example.cols[c]
  plot(fhat,cont=c(20,70),colors=c(col,col),drawpoints=F
       ,xlab="", ylab="", zlab="",size=2, ptcol="white", add=TRUE, box=FALSE, axes=FALSE)
}
rgl.snapshot(file = "images/Remove_Highest_Strahler_LH_DA2_PN_voxels.png", fmt = "png")

# Now we can look at output regions, elsewhere in the brain
lhn.axons.frag = remove_highest_strahler(lhns.chosen, prune_function = catnat:::axonic_cable.neuron)
lhon.axons.dps = dotprops(lhn.axons.frag, OmitFailures = TRUE, resample = 0.1)
lhon.axon.morph = nblast_allbyall(lhon.axons.dps, normalisation = "mean", distance = F)
clust.nat.2 = nhclust(scoremat = lhon.axon.morph, method = "ward.D2")
fsc=dendroextras::color_clusters(clust.nat.2, k=clust.no)
plot(fsc,leaflab = "none")
grp.nat <- cutree(clust.nat.2, k = clust.no)

# Calculate voxel scores
ims2 = calculate_voxels(lhn.axons.frag, grp.nat, FCWB)
lhn.axons.marked.2 = nlapply(lhn.axons,assign_supervoxel.neuron,ims=ims2)

#######################################
### Calculate Voxel Inclusion Score ###
#######################################

# Big cell types
bigs.marked = subset(lhn.dendrites.marked.1, cell.type%in%bigs[,"cell.type"])
bigs.lns.marked = nlapply(subset(most.lhns, cell.type%in%subset(bigs,type=="LN")[,"cell.type"]),assign_supervoxel.neuron,ims=ims1)
bigs.marked = c(bigs.marked,bigs.lns.marked[setdiff(names(bigs.lns.marked),names(bigs.marked))])

# uPNs
pns.termini = subset(pns.termini,!Glomerulus%in%c("VL2pP"))
pns.marked = nlapply(pns.termini,assign_supervoxel.neuron,ims=ims1)
p.df = attr(pns.marked,"df")
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
cts[cts=="1"] = "V"
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
bigs.cols = ifelse(colnames(bigs.overlap.matrix.norm)%in%subset(bigs,type=="LN")[,"cell.type"],lacroix["green"],lacroix["blue"])
# Cluster columns
bigs.colclus =hclust(dist(t(bigs.overlap.matrix.norm),method = "euclidean"), method = "ward.D2" )
bigs.colclus = dendextend::set(as.dendrogram(bigs.colclus),"branches_lwd",3)
overlap.matrix = matrix(0,ncol=clust.no,nrow=clust.no)
LHpointssummed = colSums(do.call(rbind,lapply(lhn.dendrites.marked.1,function(x) colSums(x$d[,c(8:(7+clust.no))]))))
LHnormaliser = matrix(LHpointssummed,ncol=clust.no,nrow=clust.no, byrow = TRUE)

# Cell types
cts = unique(lhn.axons.marked.2[,"cell.type"])[!is.na(unique(lhn.axons.marked.2[,"cell.type"]))]
overlap.matrix.full = matrix(0, ncol = 25, nrow = 25)
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
  overlap.matrix.full = overlap.matrix.full + (overlap.matrix.ct/length(lhons.ct))
}
colnames(overlap.matrix.full) = paste0(1:ncol(overlap.matrix.full),"o")# Columns are letters for the output regions
rownames(overlap.matrix.full) = paste0(1:nrow(overlap.matrix.full),"i")# Columns are letters for the LH regions
grp.overlap.input = cutree(pns.rowclus, k = 7, order_clusters_as_data = FALSE)
overlap.matrix = overlap.matrix.full
rownames(overlap.matrix) = paste0(grp.overlap.input[rownames(overlap.matrix)],"i")
overlap.matrix = apply(overlap.matrix, 2, function(x) tapply(x, rownames(overlap.matrix), sum))
overlap.matrix.norm = t(apply(overlap.matrix,1,function(x)x/max(x))) # Normalise axonic output
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
    i = lacroix["darkred"]
  }else if(grepl("Food|Attractive",v)){
    i = lacroix["darkgreen"]
  }else if(grepl("Egg",v)){
    i = lacroix["darkorange"]
  }else if(grepl("Humidity",v)){
    i = lacroix["navy"]
  }else if(grepl("Phero",v)){
    i = lacroix["purple"]
  }else{
    i = "darkgrey"
  }
  col.valences = c(col.valences,i)
}
names(col.valences) = valences
mod.cols = c(lacroix["green"],lacroix["cyan"],lacroix["blue"], lacroix["pink"],lacroix["orange"],lacroix["purple"],lacroix["brown"],lacroix["red"],lacroix["yellow"])
names(mod.cols) = c("Gustatory", "Inhibitory.Olfactory", "Mechanosensation", "Memory", 
                    "Multiglomerular.Olfactory", "Neuromodulatory", "Olfactory+Gustatory", 
                    "Thermosensory", "Uniglomerular.Olfactory")
all.cols = c(col.valences, mod.cols)

# Plot a final heatmap!
heatmap.col = circlize::colorRamp2(c(0, 1), c("white", "grey70"))
pdf("images/LH_supervoxel_analysis.pdf",  width = 30, height = 12)
Heatmap(mod.inputs.termini.overlap.matrix.norm, name = "",
        cluster_rows = pns.rowclus, cluster_columns = inputs.termini.colclus, 
        column_title ="",row_title ="", col = heatmap.col,
        row_names_gp = gpar(fontsize = 0, fontface = "plain",col = lacroix["green"]),column_names_gp = gpar(fontsize = 15, fontface = "plain",col=mod.cols))+
  Heatmap(pns.overlap.matrix.norm, name = "",
          cluster_rows = pns.rowclus, cluster_columns = pns.colclus, 
          column_title ="",row_title ="", col = heatmap.col,
          row_names_gp = gpar(fontsize = 0, fontface = "plain",col = lacroix["green"]),column_names_gp = gpar(fontsize = 15, fontface = "plain",col=col.valences)) +
  Heatmap(bigs.overlap.matrix.norm, name = "",
          cluster_rows = pns.rowclus, cluster_columns = bigs.colclus, 
          column_title ="",row_title ="", col = heatmap.col,
          row_names_gp = gpar(fontsize = 20, fontface = "plain",col = lacroix["green"]),column_names_gp = gpar(fontsize = 15, fontface = "plain",col=bigs.cols))
dev.off()

# Plot Heatmaps
pdf("images/LH_input_output_analysis.pdf",  width = 5, height = 7)
Heatmap(t(overlap.matrix.norm), name = "", col = heatmap.col,
        cluster_rows = colclus, cluster_columns  = FALSE,
        column_title ="",row_title ="", 
        row_names_gp = gpar(fontsize = 15, fontface = "plain",col = lacroix["green"]),
        column_names_gp = gpar(fontsize = 15, fontface = "plain",col=lacroix["pink"]))
dev.off()

## Give a rough-type assignment to each voxel.
### On the non-normalised data
voxel.input.categories = voxel.valence = apply(mod.inputs.termini.overlap.matrix, 1, function(x) colnames(mod.inputs.termini.overlap.matrix)[which.max(x)])
voxel.pn = apply(pns.overlap.matrix.norm, 1, function(x) colnames(pns.overlap.matrix)[which.max(x)])
voxel.pn.categories = valences[voxel.pn]
names(voxel.pn.categories) = names(voxel.pn)
voxel.pn.categories[is.na(voxel.pn.categories)] = "Uniglomerular.Olfactory"
voxel.valence[grepl("\\.Olfactory",voxel.valence)] = voxel.pn.categories[grepl("\\.Olfactory",voxel.valence)]
voxel.output.correspondences = apply(overlap.matrix.full, 2, function(x) rownames(overlap.matrix.full)[which.max(x)])
voxel.output.valence = voxel.valence[voxel.output.correspondences]
names(voxel.output.valence) = names(voxel.output.correspondences)

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
  col = lacroix["green"]
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
  col = lacroix["green"]
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
  col = lacroix["pink"]
  plot(fhat,cont=c(50,95),colors=c(col,col),drawpoints=F
       ,xlab="", ylab="", zlab="",size=2, ptcol="white", add=TRUE, box=FALSE, axes=FALSE)
  rgl.snapshot(paste0("/images/supervoxels_","output_",c,".png"),fmt="png")
  clear3d()
}

## We can now make a little atlas of the lateral horn input areas.
open3d(userMatrix = structure(c(0.955159366130829, -0.285889625549316, 
                                -0.0770552754402161, 0, -0.296017527580261, -0.927844762802124, 
                                -0.226886421442032, 0, -0.00663087517023087, 0.23952242732048, 
                                -0.970868229866028, 0, -0.810615502960502, -1.20646136620167, 
                                -7.49578407532869e-06, 1), .Dim = c(4L, 4L)), zoom = 0.676839649677277, 
       windowRect = c(1460L, 65L, 2585L, 1058L))
for(n in unique(voxel.valence)){
  clusters = as.numeric(gsub("i","",names(voxel.valence)[voxel.valence==n]))
  plot3d(lhr, col = "grey", alpha = 0.1, add = TRUE)
  for(c in clusters){
    fhat = attr(ims1[[c]], "fhat")
    col = all.cols[n]
    plot(fhat,cont=c(50,95),colors=c(col,col),drawpoints=F
         ,xlab="", ylab="", zlab="",size=2, ptcol="white", add=TRUE, box=FALSE, axes=FALSE)
  }
}
rgl.snapshot(paste0("images/valence_supervoxelcluster_input.png"),fmt="png")

# And again for the output zones
open3d(userMatrix = structure(c(0.998224198818207, 0.0104750925675035, 
                                -0.0586415156722069, 0, 0.00851582363247871, -0.999400317668915, 
                                -0.0335615389049053, 0, -0.0589578822255135, 0.0330026298761368, 
                                -0.997714698314667, 0, 2.96389417437521, -12.2722030161793, 2.34399720966394e-06, 
                                1), .Dim = c(4L, 4L)), zoom = 0.436296999454498, windowRect = c(1526L, 
                                                                                                44L, 2811L, 1058L))
for(n in unique(voxel.output.valence)){
  clusters = as.numeric(gsub("o","",names(voxel.output.valence)[voxel.output.valence==n]))
  plot3d(FCWB, col = "grey", alpha = 0.1, add = TRUE)
  for(c in clusters){
    fhat = attr(ims2[[c]], "fhat")
    col = all.cols[n]
    plot(fhat,cont=c(50,95),colors=c(col,col),drawpoints=F
         ,xlab="", ylab="", zlab="",size=2, ptcol="white", add=TRUE, box=FALSE, axes=FALSE)
  }
}
rgl.snapshot(paste0("images/valence_supervoxelcluster_output.png"),fmt="png")

