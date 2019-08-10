# install
if(!require('devtools')) install.packages("devtools")
if(!require('ks')) install.packages("ks")
if(!require('tidyverse')) install.packages("tidyverse")
if(!require('lhns')) devtools::install_github('jefferislab/lhns')
if(!require('natverse')) devtools::install_github('natverse/natverse')
if(!require('ComplexHeatmap')) install.packages('ComplexHeatmap')

# load
library(tidyverse)
library(lhns)
library(natverse)
library(alphashape3d)
library(ks)
library(ComplexHeatmap)
library(catnat)

# some nice colors!! Inspired by LaCroixColoR
lacroix = c("#C70E7B", "#FC6882", "#007BC3", "#54BCD1", "#EF7C12", "#F4B95A", 
            "#009F3F", "#8FDA04", "#AF6125", "#F4E3C7", "#B25D91", "#EFC7E6", 
            "#EF7C12", "#F4B95A", "#C23A4B", "#FBBB48", "#EFEF46", "#31D64D", 
            "#132157","#EE4244", "#D72000", "#1BB6AF")
names(lacroix) = c("purple", "pink",
                   "blue", "cyan",
                   "darkorange", "paleorange",
                   "darkgreen", "green",
                   "brown", "palebrown",
                   "mauve", "lightpink",
                   "orange", "midorange",
                   "darkred", "darkyellow",
                   "yellow", "palegreen", 
                   "navy","cerise",
                   "red", "marine")

# set working directory to location of this file
try(setwd(dirname(attr(body(function() {}),'srcfile')$filename)))

# Another extra function that we shall need
subtree <- function(neuron, subtree = 1){
  if(neuron$nTrees<1){
    warning("Neuron has no subtree")
  }else{
    v = unique(unlist(neuron$SubTrees[subtree]))
    neuron = nat::prune_vertices(neuron, verticestoprune = v, invert = TRUE)
  }
  neuron
}

# assign supervoxel
assign_supervoxel.neuron <- function(neuron, ims){
  df = cbind(neuron$d, matrix(0,nrow(neuron$d),length(ims)))
  for (im in 1:length(ims)){
    i = ims[[im]]
    idxs = c()
    for (point in 1:nrow(xyzmatrix(neuron))){
      idx = tryCatch(coord2ind(xyzmatrix(neuron)[point,], imdims = i), error = function(e) 0) 
      idxs = c(idxs,idx)
    }
    #if (idxs[1] == 0){
    #ids = rep(0,nrow(xyzmatrix(neuron)))
    ids=i[idxs]
    change.rows = 1:nrow(df)
    if(length(length(which(idxs==0)))>0){
      change.rows = change.rows[!change.rows%in%which(idxs==0)]
    }
    df[change.rows,as.character(im)] = ids#/max(i) # Normalise within voxel
  }
  df = df[,-c(1:ncol(neuron$d))]#/max(df[,-c(1:ncol(neuron$d))]) # Normalise by neuron
  df = cbind(df, max = apply(df[,-c(1:ncol(neuron$d))],1,function(x) which.max(x)))
  neuron$d = cbind(neuron$d,df)
  neuron$cluster = unique(df$max)[which.max(tabulate(match(df$max, unique(df$max))))]
  neuron
}

# decompose by Strahler
remove_highest_strahler <-  function(someneuronlist, 
                                     brain = FCWBNP.surf, 
                                     neuropil = NULL,
                                     prune_function = NULL){
  neurons.strahler = assign_strahler(someneuronlist, OmitFailures = TRUE, .progress = FALSE)
  neurons.fragments = neuronlist()
  for(i in 1:length(neurons.strahler)){
    id = names(neurons.strahler[i])
    x = neurons.strahler[id][[1]]
    if(!is.null(neuropil)){
      p = tryCatch( prune_in_volume(x, brain = FCWBNP.surf, neuropil = neuropil), error = function(e) NULL) # Prune to just the lateral horn
    }else{
      p=x
    }
    if(!is.null(prune_function)){
      p = tryCatch( prune_function(x), error = function(e) NULL) # Prune neuron further
    }
    if(is.null(p)){
      next
    }
    relevant.points = subset(x$d, PointNo%in%p$d$PointNo)
    p$d$strahler_order = relevant.points[match(p$d$PointNo,relevant.points$PointNo),]$strahler_order
    s = max(p$d$strahler_order)
    if(s>1){
      v = subset(p$d,strahler_order%in%c(s))[,"PointNo"] # remove highest Strahler order branch
      b = tryCatch(nat::prune_vertices(x, verticestoprune = c(v, setdiff(x$d$PointNo,p$d$PointNo)), invert = FALSE), error = function(e) NULL)
      if(is.null(b)){
        next
      }
      for(t in 1:b$nTrees){
        if(length(unlist(b$SubTrees[t]))>1){
          subt = subtree(b, subtree = t)
          subt = as.neuronlist(subt)
          attr(subt,"df") = df
          names(subt) = paste0(id,"_",t)
          neurons.fragments = c(neurons.fragments,subt)  
        }
      }
    }else{
      p = as.neuronlist(p)
      attr(p,"df") = df
      names(p) = paste0(id,"_",t)
      neurons.fragments = c(neurons.fragments,p)
    }
  }
  neurons.fragments
}

# Create continuous anatomical supervoxels
calculate_voxels <- function(someneuronlist, clusters, volume){
  ims = list()
  for(n in sort(unique(clusters))){
    l = n
    selected = someneuronlist[names(someneuronlist)%in%names(clusters)[clusters==n]]
    points = xyzmatrix(someneuronlist[names(someneuronlist)%in%names(clusters)[clusters==n]])
    H.pi <- Hpi(points,binned=TRUE) 
    fhat <- kde(points, H=H.pi, xmax = boundingbox(volume)[2,],xmin = boundingbox(volume)[1,], w = rep(1/nrow(points),nrow(points))) # Area under curve = 1
    i = im3d(fhat$estimate,BoundingBox=boundingbox(volume))
    attr(i,"contours") = fhat$cont
    attr(i,"fhat") = fhat
    ims[[n]] <- i
  }
  names(ims) <- 1:length(ims)
  ims
}
