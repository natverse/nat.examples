# install
if(!require('devtools')) install.packages("devtools")
if(!require('tidyverse')) install.packages("tidyverse")
if(!require('natverse')) devtools::install_github("natverse/natverse")
if(!require('lhns')) devtools::install_github("jefferislab/lhns")
if(!require('ks')) install.packages("ks")
if(!require('lhns')) devtools::install_github('jefferislab/lhns')
if(!require('alphashape3d')) install.packages('alphashape3d')
if(!require('ComplexHeatmap')) install.packages('ComplexHeatmap')

# load natversee functions
library(natverse)
library(tidyverse)
library(lhns)
library(alphashape3d)
library(ks)
library(ComplexHeatmap)
## What is lhns? It is a little data package for lateral horn neurons, made for Frechter et al. 2019, eLife.
## It contains some data we will play with in the next R file

# set working directory to location of this file
try(setwd(dirname(attr(body(function() {}),'srcfile')$filename)))

# load data from previoous sessions
if (exists(".RData")){
    load(".RData")
}

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
