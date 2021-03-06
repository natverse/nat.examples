if(!require('nat')) install.packages("nat")
if(!require('nat.nblast')) install.packages("nat.nblast")
if(!require('dendextend')) install.packages("dendextend")
if(!require('R.matlab')) install.packages("R.matlab")
if(!require('Morpho')) install.packages("Morpho")
if(!require('alphashape3d')) install.packages("alphashape3d")
if(!require('Rtsne')) install.packages("Rtsne")

library(nat)
library(nat.nblast)
library(dendextend)
library(Morpho)
library(alphashape3d)
library(R.matlab)
library(Rtsne)

# set working directory to location of this file
try(setwd(dirname(attr(body(function() {}),'srcfile')$filename)))

#' plot3d method for Moritz's neurons
#'
#' @param sk A neuron in Moritz's format (defined by matlab data)
#' @param ... Additional graphical parameters passed to segments3d
#' @detail assumes very little about the data - just that there are nodes in 
#'   Nx3 arrangement and edges that connect thode nodes in Ex2 arrangement.
plot3d.skel<-function(sk,col,WithNodes=FALSE,col.node='black',...){
  starts=sk$nodes[sk$edges[,1],1:3]
  stops=sk$nodes[sk$edges[,2],1:3]
  interleaved=matrix(t(cbind(starts,stops)),ncol=3,byrow=T)
  if(missing(col)) col=rainbow
  if(is.function(col)){
    g=as.ngraph(sk)
    nc=igraph::no.clusters(g)
    col=col(nc)
    cg=igraph::clusters(g)
    origvids=igraph::get.vertex.attribute(g,'label')
    startcols=col[cg$membership[match(sk$edges[,1],origvids)]]
    endcols=col[cg$membership[match(sk$edges[,2],origvids)]]
    col=as.vector(rbind(startcols,endcols))
  }
  rval=list()
  rval$segments=segments3d(interleaved,col=col,...)
  if(WithNodes){
    allnodes=unique(c(sk$edges[,1],sk$edges[,2]))
    rval$nodes=points3d(sk$nodes[allnodes,1:3],col=col.node,...)
  }
  
  invisible(segments3d(interleaved,col=col,...))
}

`xyzmatrix<-.skel`<-function(n, value, ConnectedOnly=FALSE, ...){
  if(ConnectedOnly) stop("Not yet implemented")
  if(ncol(value)!=3) stop("Expects a Nx3 matrix")
  if(nrow(value)!=nrow(n$nodes))
    stop("Replacement data must have same number of rows as existing points")
  n$nodes[,1:3]=value
  n
}

xyzmatrix.skel<-function(x,ConnectedOnly=FALSE,Transpose=FALSE,...) {
  if(ConnectedOnly) stop("Not yet implemented")
  mx=data.matrix(x$nodes[,1:3, drop=FALSE])
  if(Transpose) t(mx) else mx
}

#' Arithmetic for neuron coordinates
#'
#' If x is one number 
#' If x is a 3-vector, multiply xyz only
#' TODO Figure out how to document arithemtic functions in one go
#' @param n a neuron
#' @param x (a numeric vector to multiply neuron coords in neuron)
#' @return modified neuron
#' @export
#' @examples
#' n1<-MyNeurons[[1]]*2
#' n2<-MyNeurons[[1]]*c(2,2,2,2)
#' stopifnot(all.equal(n1,n2))
#' n3<-MyNeurons[[1]]*c(2,2,4)
`*.skel` <- function(n,x) {
	# TODO look into S3 generics for this functionality
	if(!is.numeric(x))
		stop("expects a numeric vector")
	lx=length(x)
	if(lx==1) xyzmatrix(n)<-xyzmatrix(n)*x
	else if(lx==3) xyzmatrix(n)<- t(t(xyzmatrix(n))*x)
	else stop("expects a numeric vector of length 1 or 3")
	n
}

`+.skel` <- function(n,x) {
	if(!is.numeric(x))
		stop("expects a numeric vector")
	lx=length(x)
	if(lx==1) xyzmatrix(n)<-xyzmatrix(n)+x
	else if(lx==3) xyzmatrix(n)<- t(t(xyzmatrix(n))+x)
	else stop("expects a numeric vector of length 1 or 3")
	n
}

`-.skel` <- function(n,x) n+(-x)
`/.skel` <- function(n,x) n*(1/x)

#' Convert Helmstaedter's matlab skel format into nat::neuron objects
#' 
#' @description skel objects are my direct R translation of the matlab data
#'   provided by Briggman, Helmstaedter and Denk.
#' @param x A skel format neuron to convert
#' @param ... arguments passed to as.neuron
#' @return An object of class \code{neuron}
#' @seealso \code{\link{as.neuron}}
as.neuron.skel<-function(x, ...) {
  as.neuron(as.ngraph(x), ...)
}

#' Convert Helmstaedter's matlab skel format into nat::ngraph objects
#' 
#' @details \code{ngraph} objects are thin wrappers for \code{igraph::graph} 
#'   objects.
#'   
#'   This function always removes self edges (i.e. when a node is connected to 
#'   itself), which seem to be used as a placeholder in the Helmstaedter dataset
#'   when there are no valid edges.
#'   
#'   Isolated nodes are also removed by default (i.e. the nodes that are not
#'   connected by any edges.)
#' @param remove.isolated Whether or not to remove isolated nodes from the graph
#'   (see Details).
#' @inheritParams as.neuron.skel
#' @seealso \code{\link[nat]{ngraph}}, \code{\link{as.neuron.skel}}
as.ngraph.skel<-function(x, remove.isolated=TRUE, ...) {
  self_edges=x$edges[,1]==x$edges[,2]
  if(sum(self_edges)>0){
    if(sum(self_edges)==nrow(x$edges)){
      # there are only self edges - make a dummy empty edge matrix
      x$edges=matrix(nrow=0, ncol=2)
    } else {
      x$edges[!self_edges, , drop=FALSE]
    }
  }
  ng=ngraph(x$edges, vertexlabels = seq_len(nrow(x$nodes)), xyz=x$nodes, ...)
  if(remove.isolated){
    isolated_vertices=igraph::V(ng)[igraph::degree(ng)==0]
    g=igraph::delete.vertices(graph=ng,isolated_vertices)
    ng=as.ngraph(g)
  }
}

# Load data from previous session
# load(".RData")
