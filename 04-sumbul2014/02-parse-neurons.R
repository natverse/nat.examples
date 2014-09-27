library(nat)
library(R.matlab)

message("Reading matlab data")
matlab_files=dir(patt="S[0-9]+.mat")
matlab_data=lapply(matlab_files, function(f) readMat(f)[[1]])
# nb matlab data comes in as list of lists but want to make into one long list
matlab_data=unlist(matlab_data, recursive = FALSE)

#' convert Sumbul matlab data describing neuron structure to nat::neuron object
#' 
#' @details This is how Sumbul et al describe their data format
#' arborTraces1{k}.warpedArborNodes: 3d positions of the nodes of the arbor trace
#' arborTraces1{k}.edges: node pairs defining each edge in the tree structure
#' arborTraces1{k}.id: identity of the trace, including the neuron's and the tracer's name
#' arborTraces1{k}.cluster: structural cluster the cell was assigned to in the paper
#' arborTraces1{k}.geneticLine: name of the transgenic line the cell was obtained from
sumbul2neuron<-function(x){
	points=x[[1]]
	# nat uses the convention root->child for edge directionality whereas they
	# seem to have used the reverse
	edges=x[[2]][,2:1]
  # construct a graph object from the information 
  # ngraph is a slightly customised version of an igraph object
	ng=nat::ngraph(edges,vertexlabels=seq(nrow(points)),xyz=points)
	
	as.neuron(ng,NeuronName=unlist(x['id',,]), cluster=unlist(x['cluster',,]),
	geneticLine=unlist(x['geneticLine',,]))
}

message("Converting matlab data to nat's neuron format. This could take a couple of mins")
sumbuln=nlapply(matlab_data, sumbul2neuron)

# Now make a dataframe with useful information

df=data.frame(id=sapply(sumbuln,'[[','NeuronName'), 
  cluster=sapply(sumbuln,'[[','cluster'),
  geneticLine=sapply(sumbuln,'[[','geneticLine'), 
  nTrees=sapply(sumbuln,'[[','nTrees'),
  NumSegs=sapply(sumbuln,'[[','NumSegs'),
  NumPoints=sapply(sumbuln,'[[','NumPoints'),
  StartPoint=sapply(sumbuln,'[[','StartPoint')
)

df$id=as.character(df$id)
df$tracer=factor(sub(".*_([^_]+)","\\1",df$id))
df$cell=sub("^([^_]+)_.*","\\1",df$id)

soma<-function(x) data.matrix(x$d[x$StartPoint,c("X",'Y','Z')])
somapos=t(sapply(sumbuln,soma))
colnames(somapos)=c("X",'Y','Z')
df=cbind(df,somapos)

attr(sumbuln,'df')=df

save(sumbuln,file='sumbuln.rda')

plot3d(sumbuln[1:10],soma=4,WithNodes=F)
