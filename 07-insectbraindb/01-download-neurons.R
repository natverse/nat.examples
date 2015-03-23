# Defines a function to read neurons from insectbraindb website
# See http://www.tedore.net/

#' read a neuron from the insectbraindb/tedore.net site
#' @param x Either a neuron number or a full url to the summary page of a neuron
#' @param baseurl base url for the insectbraindb website
#' @param ... additional arguments passed to read.neuron
#' @examples
#' nn=read.neurons()
read.tedore.neuron<-function(x, baseurl="http://www.tedore.net/", ...) {
	# complete url for input spec if required
	if(is.numeric(x))
		x=file.path(file.path(baseurl,"neurons"),x,"")

	# get all possible html nodes
	h=rvest::html(x)
	nn=rvest::html_nodes(h,'.dl-horizontal')

	# subset node with 3d model info
	n3d=grepl("3-D model",html_text(nn) )
	if(!any(n3d))
		stop("No 3d model available")

	# now find url for the 3d model
	n=nn[[which(n3d)]]
	url=html_attr(html_nodes(n,'a'),'href')
	url=file.path(baseurl, url)
	if(!length(url))
		stop("Unable to find unique url for 3d neuron model")
	
	read.neuron(url, ...)
}

# table of information about available neurons
neurons_table=html_table(html("http://www.tedore.net/neurons/"), header=T)[[1]]
# tidy up column / row names for R
names(neurons_table)=make.names(names(neurons_table))
rownames(neurons_table)=neurons_table$Identification.Number
# make species a factor


message("Downloading ", nrow(neurons_table), " neurons")
tedoren=nlapply(neurons_table$Identification.Number, read.tedore.neuron, .progress='text', OmitFailures=F)
names(tedoren)=neurons_table$Identification.Number

# attach data.frame then keep only good neurons / metadata
attr(tedoren, 'df')=neurons_table
actual_neurons=sapply(tedoren, is.neuron)
message("Keeping ", length(actual_neurons)," valid neurons!")
tedoren=tedoren[actual_neurons]

# make species a factor (since there are only a few distinct values)
attr(tedoren, 'df')$Species=factor(attr(tedoren, 'df')$Species)
save(tedoren, file="tedore_neurons.rda")
