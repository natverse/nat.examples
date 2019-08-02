# set the scene
## This script assumed that you have run the file "03-helmstaedter2013/00-setup.R" and "03-helmstaedter2013/01-download.R"

# get from source
library(R.matlab)
urls=file.path("http://flybrain.mrc-lmb.cam.ac.uk/si/nat/helmstaedter/",
               c('kn_e2006_ALLSKELETONS_FINAL2012.mat'))
skall = readMat(urls)

# or read from saved
skall=readMat('03-helmstaedter2013/kn_e2006_ALLSKELETONS_FINAL2012.mat')


message("Checking for presence of data (160Mb) if necessary ...")
for (url in urls){
  localfile=basename(url)
  if(file.exists(localfile)) next
  message("Downloading data ...")
  t=try(download.file(url, localfile, mode='wb'))
  if(inherits(t,'try-error')) {
    message("unable to download ", url)
    next
  }
}

#' Parse a single skeleton object cleaning up the rather messy structure
#' that comes out of readMat
parse.moritz.skel<-function(x){
  # delist
  x=x[[1]]
  stopifnot(inherits(x,'array'))
  
  vars=rownames(x)
  stopifnot(all(c('nodes','edges') %in% vars))
  
  simplevars=intersect(c('nodes','edges', 'edgeSel'), vars)
  r=sapply(simplevars, function(v) x[v,,][[1]], simplify=FALSE)
  
  othervars=setdiff(vars, simplevars)
  process_var<-function(y){
    if(inherits(y,'list')) y=y[[1]]
    if(inherits(y,'array')) apply(y, 1, unlist) else {
      if(is.numeric(y)) drop(y) else y
    }
  }
  
  r2=sapply(othervars, function(v) process_var(x[v,,]), simplify = FALSE)
  structure(c(r, r2), class=c('skel','list'))
}


# read raw matlab data
library(R.matlab)
skall=readMat('kn_e2006_ALLSKELETONS_FINAL2012.mat')

# Convert all the neurons to intermediate skeleton format
skallp=nlapply(skall$kn.e2006.ALLSKELETONS.FINAL2012, parse.moritz.skel, .progress='text')
# give the neurons names
names(skallp)=sprintf("sk%04d", seq_along(skallp))

# make a data.frame of useful metadata
df=as.data.frame(sapply(skeleton_metadata[1:3], drop, simplify = FALSE), 
                 row.names=names(skallp))
names(df)=sub("kn.e2006.ALLSKELETONS.FINAL2012.","",names(df), fixed = T)
data.frame(skallp)=df

# convert to nat's neuron representation
skalln=nlapply(skallp, as.neuron, OmitFailures = TRUE, .progress='text')
# save in native R format
save(skalln, file='skalln.rda')
# save a zip archive of SWC format neurons for all reconstructions
write.neurons(skalln, dir='skalln.swc.zip', files=names(skalln), format='swc')

