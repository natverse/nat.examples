## This script assumed that you have run the file "03-helmstaedter2013/00-setup.R"

# There is presently a problem reading the main Matlab skeleton dataset into R
# because it is uses advanced MAT format features that are not supported by
# the R.matlab package. Therefore I have provided a canned download for the time
# being

urls=file.path("http://flybrain.mrc-lmb.cam.ac.uk/si/nat/helmstaedter/",
  c('sk.uniq.rda','skeleton_metadata.rda'))

message("Checking for presence of data (30Mb) if necessary ...")
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
