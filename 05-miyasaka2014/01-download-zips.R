urls <- "http://flybrain.mrc-lmb.cam.ac.uk/si/nat/miyasaka/pruned-traces.zip"

message("Downloading data (1.6 MB) if necessary...")
for (url in urls){
  localfile <- basename(url)
  if(file.exists(localfile)) next
  message("Downloading ", url, "...")
  t <- try(download.file(url, localfile, mode='wb'))
  if(inherits(t, 'try-error')) {
    message("Unable to download ", url)
    next
  }
}

message("Unzipping")
unzip("pruned-traces.zip")
