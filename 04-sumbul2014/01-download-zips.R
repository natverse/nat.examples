# see 
# http://www.nature.com/ncomms/2014/140324/ncomms4512/full/ncomms4512.html#supplementary-information

message("Downloading supplementary data from Nature Comms website (~70Mb)")
urls=sprintf('http://www.nature.com/ncomms/2014/140324/ncomms4512/extref/ncomms4512-s%d.zip',c(2,13,14))
for (url in urls){
	localfile=basename(url)
	if(file.exists(localfile)) next
	t=try(download.file(url,localfile))
	if(inherits(t,'try-error')) {
		message("unable to download ", url)
		if(grepl("14\\.zip",url)){
			# fetch from alternative location
			url='http://zenodo.org/record/10737/files/ncomms4512-s14.zip'
			download.file(url,localfile)
		}
		next
	}
	unzip(localfile)
}
