# see 
# http://www.nature.com/ncomms/2014/140324/ncomms4512/full/ncomms4512.html#supplementary-information

message("Downloading supplementary data from Nature Comms website (~70Mb)")
urls=sprintf('http://www.nature.com/ncomms/2014/140324/ncomms4512/extref/ncomms4512-s%d.zip',c(2,13,14))
# alternative download location for third file (s14.zip->S13.mat)
urls=c(urls,'http://zenodo.org/record/10737/files/ncomms4512-s14.zip')
for (url in urls){
	localfile=basename(url)
	if(file.exists(localfile)) next
	t=try(download.file(url, localfile, mode='wb'))
	if(inherits(t,'try-error')) {
		# remove any bad download
		unlink(localfile)
		message("unable to download ", url)
		next
	}
	unzip(localfile)
}
