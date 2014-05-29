# Download metadata
###################

if(!file.exists('flycircuit.boundingbox.txt')){
  message("Downloading boundingbox metadata for tracings")
  download.file("http://www.flycircuit.tw/flycircuit.boundingbox","flycircuit.boundingbox.txt")
}
message("Processing metadata")

kcs20_neuron_names=c("fru-M-500112", "Gad1-F-900005", "Gad1-F-100010", "Gad1-F-300043", 
"fru-F-400045", "fru-F-300059", "fru-M-100078", "Gad1-F-300107", 
"Cha-F-300150", "fru-M-300145", "Gad1-F-400089", "fru-F-400017", 
"fru-F-000031", "fru-M-400058", "fru-F-200098", "fru-F-200021", 
"Gad1-F-300023", "fru-M-100014", "Gad1-F-700033", "fru-F-400181")

bbl=readLines("flycircuit.boundingbox.txt")
bbl=gsub(",","",bbl,fixed=TRUE)
names(bbl)=sub("^([^ ]+) .*","\\1",bbl)
# restrict to just these neurons
bbl=bbl[kcs20_neuron_names]
flycircuit.md=read.table(text=bbl,
	col.names=c("Neuron",'define','Lattice','LX','LY','LZ','bb','X1','X2','Y1','Y2','Z1','Z2'))

flycircuit.md=flycircuit.md[!names(flycircuit.md)%in%c("define","Lattice",'bb')]
flycircuit.md=transform(flycircuit.md,
  dx=(X2-X1)/(LX-1),dy=(Y2-Y1)/(LY-1),dz=(Z2-Z1)/(LZ-1))
flycircuit.md$gene_name=kcs20_gene_names
rownames(flycircuit.md)=kcs20_gene_names

flycircuit.md$driver=substr(flycircuit.md$gene_name,1,3)

# some information that I gleaned about the original lsm files
lsminfo=read.table("lsminfo.txt", header=TRUE)
flycircuit.md=cbind(flycircuit.md,lsminfo[rownames(flycircuit.md),])
save(flycircuit.md,file='flycircuit.md.rda')
