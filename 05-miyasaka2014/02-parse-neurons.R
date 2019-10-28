library(nat)
## This script assumed that you have run the file "05-miyasaka2014/00-setup.R", and the subsequent R files

swcs=dir("pruned-swc")
if(!length(swcs)) stop("No traces found! See 01-download-zips.R!")

df=data.frame(name=sub("\\.swc","",basename(swcs)),
  gene=factor("lhx2a", levels=c("lhx2a", "tbx21")),
  stringsAsFactors=FALSE)
df$cluster=factor(sub("([a-zA-Z]+)-.*","\\1",df$name))
df$gene[grep("-T",df$name)]="tbx21"
rownames(df)=df$name
zm=read.neurons("pruned-swc", neuronnames = df$name, df=df)

head(zm)

open3d()
# plot neurons from one cluster coloured by transgene
plot3d(zm, cluster=='mdG', col=gene, lwd=2, soma=T)

# plot all lhx2a neurons coloured by cluster
nopen3d(userMatrix = structure(c(0.999464690685272, 0.0153909111395478, 
                                 -0.0288552921265364, 0, -0.0159958377480507, 0.999654471874237, 
                                 -0.0208465103060007, 0, 0.0285244900733232, 0.0212973672896624, 
                                 0.999365568161011, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.746215641498566, 
        windowRect = c(1440L, 45L, 2942L, 1063L))
plot3d(zm, gene=='lhx2a', col=cluster, lwd=2, soma=T)
rgl.snapshot(filename ="images/nat_zfishOB_neurons.png" ,fmt = "png")

# Save all neurons
saveRDS(zm, file='neurons.rds')
