library(nat)

swcs=dir("pruned-swc")
if(!length(swcs)) stop("No traces found! See 01-download-zips.R!")

df=data.frame(name=sub("\\.swc","",basename(swcs)),
  gene=factor("lhx2a", levels=c("lhx2a", "tbx21")),
  cluster=factor(sub("([a-zA-Z]+)-.*","\\1",names(zm))),
  stringsAsFactors=FALSE)
df$gene[grep("-T",df$name)]="tbx21"
rownames(df)=df$name
zm=read.neurons("pruned-swc", neuronnames = df$name, df=df)

head(zm)

open3d()
# plot neurons from one cluster coloured by transgene
plot3d(zm, cluster=='mdG', col=gene, lwd=2, soma=T)

# plot all lhx2a neurons coloured by cluster
clear3d()
plot3d(zm, gene=='lhx2a', col=cluster, lwd=2, soma=T)
