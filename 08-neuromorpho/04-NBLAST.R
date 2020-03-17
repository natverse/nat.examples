# set the scene
## This script assumed that you have run the file "08-neuromorpho/00-setup.R" and "01-download/00-setup.R"

# load data
load("Jacobs_principal_neurons.rda")

# NBLAST
principals.dps = dotprops(principals, stepsize = 1, OmitFailures = TRUE)
principals.result = nblast_allbyall(principals.dps, normalisation = "mean")
k = 50
pdf("images/nat_neuromorpho_principals_nblast_dendrogram.pdf", width = 15, height = 10)
hclust=nhclust(scoremat = principals.result)
dkcs=colour_clusters(hclust, k = k) %>% set("branches_lwd", 2)
plot(dkcs)
dev.off()

# Curating the database for analysis with both t-SNE and PCA
labels <- principals.dps[,"species"]
nb = as.data.frame(principals.result)

# for plotting
colors = rainbow(length(unique(labels)))
names(colors) = unique(labels)

# The cell types we want to show
principals.species = sort(unique(principals.dps[,"species"]))
principals.colors = rainbow(length(principals.species))
names(principals.colors) = unique(principals.species)
principals.labels = labels[labels%in%principals.species]
in.principals = labels%in%principals.species

# Executing the algorithm on our data
tsne <- Rtsne(principals.result, dims = 2, perplexity=30, 
              verbose=TRUE, max_iter = 20000, pca  = TRUE, initial_dims = 50, 
              pca_center = TRUE, pca_scale = TRUE, check_duplicates = FALSE)
principals.tsne <-tsne
principals.tsne$Y <- tsne$Y[in.principals,]
principals.tsne$costs <-tsne$costs[in.principals]
principals.tsne$N <- sum(in.principals)

# quick plotting
plot(principals.tsne$Y, t='n', main="principals.tsne")
text(principals.tsne$Y, labels=principals.labels, col=principals.colors[principals.labels])

# make polygons around points of a cell type
principals.tsne.df <- data.frame(tSNE1 = principals.tsne$Y[,1], tSNE2 = principals.tsne$Y[,2], col=principals.colors[principals.labels],species = principals.labels )
tsne_plot <- ggplot(principals.tsne.df) + 
  geom_point(aes(x=tSNE1, y=tSNE2, color=principals.colors[principals.labels])) + 
  #geom_text_repel(aes(tSNE1, tSNE2, label = principals.labels)) + 
  #theme(legend.position="none") + 
  geom_polygon(aes(x = tSNE1, y = tSNE2, fill = species),alpha = 0.2)+
  scale_fill_manual(values =unique(principals.colors[principals.labels]))+
  ggtitle("")
tsne_plot
dev.off()

# Can't really tell different species apart, at this resolution for cell type
