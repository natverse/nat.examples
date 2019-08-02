### This script assumes you have run "/13-bridging/00-setup.R"

# Here, the aim is to  take neurons that belong to the primary neurite tract, PD2 (Frechter et al. 2019, eLife) from different datasets
## And find correspondences

# Let's get these different datasets
## PD2a1/b1 neurons (Dolan et al. 2018, Neuron) from different datasets, even from different labs! (they used to be called 11A in house)
pd2a1b1.mcfo = subset(lh.mcfo,old.cell.type=="11A")
pd2a1b1.mcfo.dps = dotprops(pd2a1b1.mcfo,resample=1)

pd2a1b1.fc = subset(lhns::most.lhns,grepl("PD2",cell.type)&skeleton.type=="FlyCircuit")
pd2a1b1.fc.dps = dotprops(pd2a1b1.fc,resample=1)

pd2a1b1.df = subset(lhns::most.lhns,grepl("PD2",cell.type)&skeleton.type=="DyeFill")
pd2a1b1.df.dps = dotprops(pd2a1b1.df,resample=1)

pd2a1b1.jfk = subset(jfw.lhns,grepl("PD2",cell.type))
pd2a1b1.jfk[,"skeleton.type"] = "DyeFill"
pd2a1b1.jfk.dps = dotprops(pd2a1b1.jfk,resample=1)

all.pd2ab.cands = c(em.pd2,pd2a1b1.mcfo,pd2a1b1.fc,pd2a1b1.df,pd2a1b1.jfk)
all.pd2ab.cands.dps = c(em.pd2.dps,pd2a1b1.mcfo.dps,pd2a1b1.fc.dps,pd2a1b1.df.dps,pd2a1b1.jfk.dps)
## These are from the package lhns, and are actually all already in the same brainspace for convenience, FCWB.
## In the example below, we'll put some of them back to where they came from, to show some different spaces.


# First, let's say we want to find a new EM neurons in the Mushroom body - the memory centre - of the fly
### NOTE: This example will not work unless you have access to a restricted CATMAID instance
### HOWEVER you can follow this pipleine with any other neuron you might find in the public CATMAID instance
FAFBv14.conn = catmaid_login(server = "https://neuropil.janelia.org/tracing/fafb/v14/", 
                             authname = "YOUR_USER_HERE", authpassword = "YOUR_PASSWORD_HERE", 
                             token = "YOUR_TOKEN_HERE") # Might not work for you
em.pd21 = read.neurons.catmaid("annotation:^Frechter PD2$") # Mostly draft, semi-complete neurons
em.pd22 = read.neurons.catmaid("annotation:^LH_PD2$")
em.pd2.fafb14 = c(em.pd21,em.pd22[!names(em.pd22)%in%names(em.pd21)])
em.pd2 = xform_brain(em.pd2.fafb14,sample=FAFB14,reference=FCWB,OmitFailures=TRUE)

# Simulate a less complete state of tracing by pruning the larger neurons
state = summary(em.pd2)
large = rownames(subset(state,branchpoints>50))
em.pd2[large] = nlapply(em.pd2[large],prune_strahler,orderstoprune=1:2)
em.pd2[,"cell.type"] = "EM fragment"
em.pd2[,"skeleton.type"] = "EM"
em.pd2[as.character(em.11a[,"skid"]),"cell.type"] = em.11a[,"cell.type"]
em.pd2.dps = dotprops(em.pd2,resample=1,OmitFailures = TRUE)

# Let's have a look at these different brainspaces

## FAFB14
nopen3d(userMatrix = structure(c(0.948390364646912, 0.192434668540955, 
                                 -0.252041101455688, 0, 0.172239527106285, -0.979956924915314, 
                                 -0.100092709064484, 0, -0.266250550746918, 0.0515153482556343, 
                                 -0.962526321411133, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, 
        windowRect = c(1440L, 45L, 2826L, 905L))
plot3d(FAFB14,alpha=0.1,col="grey")
zissou = colorRampPalette(wesanderson::wes_palette("Zissou1"))
em.pd2.fafb14[,"col"] = zissou(length(em.pd2.fafb14))
plot3d(em.pd2.fafb14,lwd=3,soma=T, col=em.pd2.fafb14[,"col"])
rgl.snapshot(filename = "images/FAFB14_PD2_EM_Sample.png", fmt = "png")

## JFRC2
nopen3d(userMatrix = structure(c(0.905018091201782, 0.231813445687294, 
                                 -0.356658041477203, 0, -0.0334783419966698, -0.797041535377502, 
                                 -0.602996289730072, 0, -0.424053758382797, 0.557662665843964, 
                                 -0.713576078414917, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, 
        windowRect = c(12L, 44L, 1398L, 904L))
plot3d(JFRC2,alpha=0.1,col="magenta")
em.pd2.jfrc2 = xform_brain(em.pd2,sample=FCWB,reference =JFRC2)
plot3d(em.pd2.jfrc2,lwd=3,soma=T, col=em.pd2.fafb14[,"col"])
rgl.snapshot(filename = "images/JFRC2_PD2_EM_Sample.png", fmt = "png")


# And now, let's move to use NBLAST!!


# Remove primary neurite before NBlasting
all.pd2ab.cands.no.pnt = primary.neurite(all.pd2ab.cands,keep.pnt = FALSE)
all.pd2ab.cands.no.pnt["JJ88"] = all.pd2ab.cands["JJ88"] # Already didn't have a primary neurite
all.pd2ab.cands.no.pnt.dps = dotprops(all.pd2ab.cands.no.pnt,OmitFailures = TRUE)

# NBlast all by all
pd2.nb = nblast_allbyall(all.pd2ab.cands.no.pnt.dps, UseAlpha = TRUE)
pd2.nb_scaled <- scale(pd2.nb)
pd2.nb = as.data.frame(pd2.nb)

# Curating the database for analysis with both t-SNE and PCA
labels <- all.pd2ab.cands.no.pnt[rownames(pd2.nb),"cell.type"]

# for plotting
colors = darjeeling(length(unique(labels)))[sample(1:length(unique(labels)),length(unique(labels)))]
names(colors) = unique(labels)
colors["EM fragment"] = "grey"
colors["PD2a1"] = "chartreuse4"
colors["PD2b1"] = "chartreuse2"

# Skeleton.types
skeleton.type = shapes =  all.pd2ab.cands.no.pnt[rownames(pd2.nb),"skeleton.type"]
shapes[shapes=="MCFO"] = 15 
shapes[shapes=="DyeFill"] = 17 
shapes[shapes=="FlyCircuit"] = 16 
shapes[shapes=="JFW"] = 18
shapes[shapes=="EM"] = 1
shapes = as.numeric(unique(shapes))
names(shapes) = unique(skeleton.type)

# Executing the algorithm on our data
tsne <- Rtsne(pd2.nb, dims = 2, perplexity=30, verbose=TRUE, max_iter = 20000, pca  = TRUE, initial_dims = 50, pca_center = TRUE,
              pca_scale = TRUE, check_duplicates = FALSE)

# Plotting
plot(tsne$Y, t='n', main="tsne")
text(tsne$Y, labels=labels, col=colors[labels])

# Make data frame for plotting
tsne.df <- data.frame(tSNE1 = tsne$Y[,1], tSNE2 = tsne$Y[,2], cell.type = labels, skeleton.type = skeleton.type)

# Find the hulls around our data for the different cell types
find_hull <- function(tsne.df) tsne.df[chull(tsne.df$tSNE1, tsne.df$tSNE2), ]
hulls <- ddply(tsne.df, "cell.type", find_hull)
hulls = subset(hulls,cell.type!="EM fragment")

# Ggplot
tsne_plot <- ggplot(tsne.df,aes(x=tSNE1, y=tSNE2, color= cell.type, fill = cell.type)) + 
  geom_point(aes(shape = skeleton.type)) + 
  #geom_text_repel(aes(tSNE1, tSNE2, label = labels)) + 
  #theme(legend.position="none") + 
  geom_polygon(data = hulls, aes(x=tSNE1, y=tSNE2, color= cell.type, fill = cell.type), alpha = 0.5) +
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  scale_shape_manual(values = shapes)+
  ggtitle("")
tsne_plot
