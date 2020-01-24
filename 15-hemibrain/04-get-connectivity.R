# This script assumed that you have run the file "15-hemibrain/00-setup.R" and "15-hemibrain/02-find-neurons.R"

## We are looking at connectivity between MBONs.
### We can get an adjaceny matrix between all of these MBONs
mbon.adj = neuprint_get_adjacency_matrix(bodyids = mbon.info$bodyid)
rownames(mbon.adj) = colnames(mbon.adj) = mbon.info$type

## Let's visualise this connectivity
my_palette <- colorRampPalette(c(lacroix[["cyan"]],lacroix[["yellow"]], lacroix[["orange"]], lacroix[["cerise"]])) (n=20)
pdf(file = "images/MBON_interconnectivity.pdf", height = 10, width = 10)
Heatmap(mbon.adj,
        col = my_palette)
dev.off()

## It is hard to parse this information mentally
### What if we just see what the compartment to compartment connectivity is?
mbon.adj.comp = mbon.adj
rownames(mbon.adj.comp) = colnames(mbon.adj.comp) = mbon.info$compartment
mbon.adj.comp = t(apply(t(mbon.adj.comp), 2, function(x) tapply(x, colnames(mbon.adj.comp), mean, na.rm = TRUE)))
mbon.adj.comp = apply(mbon.adj.comp, 2, function(x) tapply(x, rownames(mbon.adj.comp), mean, na.rm = TRUE))
pdf(file = "images/MBON_compartment_interconnectivity.pdf", height = 10, width = 10)
Heatmap(mbon.adj.comp,
        col = my_palette)
dev.off()

## What about MBON connectivity to neurons not in this set of MBONs?
### First, is there a common partner to all MBONs?
common = neuprint_common_connectivity(mbon.info$bodyid, prepost = "PRE") # could use 'post' for downstream
dim(common) # nope!
### What about all a'3 MBONs?
ap3.mbon.info = subset(mbon.info, grepl("a'3",name))
ap3.common = neuprint_common_connectivity(ap3.mbon.info$bodyid, prepost = "PRE")
dim(ap3.common) # so there are 25 common partners to the 4 MBONs
ap3.common.meta = neuprint_get_meta(colnames(ap3.common)) # what are they?
View(ap3.common.meta) # ooh of course, the dopamine input (DAN) to a'3 and other weird things, but no KCs!

## We can also get all of the neurons in the database that connect to the
### query neurons, either upstream or downstream of them
ap3.mbon.connected = neuprint_connection_table(ap3.mbon.info$bodyid, prepost = "POST")
### In which brain region are these partners?
table(ap3.mbon.connected$roi)

## Let's have a look at what the strongest downstream partners look like
ap3.mbon.connected.strong = subset(ap3.mbon.connected,weight>30, prepost = "POST")
ap3.targets = neuprint_read_neurons(ap3.mbon.connected.strong$partner)
ap3.targets = unspike(ap3.targets, threshold=1000)
nopen3d(userMatrix = structure(c(0.964227318763733, -0.0444099828600883, 
                                 -0.261329621076584, 0, 0.249727189540863, -0.178416073322296, 
                                 0.951737523078918, 0, -0.088891975581646, -0.982952415943146, 
                                 -0.160943150520325, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.676839649677277, 
        windowRect = c(1440L, 45L, 2937L, 938L))
plot3d(ap3.targets, lwd = 2)
rgl.snapshot(filename = "images/hemibrain_ap3_strong_downstream_partners.png", fmt ="png")
View(ap3.targets[,])

## Can other MBONs influence these strong partners?
### The cognate DAN for MBONs-a'3, PPL1-a'3, gets 49 (!) of its inputs from MBONs-a'3
### Which other MBONs impinge on it?
### Via whatever path?
shortest.paths = data.frame()
for(b in mbon.info$bodyid){
  sp = neuprint_get_shortest_paths(body_pre = b, body_post = "517514142")
  dupe = ifelse( is.na(which(duplicated(sp$to))[1]),nrow(sp),which(duplicated(sp$to))[1])
  shortest = sp[1:dupe,]
  if(nrow(shortest)>0 & ncol(shortest) > 0 & nrow(shortest)<5 ){
    shortest$order.to = paste(1:nrow(shortest), nrow(shortest), sep = "/")
    shortest$order.from = paste(1:nrow(shortest)-1, nrow(shortest), sep = "/")
    shortest[1,"order.from"] = mbon.info[as.character(b),"compartment"]
    shortest[nrow(shortest),"order.to"] = ap3.targets["517514142","type"]
    shortest.paths = rbind(shortest.paths, shortest) 
  }
}

# Make network
set.seed(42)
paths = aggregate(list(weight = shortest.paths$weight), list(order.from = shortest.paths$order.from,
                                                             order.to = shortest.paths$order.to),
                           sum)
n = network(paths,
            matrix.type = "edgelist",
            ignore.eval = FALSE,
            layout = "fruchtermanreingold",
            names.eval = "weight",
            directed = TRUE)
n = ggnetwork(n, cell.jitter = 0.75, arrow.gap = 0.01)

# Set colours
orders = unique(c(paths$order.to,paths$order.from))
order.cols = rep("grey30",length(orders))
names(order.cols) = orders
comps = unique(mbon.info$compartment)
mbon.cols = rep(lacroix[["pink"]],length(comps))
names(mbon.cols) = comps
mbon.cols = c(mbon.cols,
              order.cols[!names(order.cols)%in%c("PPL104",names(mbon.cols))],
              "PPL104" = lacroix[["orange"]])

# Plot
set.seed(1)
pdf(file = "images/PPL1ap3_MBON_input_graph.pdf", height = 5, width = 5)
ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(color = vertex.names),
             curvature = 0.05,
             arrow = arrow(length = unit(6, "pt"), 
                           type = "closed")) +
  geom_nodes(aes(color = vertex.names, size = 6)) +
  geom_edgetext(aes(label = weight, color = vertex.names), fill = NA) +
  geom_nodelabel_repel(aes(color = vertex.names, label = vertex.names),
                       fontface = "bold", box.padding = unit(1, "lines")) +
  scale_color_manual(values = mbon.cols) +
  scale_fill_manual(values = mbon.cols) +
  theme_blank() +  
  guides(color = FALSE, shape = FALSE, fill = FALSE, size = FALSE, linetype = FALSE) + ylab("") + xlab("")
dev.off()
