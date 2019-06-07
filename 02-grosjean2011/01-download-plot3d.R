# Download packaged projection neuron dataset accompanying
# An olfactory receptor for food-derived odours promotes male courtship in Drosophila.
# Grosjean Y, Rytz R, Farine JP, Abuin L, Cortot J, Jefferis GS, Benton R
# Nature478p236-40(2011 Sep 28)

# Load the data
load(url("http://flybrain.mrc-lmb.cam.ac.uk/si/grosjean11/MyNeuronsFCIR.rda"))
# the data object is a list of 300 neurons in a nat neuronlist object
class(MyNeurons)
class(MyNeurons[[1]])

# Open an rgl 3d display
open3d()

# Plot a single neuron
plot3d(MyNeurons[[1]])

# Clear 3d display and plot all neurons
plot3d(MyNeurons,WithNode=FALSE)
# Take a picture!
rgl.snapshot(filename ="images/nat_Grosjean_Olfactory_PNs.png" ,fmt = "png")

# Look at the metadata attached to MyNeurons
head(MyNeurons)
# you can get the whole dataframe like this
PNdf=attr(MyNeurons,'df')
# or access attached dataframe columns like this
with(MyNeurons, table(Glomerulus))

# 3d plot of neurons from olfactory glomeruli beginning DM
# coloured by glomerulus
nopen3d(userMatrix = structure(c(1, 0, 0, 0, 0, 0.342020143325668, 
                                 -0.939692620785909, 0, 0, 0.939692620785909, 0.342020143325668, 
                                 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.530321657657623, 
        windowRect = c(1480L, 85L, 3177L, 1055L)) # View neurons from a specific position and zoom, in 3D
rval=plot3d(MyNeurons, subset=grepl("^DM",Glomerulus), col=factor(Glomerulus),lwd=2, WithNodes=FALSE)
rgl.snapshot(filename ="images/nat_Grosjean_Olfactory_PNs.png" ,fmt = "png")
# make a legend so that you know which colours match which glomerulus
pdf("images/nat_Grosjean_Olfactory_PNs_key.pdf", width = 5, height = 5)
plot.new()
with(attr(rval,'df'), legend('center', legend = unique(Glomerulus), fill=unique(col)))
dev.off()

# more help for commands we have used
?plot3d.neuron
?plot3d.neuronlist
?with.neuronlist

# some other interesting commands
?find.neuron
?write.neurons
