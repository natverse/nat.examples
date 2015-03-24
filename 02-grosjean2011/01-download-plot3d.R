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
clear3d()
plot3d(MyNeurons,WithNode=FALSE)

# Look at the metadata attached to MyNeurons
head(MyNeurons)
# you can get the whole dataframe like this
PNdf=attr(MyNeurons,'df')
# or access attached dataframe columns like this
with(MyNeurons, table(Glomerulus))

# 3d plot of neurons from olfactory glomeruli beginning DM
# coloured by glomerulus
clear3d()
rval=plot3d(MyNeurons, subset=grepl("^DM",Glomerulus), col=factor(Glomerulus),
  lwd=2, WithNodes=FALSE)
# make a legend so that you know which colours match which glomerulus
with(attr(rval,'df'), legend('center', legend = unique(Glomerulus), fill=unique(col)))

# more help for commands we have used
?plot3d.neuron
?plot3d.neuronlist
?with.neuronlist

# some other interesting commands
?find.neuron
?write.neurons
