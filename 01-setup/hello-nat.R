# install
install.packages("nat")
# use
library(nat)

# plot some test data
# Drosophila Kenyon cells processed from raw data at http://flycircuit.tw

# help on tihs test data set
?kcs20
# show metadata associated with these neurons
head(kcs20)
# plot all neurons, coloured sequentially
open3d()
plot3d(kcs20)

clear3d()
plot3d(kcs20, col='red')

clear3d()
plot3d(kcs20, col='red', lwd=2)

# get help
?nat
?plot3d.neuronlist
?neuronlist

# Plot first 5 neurons
clear3d()
plot3d(kcs20[1:5])
clear3d()
plot3d(kcs20, subset=type=='gamma', col='red')
plot3d(kcs20, type=='apbp', col='blue')

# Plot all 20 neurons coloured by their Kenyon cell type
clear3d()
with(kcs20, table(type))
plot3d(kcs20, col=type)
