# install
install.packages("nat")
# use
library(nat)

# plot some test data (?kcs20 for details)
# Drosophila Kenyon cells processed from raw data at http://flycircuit.tw
head(kcs20)
plot3d(kcs20, col=type)

# get help
?nat
