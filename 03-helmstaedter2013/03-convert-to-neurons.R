# assumes you have run 00-setup.R

load("sk.uniq.rda")

library(nat)
skn=nlapply(sk.uniq, as.neuron, OmitFailures = T, .progress='text')
# save neurons in native R format
save(skn, file='skn.rda')
# save a zip archive of SWC format neurons for all reconstructions
write.neurons(skn, dir='skuniq.swc.zip', files=names(skn), format='swc')
