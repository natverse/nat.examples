# assumes you have run 00-setup.R

load("sk.uniq.rda")

library(nat)
skn=nlapply(sk.uniq, as.neuron, OmitFailures = T, .progress='text')
save(skn, file='skn.rda')
