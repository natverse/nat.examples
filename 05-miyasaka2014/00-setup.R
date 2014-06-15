if(!require('nat')) install.packages("nat")
if(!require('devtools')) install.packages("devtools")
if(!require('nat.nblast')) devtools::install_github('jefferislab/nat.nblast')

library(nat)

# set working directory to location of this file
setwd(dirname(attr(body(function() {}),'srcfile')$filename))
