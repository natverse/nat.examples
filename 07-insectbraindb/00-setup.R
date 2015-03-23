if(!require('nat')) install.packages("nat")
if(!require('rvest')) install.packages("rvest")

library(nat)
library(rvest)

# set working directory to location of this file
setwd(dirname(attr(body(function() {}),'srcfile')$filename))