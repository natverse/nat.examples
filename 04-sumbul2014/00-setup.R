if(!require('nat')) install.packages("nat")
if(!require("R.matlab")) install.packages("R.matlab")

library(nat)
library(R.matlab)

# set working directory to location of this file
setwd(dirname(attr(body(function() {}),'srcfile')$filename))