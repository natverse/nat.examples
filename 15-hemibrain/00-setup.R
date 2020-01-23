# install
if(!require('devtools')) install.packages("devtools")
if(!require('natverse')) devtools::install_github("natverse/natverse")

# load
library(natverse)
library(neuprintr)
library(dendroextras)

# set working directory to location of this file
try(setwd(dirname(attr(body(function() {}),'srcfile')$filename)))
