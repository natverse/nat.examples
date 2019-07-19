if(!require('nat')) install.packages("nat")
if(!require('rvest')) install.packages("rvest")
if(!require('nat.nblast')) install.packages("nat.nblast")
if(!require('dendextend')) install.packages("dendextend")
if(!require('devtools')) install.packages("devtools")
if(!require('neuromorphr')) devtools::install_github("jefferislab/neuromorphr")

library(nat)
library(nat.nblast)
library(rvest)
library(neuromorphr)

# set working directory to location of this file
try(setwd(dirname(attr(body(function() {}),'srcfile')$filename)))
