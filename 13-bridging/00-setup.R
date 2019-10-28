# install
if(!require('devtools')) install.packages("devtools")
if(!require('devtools')) install.packages("wesanderson")
if(!require('tidyverse')) install.packages("tidyverse")
if(!require('natverse')) devtools::install_github("natverse/natverse")
if(!require('lhns')) devtools::install_github("jefferislab/lhns")
if(!require('catnat')) devtools::install_github("jefferislab/catnat")

# load natversee functions
library(natverse)
library(tidyverse)
library(lhns)
library(catnat)
library(wesanderson) # For nice, cinema-inspired colours
## What is lhns? It is a little data package for lateral horn neurons, made for Frechter et al. 2019, eLife.
## It contains some data we will play with in the next R file
## What is catnat? It is an experimental package for natverse functions.

# set working directory to location of this file
try(setwd(dirname(attr(body(function() {}),'srcfile')$filename)))

# load data from previous sessions
if (exists(".RData")){
  load(".RData")
}


