# install
if(!require('devtools')) install.packages("devtools")
if(!require('tidyverse')) install.packages("tidyverse")
if(!require('natverse')) devtools::install_github("natverse/natverse")
if(!require('lhns')) devtools::install_github("jefferislab/lhns")

# load natversee functions
library(natverse)
library(tidyverse)
library(lhns)
## What is lhns? It is a little data package for lateral horn neurons, made for Frechter et al. 2019, eLife.
## It contains some data we will play with in the next R file

# set working directory to location of this file
try(setwd(dirname(attr(body(function() {}),'srcfile')$filename)))


