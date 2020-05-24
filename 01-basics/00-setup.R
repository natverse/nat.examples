# install
if(!require("natmanager")) install.packages("natmanager")
# See http://natverse.org/install/ for details / troubleshooting

if(!require('natverse')) natmanager::install("natverse")

# This is a big download. 
if(!require('lhns')) devtools::install_github("jefferislab/lhns")
# If you have trouble installing lhns with `install_github``

## in the terminal. Choose a directory that you like
# cd  ~/natverse
## use --depth to download only the latest revision (1GB instead of 2B for whole history)
# git clone --depth 1 https://github.com/jefferislab/lhns

## This still took me serveral mins to download, but it is much less fragile that the install_github mechanism and you get some indication of progress

## Then, in R (adjusting path to git checkout of lhns as appropriate)
# remotes::install_local("~/natverse/lhns")


# load natverse functions
library(natverse)
library(tidyverse)
library(lhns)
## What is lhns? It is a little data package for lateral horn neurons, made for Frechter et al. 2019, eLife.
## It contains some data we will play with in the next R file

# set working directory to location of this file
setwd(here::here('01-basics'))
