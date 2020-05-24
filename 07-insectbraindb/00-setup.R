# install
if(!require("natmanager")) install.packages("natmanager")
if(!require('natverse')) natmanager::install('natverse')
if(!require('Rvcg')) install.packages("Rvcg")
if(!require('ggplot2')) install.packages("ggplot2")

# load natversee functions
library(natverse)
library(insectbrainr)
library(Rvcg)
library(ggplot2)
library(ggpubr)

# set working directory to location of this file
setwd(here::here('07-insectbraindb'))

