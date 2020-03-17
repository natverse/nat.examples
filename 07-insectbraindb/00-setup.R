# install
if(!require('devtools')) install.packages("devtools")
if(!require('tidyverse')) install.packages("tidyverse")
if(!require('natverse')) devtools::install_github("natverse/natverse")
if(!require('Rvcg')) install.packages("Rvcg")
if(!require('ggplot2')) install.packages("ggplot2")

# load natversee functions
library(natverse)
library(insectbrainr)
library(Rvcg)
library(ggplot2)
library(ggpubr)

# set working directory to location of this file
try(setwd(dirname(attr(body(function() {}),'srcfile')$filename)))

# load data from previoous sessions
if (exists(".RData")){
    load(".RData")
}

