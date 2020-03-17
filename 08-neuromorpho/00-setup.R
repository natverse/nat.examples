# install
if(!require('devtools')) install.packages("devtools")
if(!require('tidyverse')) install.packages("tidyverse")
if(!require('natverse')) devtools::install_github("natverse/natverse")
if(!require('rvest')) install.packages("rvest")
if(!require('ggplot2')) install.packages("ggplot2")
if(!require('Rtsne')) install.packages("Rtsne")

# load natversee functions
library(natverse)
library(ggplot2)
library(Rtsne)

# set working directory to location of this file
try(setwd(dirname(attr(body(function() {}),'srcfile')$filename)))
