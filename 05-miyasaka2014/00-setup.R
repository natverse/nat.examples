# install
if(!require('devtools')) install.packages("devtools")
if(!require('tidyverse')) install.packages("tidyverse")
if(!require('natverse')) devtools::install_github("natverse/natverse")

# load natversee functions
library(natverse)
library(tidyverse)

# set working directory to location of this file
try(setwd(dirname(attr(body(function() {}),'srcfile')$filename)))
