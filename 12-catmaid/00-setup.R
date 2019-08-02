# install
if(!require('devtools')) install.packages("devtools")
if(!require('tidyverse')) install.packages("tidyverse")
if(!require('natverse')) devtools::install_github("natverse/natverse")

# load natversee functions
library(natverse)
library(tidyverse)

# set working directory to location of this file
try(setwd(dirname(attr(body(function() {}),'srcfile')$filename)))

# load data from previoous sessions
if (exists(".RData")){
  load(".RData")
}

# Extra functions

process_lhn_name <- function(x) {
  res=stringr::str_match(x, "([AP][DV][1-9][0-9]{0,1})([a-z])([1-9][0-9]{0,2})")
  df=data.frame(pnt=res[,2], anatomy.group=paste0(res[,2], res[,3]), cell.type=res[,1],stringsAsFactors = F)
  isna=is.na(res[,2]) | is.na(res[,3])
  df[['anatomy.group']][isna]=NA
  df
}
