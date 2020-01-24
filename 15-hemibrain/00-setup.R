# install
if(!require('devtools')) install.packages("devtools")
if(!require('natverse')) devtools::install_github("natverse/natverse")
if(!require('ComplexHeatmap')) install.packages("ComplexHeatmap")
if(!require('devtools')) install.packages("ggnetwork")

# load
library(natverse)
library(neuprintr)
library(dendroextras)
library(ComplexHeatmap)
library(ggnetwork)

# set working directory to location of this file
try(setwd(dirname(attr(body(function() {}),'srcfile')$filename)))

# Colours
# some nice colors!! Inspired by LaCroixColoR
lacroix = c("#C70E7B", "#FC6882", "#007BC3", "#54BCD1", "#EF7C12", "#F4B95A", 
            "#009F3F", "#8FDA04", "#AF6125", "#F4E3C7", "#B25D91", "#EFC7E6", 
            "#EF7C12", "#F4B95A", "#C23A4B", "#FBBB48", "#EFEF46", "#31D64D", 
            "#132157","#EE4244", "#D72000", "#1BB6AF")
names(lacroix) = c("purple", "pink",
                   "blue", "cyan",
                   "darkorange", "paleorange",
                   "darkgreen", "green",
                   "brown", "palebrown",
                   "mauve", "lightpink",
                   "orange", "midorange",
                   "darkred", "darkyellow",
                   "yellow", "palegreen", 
                   "navy","cerise",
                   "red", "marine")
