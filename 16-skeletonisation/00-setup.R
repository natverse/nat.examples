# install
if(!require('remotes')) install.packages("remotes")
if(!require('natverse')) remotes::install_github("natverse/natverse")
if(!require('fafbseg')) remotes::install_github("natverse/fafbseg")
if(!require('alphashape3d')) install.packages("alphashape3d")
if(!require('reticulate')) install.packages("reticulate")
if(!require('Rvcg')) install.packages("Rvcg")

# load
library(natverse)
library(fafbseg)
library(alphashape3d)
library(Rvcg)

# In order to work with the skeletonisation function, skeletor
## You need to show R how to play with python and P. Schlegel's skeletor library
### Detailed instructions can be found in a vignette here: http://natverse.org/fafbseg/articles/articles/installing-cloudvolume-meshparty.html
### You need to have blender for this to work  https://www.blender.org/download/
#### The next .R file will not run without it!

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

# Get a correctly orientated view of the fibre
fibre_view <- function(){
  nopen3d(userMatrix = structure(c(-0.0326365306973457, 0.0871003195643425, 
                                   -0.995664954185486, 0, -0.985228538513184, -0.170360401272774, 
                                   0.0173917300999165, 0, -0.168106958270073, 0.981525063514709, 
                                   0.0913735032081604, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.584679365158081, 
          windowRect = c(1440L, 101L, 3155L, 1027L))
}
