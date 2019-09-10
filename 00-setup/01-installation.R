# Install the nat tools
if(require(devtools)){ install.packages("devtools")}
if(require(igraph)){ install.packages("igraph")}
if(require(Morpho)){ install.packages("Morpho")}
if(require(ggpubr)){install.packages("ggpubr")}
if(require(reshape2)){install.packages("reshape2")}
if(require(natverse)){devtools::install_github("natverse/natverse")}
if(require(catnat)){devtools::install_github("jefferislab/catnat")}

