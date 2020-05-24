# Install the nat tools
if(!require("natmanager")) install.packages("natmanager")
# See http://natverse.org/install/ for details / troubleshooting
natmanager::install("natverse")

if(!require("ggpubr")) install.packages("ggpubr")
if(!require("reshape2")) install.packages("reshape2")
if(!require("catnat")) remotes::install_github("jefferislab/catnat")
