# set the scene
## This script assumed that you have run the file "08-openworm/00-setup.R"

### A folder of NeuroML files was downloaded from a openworm (http://openworm.org/) Github repository (https://github.com/openworm/CElegansNeuroML) and included in this repository
# get data
celegans.neurons = read.neurons("data/generatedNeuroML/")

# If you want to read specific neurons straight from Github, you can do this as so:
vds=paste0("VD", 1:13)
vdurls=paste0("https://raw.githubusercontent.com/openworm/CElegansNeuroML/",
              "103d500e066125688aa7ac5eac7e9b2bb4490561/CElegans/generatedNeuroML/",vds,
              ".morph.xml")
vdnl=read.neurons(vdurls, neuronnames=vds)
plot3d(vdnl)
