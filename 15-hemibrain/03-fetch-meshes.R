# This script assumed that you have run the file "15-hemibrain/00-setup.R"

## Let's get some neuropil volume data from neuPrint!
### Specifically, in the last R file we chose to look at MBONs
### So let's fetch their meshes for the mushroombody!
### First, we gotta see what is available:
rois = sort(neuprint_ROIs())
rois

# Read in the MB mesh!
mb.mesh = neuprint_ROI_mesh(roi = "MB(R)")
nopen3d(userMatrix = structure(c(0.947127282619476, 0.222770735621452, 
                                 -0.230919197201729, 0, 0.247805655002594, -0.0506941750645638, 
                                 0.967482566833496, 0, 0.203820571303368, -0.973551869392395, 
                                 -0.103217624127865, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.58467960357666, 
        windowRect = c(1460L, 65L, 2877L, 996L)) # set view
plot3d(mb.mesh, add = TRUE, alpha = 0.1, col = lacroix[["pink"]])

# And highlight compartment alpha prime 3
ap3.mesh = neuprint_ROI_mesh(roi = "a'3(R)")
plot3d(ap3.mesh, add = TRUE, alpha = 0.5, col = lacroix[["purple"]])
rgl.snapshot(filename = "images/hemibrain_MB_ap3.png", fmt ="png")

## Maybe get the whole hemibrai mesh?
## hemibrain = neuprint_ROI_mesh(roi = "hemibrain")

# And with some neurons!
nopen3d(userMatrix = structure(c(0.98370349407196, -0.104569993913174, 
                                 -0.146263062953949, 0, 0.147007316350937, -0.000597098842263222, 
                                 0.989135324954987, 0, -0.103521309792995, -0.994517505168915, 
                                 0.0147849693894386, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.644609212875366, 
        windowRect = c(20L, 65L, 1191L, 843L))
plot3d(mb.mesh, add = TRUE, alpha = 0.1, col = "grey")
plot3d(mbons, col = sample(lacroix,length(mbons), replace = TRUE), lwd = 2)
rgl.snapshot(filename = "images/hemibrain_mbons_mb.png", fmt ="png")
