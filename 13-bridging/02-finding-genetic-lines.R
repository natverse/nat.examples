### This script assumes you have run "/13-bridging/00-setup.R"

# Here, the aim is to take an EM neuron an find a match for it in terms of locating a genetic resource.

# First, let's say we want to find a new EM neurons in the Mushroom body - the memory centre - of the fly
## We can plot the mushroom body, and jump to its location in a live EM CATMAID instance!!
## First, let's connect to the public FAFB instance (Zheng et al. 2018) hosted publicly by Virtual Fly Brain
adult.conn = catmaid_login(server="https://catmaid-fafb.virtualflybrain.org/")
nopen3d(userMatrix = structure(c(0.990703582763672, 0.0126308705657721, 
                                 -0.135451108217239, 0, -0.00503587257117033, -0.9915931224823, 
                                 -0.129299566149712, 0, -0.135945424437523, 0.128779590129852, 
                                 -0.982310831546783, 0, 67.1707067927694, -81.6641592721123, 10.5346586315227, 
                                 1), .Dim = c(4L, 4L)), zoom = 0.242946609854698, windowRect = c(1440L, 
                                                                                                 45L, 2826L, 905L))
plot3d(JFRC2,alpha=0,col="magenta")
plot3d(subset(JFRC2NP.surf,c("MB_VL_R","MB_ML_R","MB_PED_R")),alpha=0.1,col="purple")
rgl.snapshot(filename = "images/JFRC2_MB_Lobes.png", fmt = "png")

# Jump!!
open_fafb(subset(JFRC2NP.surf,c("MB_VL_R","MB_ML_R","MB_PED_R")))

# Now, let's say we have traced a new neuron, in this case a mushroom body output neuron - an MBON.
## We can read it into r and bridge it to the brainspace that contains the most genetic line, JFRC2
### NOTE: This example will not work unless you have access to a restricted CATMAID instance
### HOWEVER you can follow this pipleine with any other neuron you might find in the public CATMAID instance
FAFBv14.conn = catmaid_login(server = "https://neuropil.janelia.org/tracing/fafb/v14/", 
                             authname = "YOUR_USER_HERE", authpassword = "YOUR_PASSWORD_HERE", 
                             token = "YOUR_TOKEN_HERE") # Might not work for you
em.mbonax = fetchn_fafb("annotation:^MBON aX Right ASB$",mirror=FALSE,reference = JFRC2)
em.mbonax = nlapply(em.mbonax,unspike,threshold=50)
em.mbonax = nlapply(em.mbonax,smooth_neuron,method="gauss",sigma=2)

# Plot the new MBON
## Front on
nopen3d(userMatrix = structure(c(0.990703582763672, 0.0126308705657721, 
                                 -0.135451108217239, 0, -0.00503587257117033, -0.9915931224823, 
                                 -0.129299566149712, 0, -0.135945424437523, 0.128779590129852, 
                                 -0.982310831546783, 0, 67.1707067927694, -81.6641592721123, 10.5346586315227, 
                                 1), .Dim = c(4L, 4L)), zoom = 0.242946609854698, windowRect = c(1440L, 
                                                                                                 45L, 2826L, 905L))
plot3d(JFRC2,alpha=0,col="magenta")
plot3d(subset(JFRC2NP.surf,c("MB_VL_R","MB_ML_R","MB_PED_R")),alpha=0.1,col="purple")
plot3d(em.mbonax,col="black",lwd=3,soma=T,WithConnectors=TRUE)
rgl.snapshot(filename = "images/JFRC2_MBONax.png", fmt = "png")
## Side on
nopen3d(userMatrix = structure(c(-0.206468492746353, 0.0996669754385948, 
                                 -0.973364174365997, 0, -0.117232859134674, -0.990152478218079, 
                                 -0.0765188336372375, 0, -0.971404790878296, 0.0983114689588547, 
                                 0.216119319200516, 0, -47.731201171875, -82.9554214477539, -77.7535095214844, 
                                 1), .Dim = c(4L, 4L)), zoom = 0.242946609854698, windowRect = c(1440L, 
                                                                                                 45L, 2826L, 905L))
plot3d(JFRC2,alpha=0,col="magenta")
plot3d(subset(JFRC2NP.surf,c("MB_VL_R","MB_ML_R","MB_PED_R")),alpha=0.1,col="purple")
plot3d(em.mbonax,col="black",lwd=3,soma=T,WithConnectors=TRUE)
rgl.snapshot(filename = "images/JFRC2_MBONax_side.png", fmt = "png")

# Now let's attempt to find a GAL4 GMR line for it
em.mbonax.search = fetchn_fafb("annotation:^MBON aX Right ASB$",mirror=FALSE,reference = FCWB)
gmrdps=readRDS("data/gmrdps.rds")
doMC::registerDoMC(7)
x=nblast_fafb("2062703", reference=FCWB, mirror=F, db=gmrdps, UseAlpha = T)
summary(x, db=gmrdps) # the results are sorted according to `muscore`, which is the mean forward and reverse score.

# Scan through the hits
nopen3d()
for(i in 1:100){
  clear3d()
  message(i)
  plot3d(x, db=gmrdps, hits=i)
  progress = readline("Next? ")
}
# Maybe 75A02, 12G09, 34E09,91F03

# Let's plot our good matches
nopen3d(userMatrix = structure(c(0.998708069324493, -0.0172212868928909, 
                                 -0.0478079840540886, 0, -0.018259335309267, -0.999605298042297, 
                                 -0.0213626269251108, 0, -0.0474211499094963, 0.0222078841179609, 
                                 -0.998628079891205, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, 
        windowRect = c(1469L, 47L, 2855L, 907L))
plot3d(JFRC2,col="magenta",alpha=0.1)
plot3d(em.mbonax,col="black",lwd=3,soma=T,WithConnectors = TRUE)
gmr.matches = xform_brain(gmrdps[c("75A02", "12G09", "34E09","91F03")],sample=FCWB,reference=JFRC2)
plot3d(gmr.matches,lwd=3)
rgl.snapshot(filename = "images/MBONax_EM_GMRmatch.png", fmt = "png")
