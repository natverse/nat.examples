### This script assumes you have run "/13-bridging/00-setup.R"

# Here, the aim is to take an EM neuron an find a match for it in terms of locating a genetic resource.

# First, let's say we want to find a new EM neurons in the Mushroom body - the memory centre - of the fly
### NOTE: This example will not work unless you have access to a restricted CATMAID instance
### HOWEVER you can follow this pipleine with any other neuron you might find in the public CATMAID instance
FAFBv14.conn = catmaid_login(server = "https://neuropil.janelia.org/tracing/fafb/v14/", 
                             authname = "YOUR_USER_HERE", authpassword = "YOUR_PASSWORD_HERE", 
                             token = "YOUR_TOKEN_HERE") # Might not work for you

# Find interesting neuron tracing
new.pn = fetchn_fafb(4954525,mirror=FALSE,reference=JFRC2)
new.pn = nlapply(new.pn,smooth_neuron,method="gauss",sigma=2)
new.pn.fcwb = fetchn_fafb(4954525,mirror=FALSE,reference=FCWB)

# NBLAST against GMR GAL4s:
gmrdps=readRDS("data/gmrdps.rds")
x=nblast_fafb(4954525, reference=FCWB, mirror=F, db=gmrdps, UseAlpha = T)
summary(x, db=gmrdps) # the results are sorted according to `muscore`, which is the mean forward and reverse score.

# Scan through the hits
nopen3d()
for(i in 1:100){
  clear3d()
  message(i)
  plot3d(x, db=gmrdps, hits=i)
  progress = readline("Next? ")
}
example.matches = c("75B02", "75B08","28D07", "37F11", "59E06", "28C03", "83A11", "28C06", 
                    "94C04", "52D09") # last one is pretty Clean

# Plot neuron alone

## FAFB14
nopen3d(userMatrix = structure(c(0.948390364646912, 0.192434668540955, 
                                 -0.252041101455688, 0, 0.172239527106285, -0.979956924915314, 
                                 -0.100092709064484, 0, -0.266250550746918, 0.0515153482556343, 
                                 -0.962526321411133, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, 
        windowRect = c(1440L, 45L, 2826L, 905L))
plot3d(FAFB14,alpha=0.1,col="grey")
new.pn.fafb14 = read.neurons.catmaid("4954525")
plot3d(new.pn.fafb14,lwd=3,soma=T, col= "black")
rgl.snapshot(filename = "images/FAFB14_NewPN_EM.png", fmt = "png")

## JFRC2
nopen3d(userMatrix = structure(c(0.905018091201782, 0.231813445687294, 
                                 -0.356658041477203, 0, -0.0334783419966698, -0.797041535377502, 
                                 -0.602996289730072, 0, -0.424053758382797, 0.557662665843964, 
                                 -0.713576078414917, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, 
        windowRect = c(12L, 44L, 1398L, 904L))
plot3d(JFRC2,alpha=0.1,col="magenta")
plot3d(new.pn,lwd=3,soma=T, col= "black")
rgl.snapshot(filename = "images/JFRC2_NewPN_EM.png", fmt = "png")

# Plot GMR match
nopen3d(userMatrix = structure(c(0.998708069324493, -0.0172212868928909, 
                                 -0.0478079840540886, 0, -0.018259335309267, -0.999605298042297, 
                                 -0.0213626269251108, 0, -0.0474211499094963, 0.0222078841179609, 
                                 -0.998628079891205, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, 
        windowRect = c(1469L, 47L, 2855L, 907L))
plot3d(JFRC2,col="magenta",alpha=0.1)
plot3d(new.pn,col="black",lwd=3,soma=T,WithConnectors = FALSE)
gmr.match = xform_brain(gmrdps["75B02"],sample=FCWB,reference=JFRC2)
plot3d(gmr.match,lwd=3,col="chartreuse")

# Plot split match
nopen3d(userMatrix = structure(c(0.998708069324493, -0.0172212868928909, 
                                 -0.0478079840540886, 0, -0.018259335309267, -0.999605298042297, 
                                 -0.0213626269251108, 0, -0.0474211499094963, 0.0222078841179609, 
                                 -0.998628079891205, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, 
        windowRect = c(1469L, 47L, 2855L, 907L))
plot3d(JFRC2,col="magenta",alpha=0.1)
plot3d(new.pn,col="black",lwd=3,soma=T,WithConnectors = FALSE)
split.match = xform_brain(lh.mcfo["L2444#1"],sample=FCWB,reference=JFRC2)
plot3d(split.match,lwd=3,col="chartreuse4",soma=T)
rgl.snapshot(filename = "images/NewPN_EM_SplitGal4match.png", fmt = "png")

# Generate histogram
result = nblast(dotprops(new.pn.fcwb,resample=1),gmrdps,normalised = TRUE,UseAlpha=TRUE)
result = data.frame(result)
p<-ggplot(result, aes(x=X4954525)) +
  geom_histogram(position="identity", alpha=0.5)+
  scale_color_grey()+scale_fill_grey() +
  theme_classic() +
  guides(colour=FALSE, fill = FALSE)
p
rgl.snapshot(filename = "images/NewPN_EM_GMRmatch.png", fmt = "png")
