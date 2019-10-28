### This script assumes you have run "/13-bridging/00-setup.R"

# Here, the aim is to show examples of PD2a1/b1 neurons (Dolan et al. 2018, Neuron) from different datasets, even from different labs! (they used to be called 11A in house)

# Let's get these different datasets
pd2a1b1.segs = subset(lh.splits.dps,old.cell.type=="11A")

pd2a1b1.mcfo = subset(lh.mcfo,old.cell.type=="11A")

pd2a1b1.fc = subset(lhns::most.lhns,grepl("PD2",cell.type)&skeleton.type=="FlyCircuit")

pd2a1b1.df = subset(lhns::most.lhns,grepl("PD2",cell.type)&skeleton.type=="DyeFill")

pd2a1b1.jfk = subset(jfw.lhns,grepl("PD2",cell.type))

all.pd2ab.cands = c(em.pd2,pd2a1b1.mcfo,pd2a1b1.fc,pd2a1b1.df,pd2a1b1.jfk)
## These are from the package lhns, and are actually all already in the same brainspace for convenience, FCWB.
## In the example below, we'll put some of them back to where they came from, to show some different spaces.

## First, segmented confocal data. These were original read into R from .nrrd files. These are not single neuron data. Ref Dolan et al. 2019, eLife.
pd2a1b1.segs.jfrc2 = xform_brain(pd2a1b1.segs,sample = FCWB,reference=JFRC2)
nopen3d(userMatrix = structure(c(0.905018091201782, 0.231813445687294, 
                                 -0.356658041477203, 0, -0.0334783419966698, -0.797041535377502, 
                                 -0.602996289730072, 0, -0.424053758382797, 0.557662665843964, 
                                 -0.713576078414917, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, 
        windowRect = c(12L, 44L, 1398L, 904L))
plot3d(JFRC2,alpha=0.1,col="magenta")
plot3d(pd2a1b1.segs.jfrc2[1],col="chartreuse2",lwd=3)
rgl.snapshot(filename = "images/PD2a1b1_confocal_segmentation.png", fmt = "png")

## Next, single cell data from ulti-colour flpout on the same lines (Nern et al. 2015). 
## These were original read into R from .SWC files, generated in Amira. These are not single neuron data. Ref Dolan et al. 2019, eLife.
c11a.mcfo = subset(lh.mcfo,old.cell.type=="11A")
c11a.mcfo.jfrc2 = xform_brain(c11a.mcfo,sample = FCWB,reference=JFRC2)
pd2a1.mcfo.exemplar = c11a.mcfo.jfrc2["L989#6"]
pd2b1.mcfo.exemplar =  c11a.mcfo.jfrc2["L989#2"]
nopen3d(userMatrix = structure(c(0.905018091201782, 0.231813445687294, 
                                 -0.356658041477203, 0, -0.0334783419966698, -0.797041535377502, 
                                 -0.602996289730072, 0, -0.424053758382797, 0.557662665843964, 
                                 -0.713576078414917, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, 
        windowRect = c(12L, 44L, 1398L, 904L))
plot3d(JFRC2,alpha=0.1,col="magenta")
plot3d(pd2a1.mcfo.exemplar,col="chartreuse4",lwd=3,soma=F)
plot3d(pd2b1.mcfo.exemplar,col="chartreuse",lwd=3, soma=F)
rgl.snapshot(filename = "images/PD2a1b1_confocal_MCFO.png", fmt = "png")

# JFRC2
nopen3d(userMatrix = structure(c(0.905018091201782, 0.231813445687294, 
                                 -0.356658041477203, 0, -0.0334783419966698, -0.797041535377502, 
                                 -0.602996289730072, 0, -0.424053758382797, 0.557662665843964, 
                                 -0.713576078414917, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, 
        windowRect = c(12L, 44L, 1398L, 904L))
plot3d(JFRC2,alpha=0.1,col="magenta")
darjeeling = colorRampPalette(wesanderson::wes_palette("Darjeeling1"))
pd2.light = c(pd2a1b1.mcfo,pd2a1b1.fc,pd2a1b1.jfk)
pd2.light[,"col"] = darjeeling(length(pd2.light))
pd2.light = xform_brain(pd2.light,sample=FCWB,reference =JFRC2)
plot3d(pd2.light,lwd=3,soma=T, col=pd2.light[,"col"])
rgl.snapshot(filename = "images/JFRC2_PD2_Light_Library.png", fmt = "png")

# Also the separate data sets in FCWB
nopen3d(userMatrix = structure(c(0.905018091201782, 0.231813445687294, 
                                 -0.356658041477203, 0, -0.0334783419966698, -0.797041535377502, 
                                 -0.602996289730072, 0, -0.424053758382797, 0.557662665843964, 
                                 -0.713576078414917, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, 
        windowRect = c(12L, 44L, 1398L, 904L))
plot3d(FCWB,alpha=0.1,col="grey")
plot3d(pd2a1b1.fc,lwd=3,soma=T, col=darjeeling(length(pd2a1b1.fc)))
rgl.snapshot(filename = "images/FCWB_PD2_FlyCircuit.png", fmt = "png")

# Also the separate data sets in IS2
nopen3d(userMatrix = structure(c(0.905018091201782, 0.231813445687294, 
                                 -0.356658041477203, 0, -0.0334783419966698, -0.797041535377502, 
                                 -0.602996289730072, 0, -0.424053758382797, 0.557662665843964, 
                                 -0.713576078414917, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.7, 
        windowRect = c(12L, 44L, 1398L, 904L))
plot3d(IS2,alpha=0.1,col="grey")
pd2a1b1.df.is2 = xform_brain(pd2a1b1.df, sample = FCWB, reference = IS2)
plot3d(pd2a1b1.df.is2,lwd=3,soma=T, col=darjeeling(length(pd2a1b1.df.is2)))
rgl.snapshot(filename = "images/IS2_PD2_DyeFills.png", fmt = "png")

# Also the separate data sets in JFRC2013
nopen3d(userMatrix = structure(c(0.905018091201782, 0.231813445687294, 
                                 -0.356658041477203, 0, -0.0334783419966698, -0.797041535377502, 
                                 -0.602996289730072, 0, -0.424053758382797, 0.557662665843964, 
                                 -0.713576078414917, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, 
        windowRect = c(12L, 44L, 1398L, 904L))
plot3d(JFRC2013,alpha=0.1,col="grey")
pd2a1b1.mcfo.jfrc2013 = xform_brain(pd2a1b1.mcfo,sample=FCWB,reference = JFRC2013)
plot3d(pd2a1b1.mcfo.jfrc2013,lwd=3,soma=T, col=darjeeling(length(pd2a1b1.mcfo.jfrc2013)))
rgl.snapshot(filename = "images/JFRC2013_PD2_MCFO.png", fmt = "png")

# Also the separate data sets in JFRC2013
clear3d();
plot3d(JFRC2,alpha=0.1,col="magenta")
plot3d(pd2a1b1.jfk,lwd=3,soma=FALSE, col=darjeeling(length(pd2a1b1.jfk)))
rgl.snapshot(filename = "images/JFRC2_PD2_JFW.png", fmt = "png")
# So far we have dealt with light-level data. We can also lok at synaptic, EM Data for these same neurons.
## Connect to the public FAFB instance (Zheng et al. 2018) hosted publicly by Virtual Fly Brain
adult.conn = catmaid_login(server="https://catmaid-fafb.virtualflybrain.org/")

# Let's see some Dolan et al. neurons in their native EM space
nopen3d(userMatrix = structure(c(0.948390364646912, 0.192434668540955, 
                                 -0.252041101455688, 0, 0.172239527106285, -0.979956924915314, 
                                 -0.100092709064484, 0, -0.266250550746918, 0.0515153482556343, 
                                 -0.962526321411133, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, 
        windowRect = c(1440L, 45L, 2826L, 905L))
plot3d(FAFB14,alpha=0.1,col="grey")
rgl.snapshot(filename = "/images/FAFB14_Empty.png", fmt = "png")
emlhns.fafb14 = read.neurons.catmaid("annotation:Dolan et al.")
zissou = colorRampPalette(wesanderson::wes_palette("Zissou1"))
emlhns.fafb14[,"col"] = zissou(length(emlhns.fafb14))
plot3d(emlhns.fafb14,lwd=3,soma=T, col=emlhns.fafb14[,"col"])
rgl.snapshot(filename = "images/FAFB14_EMLHNs.png", fmt = "png")

## Now where to we want to put this data? These LM brains:
nopen3d(userMatrix = structure(c(0.905018091201782, 0.231813445687294, 
                                 -0.356658041477203, 0, -0.0334783419966698, -0.797041535377502, 
                                 -0.602996289730072, 0, -0.424053758382797, 0.557662665843964, 
                                 -0.713576078414917, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, 
        windowRect = c(12L, 44L, 1398L, 904L))
plot3d(JFRC2,alpha=0.1,col="magenta")
rgl.snapshot(filename = "images/JFRC2_Empty.png", fmt = "png")
clear3d();
plot3d(JFRC2013,alpha=0.1,col="magenta")
rgl.snapshot(filename = "images/JFRC2-13_Empty.png", fmt = "png")

## So now PD2a1/b1 neurons transformed into JFRC2 tewmplate space (Dolan et al. 2018)
em.pd2 = fetchn_fafb("name:^PD2a1/b1", mirror = F, reference = JFRC2) # Already bridges to JFRC2 straight from a CATMAID instance!

# There are a few neurons here. For easy viewing, let's choose out favourites.
pd2a1.EM.exemplar = subset(em.11a,grepl("PD2a1#1", name))
pd2b1.EM.exemplar =  subset(em.11a,grepl("PD2b1#1", name))

# And go in for the plotting
nopen3d(userMatrix = structure(c(0.905018091201782, 0.231813445687294, 
                                 -0.356658041477203, 0, -0.0334783419966698, -0.797041535377502, 
                                 -0.602996289730072, 0, -0.424053758382797, 0.557662665843964, 
                                 -0.713576078414917, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, 
        windowRect = c(12L, 44L, 1398L, 904L))
plot3d(JFRC2,alpha=0.1,col="magenta")
plot3d(pd2b1.EM.exemplar,col="grey50",lwd=3,soma=F,WithConnectors = FALSE)
plot3d(pd2a1.EM.exemplar,col="black",lwd=3, soma=F,WithConnectors = FALSE)
rgl.snapshot(filename = "images/PD2a1b1_EM.png", fmt = "png")

# But let's also see the synapses - that's the point of EM data after all!!
clear3d()
plot3d(JFRC2,alpha=0,col="magenta")
plot3d(pd2a1.EM.exemplar,col="chartreuse4",lwd=3, soma=F,WithConnectors = FALSE)
syn.in = get.synapses(pd2a1.EM.exemplar,"POST");spheres3d(xyzmatrix(syn.in),radius = 0.5,col="cyan")
syn.out = get.synapses(pd2a1.EM.exemplar,"PRE");spheres3d(xyzmatrix(syn.out),radius = 0.75,col="red")
rgl.snapshot(filename = "images/PD2a1_EM.png", fmt = "png")

# And again, for PD2b1
clear3d()
plot3d(JFRC2,alpha=0,col="magenta")
plot3d(pd2b1.EM.exemplar,col="chartreuse2",lwd=3, soma=F,WithConnectors = FALSE)
syn.in = get.synapses(pd2b1.EM.exemplar,"POST");spheres3d(xyzmatrix(syn.in),radius = 0.5,col="cyan")
syn.out = get.synapses(pd2b1.EM.exemplar,"PRE");spheres3d(xyzmatrix(syn.out),radius = 0.75,col="red")
rgl.snapshot(filename = "images/PD2b1_EM.png", fmt = "png")

# So now we can show all of the neurons together, ina single space!
# Fix colours for anatomy groups
darjeeling = colorRampPalette(wesanderson::wes_palette("Darjeeling1"))
lhns::most.lhns[,"col"] = darjeeling(length(unique(lhns::most.lhns[,"anatomy.group"])))[as.factor(lhns::most.lhns[,"anatomy.group"])]
most.real.lhns = subset(lhns::most.lhns,good.trace==TRUE&cell.type!="notLHproper")
most.real.lhns = xform_brain(most.real.lhns,sample=FCWB,reference=JFRC2)
nopen3d(userMatrix = structure(c(0.998623549938202, -0.0109847662970424, 
                                 -0.0512875616550446, 0, -0.0147177278995514, -0.9972243309021, 
                                 -0.0729851648211479, 0, -0.050343431532383, 0.0736394673585892, 
                                 -0.996013343334198, 0, 3.54936680187797, -2.4644332365188, 9.03038116879839, 
                                 1), .Dim = c(4L, 4L)), zoom = 0.481017380952835, windowRect = c(1468L, 
                                                                                                 45L, 2835L, 750L))
plot3d(JFRC2,alpha=0.1,col="magenta")
plot3d(most.real.lhns,col = most.real.lhns[,"col"],soma=T,lwd=3)
rgl.snapshot(paste0("/images/All_Library_LHNs.png"),fmt="png")



