# set the scene
## This script assumed that you have run the file "11-frechter2019/00-setup.R"

# get data
# Show neurons
opns = subset(lhns::most.lhins, grepl("Olf", modality))
nopen3d(userMatrix = structure(c(0.905018091201782, 0.231813445687294, 
                                 -0.356658041477203, 0, -0.0334783419966698, -0.797041535377502, 
                                 -0.602996289730072, 0, -0.424053758382797, 0.557662665843964, 
                                 -0.713576078414917, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, 
        windowRect = c(12L, 44L, 1398L, 904L))
plot3d(FCWB,alpha=0.1,col="lightgrey")
plot3d(opns, lwd = 2, soma = TRUE)
rgl.snapshot(filename ="images/nat_adultfly_neurons.png" ,fmt = "png")

# Though, this is what it looked like before we mirrored it
nopen3d(userMatrix = structure(c(0.905018091201782, 0.231813445687294, 
                                 -0.356658041477203, 0, -0.0334783419966698, -0.797041535377502, 
                                 -0.602996289730072, 0, -0.424053758382797, 0.557662665843964, 
                                 -0.713576078414917, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.556837737560272, 
        windowRect = c(12L, 44L, 1398L, 904L))
plot3d(FCWB,alpha=0.1,col="lightgrey")
plot3d(subset(opns, LH_side == "R"), lwd = 2, soma = TRUE)
plot3d(mirror_brain(subset(opns, LH_side == "L"), brain = FCWB), lwd = 2, soma = TRUE)
rgl.snapshot(filename ="images/nat_adultfly_neurons_leftright.png" ,fmt = "png")

## Examine a right side and a left side neuron
da1.r = subset(opns, glomerulus == "DA1" & LH_side =="R")
da1.l = subset(opns, glomerulus == "DA1" & LH_side =="L")
da1.lm = mirror_brain(x = da1.l, brain = FCWB)
plot3d(da1.r, lwd = 2, soma = TRUE, col = "red")
plot3d(da1.lm, lwd = 2, soma = TRUE, col = "cyan")
plot3d(FCWB,alpha=0.1,col="lightgrey")
rgl.snapshot(filename ="images/nat_adultfly_da1_leftright.png" ,fmt = "png")
clear3d()
plot3d(da1.r, lwd = 2, soma = TRUE, col = "red")
plot3d(da1.l, lwd = 2, soma = TRUE, col = "cyan")
plot3d(FCWB,alpha=0.1,col="lightgrey")
rgl.snapshot(filename ="images/nat_adultfly_da1_mirrored.png" ,fmt = "png")


## We can also plot every one of the Frechter at al. 2019 lateral horn cell types
### With its name
### And a little model of the lateral horn (LH)
### Let's try that
dir.create("images/LHNs/")
lhons = subset(lhns::most.lhns,type=="ON")
lhlns = subset(lhns::most.lhns,type=="LN")
lhons.cols = colorRampPalette(c("darkblue","deepskyblue1"))(length(unique(lhons[,"cell.type"])))
lhlns.cols = colorRampPalette(c("darkgreen","chartreuse1"))(length(unique(lhlns[,"cell.type"])))
cols = c(lhons.cols,lhlns.cols)
names(cols) = c(sort(unique(lhons[,"cell.type"])),sort(unique(lhlns[,"cell.type"])))
letter.location = c(150, 34, 78.6946086)
letter.location2 = c(125, 29, 78.6946086)
nopen3d(userMatrix = structure(c(0.988249838352203, 0.00959446653723717, 
                                 -0.152546763420105, 0, -0.0010739240096882, -0.997567236423492, 
                                 -0.0696997791528702, 0, -0.152844399213791, 0.0690445750951767, 
                                 -0.985835373401642, 0, 62.5438867677425, -57.9483397676875, 3.67746390670118, 
                                 1), .Dim = c(4L, 4L)), zoom = 0.376889795064926, windowRect = c(1440L, 
                                                                                                 45L, 2704L, 971L))
for(ct in sort(unique(lhns::most.lhns[,"cell.type"]),decreasing=FALSE)){
  plot3d(FCWB,alpha=0)
  plot3d(subset(FCWBNP.surf,"LH_R"),alpha = 0.1,col="grey")
  text3d(letter.location2,texts = ct,col="black",font=2,cex=5)
  plot3d(subset(lhns::most.lhns,cell.type==ct)[1],lwd=4,soma=T,col=cols[ct])
  rgl.snapshot(paste0("images/LHNs/LHN_",ct,".png"),fmt="png")
  clear3d()
}

