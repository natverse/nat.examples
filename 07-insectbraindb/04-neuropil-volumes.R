## This script assumed that you have run the file "07-insectbraindb/00-setup.R"

## Righto, which species does the database host?
species_info = insectbraindb_species_info()
insect.species = unique(species_info$scientific_name)

## Okay, so a bunch of bees, wasps, moths, beetles and 'worms'
print(unique(species_info$common_name))

## Let's try to read every brain the repository has
insect.brains = list()
for(insect in insect.species){
  message("Pulling ", insect, " brain mesh")
  for(sex in c("UNKNOWN", "MALE", "FEMALE")){
    message(sex)
    insect.brain = insectbraindb_read_brain(species = insect, brain.sex  = sex, progress = TRUE)
    insect.sex = paste(insect, sex, sep = " ")
    if(!is.null(insect.brain)){
      insect.brains[[insect.sex]] = insect.brain
    }
  }
}
## You will notice that we failed to read two brains there. It seems there is no reconstruction yet for Bombus terrestris, and that the Helicoverpa armigera requires a user with an account on insectbraindb.org 

## Because you donwload .obj files to the same temporary directory, which persists as long as the R session, re-reading the same brain again is quicker - because you have already downloaded the underlying files! Try re-running the above code in your R console to see (Markdown behaviour can be different).

## Let's also chuck in a Drosophila melanogaster brain for good measure 
if(!require('devtools')) install.packages("devtools")
if(!require('nat.flybrains')) devtools::install_github("jefferislab/nat.flybrains")
JFRC2NP.surf = nat.flybrains::JFRC2NP.surf
JFRC2NP.surf$scientific_name = "Drosophila melanogaster"
JFRC2NP.surf$common_name = "Vinegar fly" # Not fruit fly
JFRC2NP.surf$sex = "UNKNOWN" # Actually, it is intersex
insect.brains[["Drosophila melanogaster UNKNOWN"]] = JFRC2NP.surf
## It is the best of the insects after all

### We need to center all of the brains at zero, in order to visually compare them
insect.brains = lapply(insect.brains, function(ib){
  center = colMeans(xyzmatrix(ib))
  xyzmatrix(ib) = xyzmatrix(ib) - matrix(center,ncol=3, nrow = nrow(xyzmatrix(ib)), byrow = TRUE)
  ib
})

## In order to see everything in thei relative sizes, we need to plot the bounding box of the biggest brain
insect.brain.points = do.call(rbind,lapply(insect.brains,nat::xyzmatrix))
bbx = nat::boundingbox(insect.brain.points)
bbxm = as.mesh3d(bbx)

## Hold right click to pan
nat::nopen3d(userMatrix = structure(c(0.998683869838715, -0.00999264605343342, 
                                      -0.050302866846323, 0, -0.0184538997709751, -0.985155284404755, 
                                      -0.170669943094254, 0, -0.0478506684303284, 0.171373650431633, 
                                      -0.984043300151825, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.600000023841858, 
             windowRect = c(1440L, 45L, 3209L, 1063L))
for(ib in insect.brains){
  clear3d()
  plot3d(bbxm, alpha = 0, add = TRUE)
  message(ib$scientific_name, " the ", ib$common_name)
  plot3d(ib)
  progress = readline(prompt = "Press any key for next brain ... ")
  clear3d();    plot3d(bbxm, alpha = 0, add = TRUE); plot3d(ib, alpha = 0.3)
  rgl.snapshot(file = paste0("images/insectbrains/",paste(unlist(strsplit(ib$common_name," ")),collapse="_"), "_",ib$sex, ".png"), fmt = "png")
}
### Good stuff

## Seems like not all of these brains are fully fleshed out. Let's remove the partial reconstructions
insect.brains = insect.brains[!grepl("Agrotis segetum",names(insect.brains))]

## We can also subset these brains by a particular brain area. Everyone likes the antennal lobe, maybe it's the best studied bit of insect neuro-anatomy. So let's have a gander at that.

## Let's have a look at how standardised the neuropil names are, across these species
neuropils = lapply(insect.brains, function(brain) sort(brain$neuropil_full_names))
als = sapply(neuropils, function(np) "Antennal Lobe"%in%np)
als
## Okay, so we have it in all, and indeed it has a capital A and a capital L

## Brains that contain a antennal lobe neuropil object
with.al = sapply(insect.brains, function(ib) sum(grepl("^Antennal Lobe$", ib$neuropil_full_names))>0)
insect.brains.with.al = insect.brains[with.al]
insect.brains.with.al[["Drosophila melanogaster UNKNOWN"]] = JFRC2NP.surf # Got dropped out, as it does not have all the entries of a insectbraindb read brain

## Fab. Now what we really might like to know if how the volume of the antennal lobe differs across species. So let's put our analysis hat on and take a look at that.

## There are a few ways of calculating a volume for a mesh. I am going to be lazy and use the package alphashape3d.
## We will also use pbapply, because volume calculation can take some time, and it is nice to know how things are going
if(!require('alphashape3d')) install.packages("alphashape3d")
if(!require('pbapply')) install.packages("pbapply")

## Create function to calculate volume
calculate_volume <- function(brain, neuropil = NULL, alpha = 10, scale = 1e-9, method = c("convhull","ashape"), resample = 0.1, ...){
  method = match.arg(method)
  if(is.null(neuropil)){
    brain = subset(brain, neuropil)
  }
  if(resample){
    brain = vcgUniformRemesh(as.mesh3d(brain), voxelSize = resample, offset = 0, discretize = FALSE,
                             multiSample = FALSE, absDist = FALSE, mergeClost = FALSE,
                             silent = FALSE)
  }
  points = unique(nat::xyzmatrix(subset(brain, neuropil)))
  if(method=="ashape"){
    a = alphashape3d::ashape3d(points, alpha = alpha, pert = TRUE)
    volume_ashape3d(a)*scale # convert to mm3
  }else{
    geometry::convhulln(points, output.options = "FA")$vol*scale
  }
}

## This was of calculating volumes doesn't seem terribly accurate though, very sensitive to alpha value?
## It would be better to use this R package, however, its installation is a little involved, esp.p on Windows
## It is probably worth the time though, here it is: https://github.com/zarquon42b/RvtkStatismo

## Once installed, we'll need to load it
library(RvtkStatismo)

## Get the volumes for the whole brain
insect.brain.volumes = pbapply::pbsapply(insect.brains.with.al, calculate_volume, neuropil = NULL)

## Get the volumes for just the antennal lobe
insect.al.volumes = pbapply::pbsapply(insect.brains.with.al, calculate_volume, neuropil = "^AL_left|^AL_right|^AL_noside|^AL_R$|^AL_L$")

## Hmm, that took a while

## Assemble data.frame
species.with.al = sapply(insect.brains.with.al, function(x) x$common_name)
sex.with.al = unlist(sapply(insect.brains.with.al, function(x) x$sex[[1]]))
df = data.frame(species = c(species.with.al, species.with.al),
                sex = c(sex.with.al, sex.with.al),
                volume = c(insect.brain.volumes, insect.al.volumes),
                neuropil = c(rep("whole",length(insect.brain.volumes)),rep("AL",length(insect.al.volumes)))
)

## Plot!
pdf("images/insectbraindb_volume_and_AL_volume_comparison.pdf", width = 10, height = 7)
ggplot2::ggplot(df, aes(x=species, y=volume, color = neuropil, group = neuropil, shape=sex)) +
  geom_jitter(position=position_dodge(0.2))+
  theme_classic()
dev.off()

## Hmm, let's see total brain volume against antennal lobe volume
df2 = data.frame(species = species.with.al,
                 sex = sex.with.al,
                 volume = insect.brain.volumes,
                 al.volume = insect.al.volumes
)
pdf("images/insectbraindb_volume_vs_AL_volume.pdf", width = 10, height = 7)
ggplot2::ggplot(df2, aes(x=volume, y=al.volume, color = species)) +
  geom_point() + 
  geom_smooth(data = df2, aes(x=volume, y=al.volume), color = "black", method=lm, se = FALSE)+
  theme_classic()
dev.off()

## Also note, that some of these brains, such as that for Agrotis segetum, are not actually complete seemingly.

## It is hard to interpret this sort of result. Perhaps we are better served by trying to see an olfactory - vision trade-off, and sompare the AL size with the size of the optic lobes? The optic lobes are a bit more complicated

## Which insects have which major optic neuropils?
lamina = sapply(neuropils, function(np) sum(grepl("Lamina", np)))
lamina
lobula = sapply(neuropils, function(np) sum(grepl("Lobula", np)))
lobula
medulla = sapply(neuropils, function(np) sum(grepl("Medulla", np)))
medulla
## Looks like the lamina is missing from many of these brains. We'll combine lobula and medulla volumes for our calculation

## Brains that contain a optic lobe neuropil objects, though not necessarily all of them
with.optic = sapply(insect.brains, function(ib) sum(grepl("^Lobula|^Medulla", ib$neuropil_full_names))>0)
insect.brains.with.optic = insect.brains[with.optic]

## Calculate the neuropil volumes
insect.optic.volumes = pbapply::pbsapply(insect.brains.with.optic, function(x)
  calculate_volume(subset(x, x$RegionList[grepl("^Lobula|^Medulla",x$neuropil_full_names)]))
)

## And our friend, Drosophila melanogaster
insect.brains.with.optic[["Drosophila melanogaster UNKNOWN"]] = JFRC2NP.surf # Got dropped out, as it does not have all the entries of a insectbraindb read brain
insect.optic.volumes[["Drosophila melanogaster UNKNOWN"]] = calculate_volume(subset(JFRC2NP.surf, "LOP_|ME_|LO_"), alpha = 3)

## Assemble data.frame
species.with.optic = sapply(insect.brains.with.optic, function(x) x$common_name)
sex.with.optic = unlist(sapply(insect.brains.with.optic, function(x) x$sex[[1]]))
df3 = data.frame(species = species.with.optic,
                 al.volume = insect.al.volumes[names(insect.optic.volumes)],
                 optic.volume = insect.optic.volumes
)

## Plot!
pdf("images/insectbraindb_optic_volume_vs_AL_volume.pdf", width = 10, height = 7)
ggplot2::ggplot(df3, aes(x=optic.volume, y=al.volume, color = species)) +
  geom_point() + 
  geom_smooth(data = df3, aes(x=optic.volume, y=al.volume), color  ="black", method=lm, se = FALSE)+
  theme_classic()
dev.off()

## Maybe we should normalise by total brain volume?
df4 = data.frame(species = species.with.optic,
                 al.norm.volume = insect.al.volumes[names(insect.optic.volumes)]/insect.brain.volumes[names(insect.optic.volumes)],
                 optic.norm.volume = insect.optic.volumes/insect.brain.volumes[names(insect.optic.volumes)]
)

## Plot!
pdf("images/insectbraindb_normalised_optic_volume_vs_normalised_AL_volume.pdf", width = 10, height = 7)
ggplot2::ggplot(df4, aes(x=optic.norm.volume, y=al.norm.volume, color = species)) +
  geom_point() + 
  geom_smooth(data = df4, aes(x=optic.norm.volume, y=al.norm.volume), color = "black", method=lm, se = FALSE)+
  theme_classic()
dev.off()

# Let's now plot AL size against Calyx size
calyx = c("Calyx",
          "Calyx basal zone",
          "Lateral calyx basal ring",
          "Lateral calyx collar",
          "Lateral calyx inner zone",
          "Lateral calyx lip",
          "Lateral calyx outer zone",
          "Medial calyx basal ring",
          "Medial calyx collar",
          "Medial calyx inner zone",
          "Medial calyx lip",
          "Medial calyx outer zone")
mblobes = c("Heel",
            "Medial gamma lobe",
            "Medial lobe",
            "Pedunculus",
            "Vertical gamma lobe",
            "Vertical lobe",
            "Y-lobes",
            "Dorsal lobelet",
            "Ventral lobelet",
            "Y-tract"
)

## Brains that contain a optic lobe neuropil objects, though not necessarily all of them
with.mb = sapply(insect.brains, function(ib) sum(ib$neuropil_full_names%in%calyx)>0 & sum(ib$neuropil_full_names%in%mblobes)>0 )
insect.brains.with.mb = insect.brains[with.mb]

## Calculate the neuropil volumes
insect.calyx.volumes = pbapply::pbsapply(insect.brains.with.mb, function(x)
  calculate_volume(subset(x, x$RegionList[x$neuropil_full_names%in%calyx])))
insect.mblobes.volumes = pbapply::pbsapply(insect.brains.with.mb, function(x)
  calculate_volume(subset(x, x$RegionList[x$neuropil_full_names%in%mblobes])))
insect.calyx.volumes[["Drosophila melanogaster UNKNOWN"]] = calculate_volume(subset(JFRC2NP.surf, "CA_"), alpha = 3)
insect.mblobes.volumes[["Drosophila melanogaster UNKNOWN"]] = calculate_volume(subset(JFRC2NP.surf, "MB_"), alpha = 3)
insect.brains.with.mb[["Drosophila melanogaster UNKNOWN"]] = JFRC2NP.surf

## Assemble data.frame
insect.brains.with.mb = insect.brains.with.mb[names(insect.brains.with.mb)!="Megalopta genalis FEMALE"]
species.with.mb = unlist(sapply(insect.brains.with.mb, function(x) x$common_name))
sex.with.mb = unlist(sapply(insect.brains.with.mb, function(x) x$sex[[1]]))
df5 = data.frame(species = species.with.mb,
                 al.norm.volume = insect.al.volumes[names(insect.brains.with.mb)]/insect.brain.volumes[names(insect.brains.with.mb)],
                 optic.norm.volume = insect.optic.volumes[names(insect.brains.with.mb)]/insect.brain.volumes[names(insect.brains.with.mb)],
                 calyx.norm.volume = insect.calyx.volumes[names(insect.brains.with.mb)]/insect.brain.volumes[names(insect.brains.with.mb)],
                 mblobes.norm.volume = insect.mblobes.volumes[names(insect.brains.with.mb)]/insect.brain.volumes[names(insect.brains.with.mb)],
                 total.volume = insect.brain.volumes[names(insect.brains.with.mb)]
)

## Plot!!
pdf("images/insectbraindb_normalised_mb_lobes_volume_vs_normalised_AL_volume.pdf", width = 10, height = 7)
ggplot2::ggplot(df5, aes(x=mblobes.norm.volume, y=al.norm.volume, color = species)) +
  geom_point() + 
  geom_smooth(data = df5, aes(x=mblobes.norm.volume, y=al.norm.volume), color = "black", method=lm, se = FALSE)+
  theme_classic()
dev.off()
pdf("images/insectbraindb_normalised_mb_calyx_volume_vs_normalised_AL_volume.pdf", width = 10, height = 7)
ggplot2::ggplot(df5, aes(x=calyx.norm.volume, y=al.norm.volume, color = species)) +
  geom_point() + 
  geom_smooth(data = df5, aes(x=calyx.norm.volume, y=al.norm.volume), color = "black", method=lm, se = FALSE)+
  theme_classic()
dev.off()
pdf("images/insectbraindb_normalised_mb_volume_vs_normalised_AL_volume.pdf", width = 10, height = 7)
ggplot2::ggplot(df5, aes(x=calyx.norm.volume+mblobes.norm.volume, y=al.norm.volume+optic.norm.volume, color = species)) +
  geom_point() + 
  geom_smooth(data = df5, aes(x=calyx.norm.volume+mblobes.norm.volume, y=al.norm.volume+optic.norm.volume), color = "black", method=lm, se = FALSE)+
  theme_classic()
dev.off()

## And let's see those neuropils
nat::nopen3d(userMatrix = structure(c(0.998683869838715, -0.00999264605343342, 
                                      -0.050302866846323, 0, -0.0184538997709751, -0.985155284404755, 
                                      -0.170669943094254, 0, -0.0478506684303284, 0.171373650431633, 
                                      -0.984043300151825, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.600000023841858, 
             windowRect = c(1440L, 45L, 3209L, 1063L))
for(i in names(insect.brains)){
  clear3d()
  ib = insect.brains[[i]]
  if(i %in% names(insect.brains.with.al)){
    message(ib$scientific_name, " the ", ib$common_name)
    plot3d(subset(ib, "^AL_left|^AL_right|^AL_noside|^AL_R$|^AL_L$"), col = "red")
    plot3d(ib, col = "lightgrey", alpha = 0.1)
    rgl.snapshot(file = paste0("images/insectbrainsAL/AL_",paste(unlist(strsplit(ib$common_name," ")),collapse="_"), "_",ib$sex, ".png"), fmt = "png")
    clear3d()
  }
  if(i %in% names(insect.brains.with.optic)){
    message(ib$scientific_name, " the ", ib$common_name)
    tryCatch(plot3d(subset(ib, ib$RegionList[if(is.null(ib$neuropil_full_names)){
      grepl("^LOP_|^LO_", ib$RegionList)
    }else{
      grepl("^Lobula",ib$neuropil_full_names)
    }]), col = "cyan"),
    error = function(e) NULL)
    tryCatch(plot3d(subset(ib, ib$RegionList[if(is.null(ib$neuropil_full_names)){
      grepl("^ME_", ib$RegionList)
    }else{
      grepl("^Medulla",ib$neuropil_full_names)
    }]), col = "blue"),
    error = function(e) NULL)
    plot3d(ib, col = "lightgrey", alpha = 0.1)
    rgl.snapshot(file = paste0("images/insectbrainsOL/OL_",paste(unlist(strsplit(ib$common_name," ")),collapse="_"), "_",ib$sex, ".png"), fmt = "png")
    clear3d()
  }
  if(i %in% names(insect.brains.with.mb)){
    message(ib$scientific_name, " the ", ib$common_name)
    plot3d(subset(ib, ib$RegionList[if(is.null(ib$neuropil_full_names)){
      grepl("^MB_CA_", ib$RegionList)
    }else{
      ib$neuropil_full_names%in%calyx
    }]), col = "purple")
    plot3d(subset(ib, ib$RegionList[if(is.null(ib$neuropil_full_names)){
      grepl("^MB_M|^MB_P|MB_V", ib$RegionList)
    }else{
      ib$neuropil_full_names%in%mblobes
    }]), col = "magenta")
    plot3d(ib, col = "lightgrey", alpha = 0.1)
    rgl.snapshot(file = paste0("images/insectbrainsMB/MB_",paste(unlist(strsplit(ib$common_name," ")),collapse="_"), "_",ib$sex, ".png"), fmt = "png")
    clear3d()
  }
  plot3d(bbxm, alpha = 0, add = TRUE)
  #rgl.snapshot(file = paste0("images/insectbrainsNP/NP_",paste(unlist(strsplit(ib$common_name," ")),collapse="_"), "_",ib$sex, ".png"), fmt = "png")
  #progress = readline(prompt = "Press any key for next brain ... ")
}

### Let's just plot the neuropil volumes side by side
m = reshape2::melt(df5, id.vars = c("species","total.volume"))
colnames(m) = c("species","total.volume","neuropil", "volume")
pdf("images/insectbraindb_normalised_neuropils_bar_plot.pdf", width = 8, height = 10)
ggdotchart(m, x = "species", y = "volume",
          fill = "neuropil", color = "neuropil", 
          add = "segment",
          add.params = list(color = "lightgray", size = 1.5),
          rotate = TRUE,
          dot.size = 6,                                 # Large dot size
          # position = position_dodge(0.1),
          ggtheme = theme_pubclean())+
  scale_color_manual(values = c(al.norm.volume = "red",
                                optic.norm.volume = "cyan",
                                calyx.norm.volume = "purple",
                                mblobes.norm.volume = "magenta"))+
  theme(legend.position = "none") + xlab("") + ylab("")
dev.off()
pdf("images/insectbraindb_normalised_neuropils_vs_total_volume.pdf", width = 10, height = 7)
ggscatter(m, x = "total.volume", y = "volume",
          fill = "neuropil", color = "neuropil", 
          add = "reg.line",                         # Add regression line
          conf.int = TRUE                          # Add confidence interval
)+stat_cor(aes(color = neuropil), label.x = 3)           # Add correlation coefficient
dev.off()

## Those budworms really skew things!

## Of course, to be sure, we would need to control by reconstruction method and work harder to establish the exact cell type correspondences.

