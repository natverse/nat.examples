# This script assumed that you have run the file "16-skeletonisation/00-setup.R"

# Load the data on Planarian muscle
planarian.muscle = readobj::read.obj("data/19688.obj", convert.rgl  = TRUE)

# What do we have here?
print(names(planarian.muscle))

# A cell body fibre, a scalebar and a nucleus!
## We just want the fibre. Let's have a look at it,
planarian.fibre = planarian.muscle[[1]]
fibre_view()
plot3d(planarian.fibre, add = TRUE, col = lacroix["red"])
rgl.snapshot(filename = "images/original_planarian_muscle.png")

# Interesting
## To skeletonise, we need to make a complete cohesive mesh
### Let's do that with the alphashape package
alpha = 0.75
planarian.fibre.alpha = ashape3d(xyzmatrix(planarian.fibre), alpha = alpha, pert = TRUE)
plot(planarian.fibre.alpha, transparency = 0.1)
rgl.snapshot(filename = "images/planarian_muscle_alphashape_messy.png")
### Ah put look! There are lots of internal faces :(

# Looking good in genral though
### Now we have an issue that actually we have a lot of internal faces here
#### They may prove problematic, so let's purge them
planarian.fibre.alpha$triang = subset(planarian.fibre.alpha$triang,  planarian.fibre.alpha$triang[,"on.ch"]%in%c(0) | planarian.fibre.alpha$triang[,paste0("fc:",alpha)]%in%c(3) & !planarian.fibre.alpha$triang[,paste0("fc:",alpha)]%in%c(0))
planarian.fibre.alpha$edge = subset(planarian.fibre.alpha$edge, planarian.fibre.alpha$edge[,"on.ch"]%in%c(0) | planarian.fibre.alpha$edge[,paste0("fc:",alpha)]%in%c(3) & !planarian.fibre.alpha$edge[,paste0("fc:",alpha)]%in%c(0))
planarian.fibre.alpha$vertex = subset(planarian.fibre.alpha$vertex, planarian.fibre.alpha$vertex[,"on.ch"]%in%c(0) | planarian.fibre.alpha$vertex[,paste0("fc:",alpha)]%in%c(3) & !planarian.fibre.alpha$vertex[,paste0("fc:",alpha)]%in%c(0))
plot(planarian.fibre.alpha, transparency = 0.1)
rgl.snapshot(filename = "images/planarian_muscle_alphashape_clean.png")

# Now save .obj
fibre_view()
plot(planarian.fibre.alpha, transparency = 0.1)
writeOBJ("data/planarian_muscle.obj")

# Make into a mesh3d object, for easier use
planarian.fibre.mesh = readobj::read.obj("data/planarian_muscle.obj",convert.rgl = TRUE)[[1]]
clear3d()
plot3d(planarian.fibre.mesh, add = TRUE, col = lacroix["pink"], alpha = 0.25)
rgl.snapshot(filename = "images/planarian_muscle.png")

# Cool! A cohesive mesh
### Let's clean this mesh a bit though
#### To make skeletonisation later, easier
planarian.fibre.mesh.cleaned = vcgClean(planarian.fibre.mesh, sel = 0:6, iterate = TRUE)
planarian.fibre.mesh.cleaned = vcgIsolated(planarian.fibre.mesh.cleaned)
volume = tryCatch(Rvcg::vcgVolume(planarian.fibre.mesh.cleaned), error = function(e) NULL)
message("Estimated volume is: ", volume," units cubed")
clear3d()
plot3d(planarian.fibre.mesh.cleaned, add = TRUE, col = lacroix["pink"], alpha = 0.25)
rgl.snapshot(filename = "images/planarian_muscle.png")

## Now it is time to save this as a new.obj file
### And then use that to run our python-based skeletonisation.
fibre_view()
shade3d(planarian.fibre.mesh.cleaned)
writeOBJ("data/planarian_muscle_cleaned.obj")

# Now run skeletor to skeletonise it!
## The parameters may really take some playing with
planarian.muscle.skeleton = skeletor(obj = "data/planarian_muscle_cleaned.obj",
                                     validate = TRUE,
                                     radius = FALSE,
                                     method = "edge_collapse",
                                     shape_weight = 1,
                                     sample_weight = 0.5,
                                     mesh3d = FALSE,
                                     reroot = FALSE,
                                     clean = FALSE,
                                     iter_lim = 100)

# Let's see the meshes with its skeleton
fibre_view()
plot3d(planarian.muscle.skeleton, lwd = 4, col = lacroix["brown"])
plot3d(planarian.fibre.mesh.cleaned, add = TRUE, col = lacroix["pink"], alpha = 0.1)
rgl.snapshot(filename = "images/planarian_muscle_skeleton.png")

## Some stats
summary(planarian.muscle.skeleton)

## What if we just want the longest path?
simp = nlapply(planarian.muscle.skeleton, simplify_neuron, n = 0, invert = FALSE)
clear3d()
plot3d(simp, lwd = 4, col = lacroix["red"])
plot3d(planarian.muscle.skeleton, lwd = 4, col = lacroix["brown"])
rgl.snapshot(filename = "images/planarian_muscle_chosen.png")
leng = summary(simp)$cable.length
message("The length of this muscle fibre is ", leng, " microns")

## What if we do not want the cable length with the fork?
### If that's the case we can use nat::prune_online to interactively remove cable you do not want
#### I.e. pruned = prune_online(planarian.muscle.skeleton)
#### summary(pruned)



