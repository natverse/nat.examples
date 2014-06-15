zm <- readRDS('neurons.rds')

# Mirror all neurons past middle of x range
mid_x <- 248.857/2
neurons_to_mirror <- names(root_xs[root_xs > mid_x])
zm.canonical <- zm
zm.canonical[neurons_to_mirror] <- mirror(zm[neurons_to_mirror], warpfile="mirroring_registration.list/", mirrorAxisSize=248.857)

# Save results
saveRDS(zm.canonical, file="neurons_canonical.rds")
