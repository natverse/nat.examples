### This script assumes you have run "/00-basics/00-setup.R"

# Let's see if we can cut through the main tract of some neurons, and then see how they are poistioned in that tract
## For extyra fun, we can compate a dataset derived from light microscopy, and one from electron micrographs.
## We can then see the benefits from a lot of data from a single brain, rather than sparse registered data from many.

# Connect to the public FAFB instance (Zheng et al. 2018) hosted publicly by Virtual Fly Brain
adult.conn = catmaid_login(server="https://catmaid-fafb.virtualflybrain.org/")

# Get all EM D. melanogaster uniglomerular olfactory projeciton neurons
em.pns = read.neurons.catmaid("annotation:^PN glomerulus")

# Find presynapses
pres = xyzmatrix(get.synapses(em.pns, "PRE"))
pin = points_in_neuropil(pres, brain = FAFB14NP.surf, alpha = 30000)
pinm = aggregate(pin$neuropil, list(pin$neuropil), length)
colnames(pinm) = c("neuropil", "presynapses")
pinm = subset(pinm, neuropil != 0)
pinm = subset(pinm, presynapses > 50)

# Plot!
pdf("images/nat_basic_EM_presynapses_in_neuropils.pdf", width = 20, height = 5)
ggplot(pinm, aes(x=neuropil, y=presynapses, fill = transmitter)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  ylab("no. presynapses") +
  theme(legend.position="none")  + 
  theme(axis.text.x = element_text(size=20,angle = 90, hjust = 1),axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20))
dev.off()

# Show neuropils of the fly brain
em.pns.fcwb = xform_brain(em.pns, sample = FAFB14, reference = FCWB) # bridge from EM to LM!
nopen3d(userMatrix = structure(c(0.998409867286682, -0.00922653265297413, 
                                 -0.0556101016700268, 0, -0.0129726231098175, -0.997642815113068, 
                                 -0.0673833340406418, 0, -0.0548573285341263, 0.0679976567625999, 
                                 -0.996176183223724, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.481017380952835, 
        windowRect = c(1460L, 65L, 3255L, 1048L))
plot3d(FCWB, alpha = 0.1)
plot3d(subset(FCWBNP.surf, "LH"), col = "green", alpha = 0.2)
plot3d(subset(FCWBNP.surf, "AL_"), col = "red", alpha = 0.2)
plot3d(subset(FCWBNP.surf, "MB_CA_"), col = "purple", alpha = 0.2)
plot3d(em.pns.fcwb, col = "chartreuse", lwd = 3, soma = TRUE)
rgl.snapshot(filename = "/images/FCWB_EM_OPNs_neuropils.png", fmt = "png")
