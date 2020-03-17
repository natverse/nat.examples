# set the scene
## This script assumed that you have run the file "08-neuromorpho/00-setup.R"

# search for neurons
neocortex.df = neuromorpho_search(search_terms = c("brain_region:neocortex","archive:Jacobs"))
neocortex.df = subset(neocortex.df, species != "cricket")
principal.df = subset(neocortex.df, grepl("principal",neocortex.df$cell_type))
interneuron.df = subset(neocortex.df, grepl("interneuron",neocortex.df$cell_type))

# how many species?
table(principal.df$species)

# read neurons
principals = neuromorpho_read_neurons(neuron_name = principal.df$neuron_name, batch.size = 5, nat = TRUE, progress = TRUE, OmitFailures = TRUE)
interneurons = neuromorpho_read_neurons(neuron_name = interneuron.df$neuron_name, batch.size = 5, nat = TRUE, progress = TRUE, OmitFailures = TRUE)
principals = principals[principals[,"species"]!="cricket"]

# de-capitalise the species names
principals[,"species"] = tolower(principals[,"species"])
interneurons[,"species"] = tolower(interneurons[,"species"])

# save data
save(principals, file = "Jacobs_principal_neurons.rda")
save(interneurons, file = "Jacobs_interneurons.rda")

