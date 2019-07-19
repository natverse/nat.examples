# search for neurons
neocortex.df = neuromorpho_search(search_terms= c("brain_region:neocortex","archive:Jacobs"))
principal.df = subset(neocortex.df, grepl("principal",neocortex.df$cell_type))
interneuron.df = subset(neocortex.df, grepl("interneuron",neocortex.df$cell_type))

# read neurons
principals = neuromorpho_read_neurons(neuron_name = principal.df$neuron_name, batch.size = 5, nat = TRUE, progress = TRUE, OmitFailures = TRUE)
interneurons = neuromorpho_read_neurons(neuron_name = interneuron.df$neuron_name, batch.size = 5, nat = TRUE, progress = TRUE, OmitFailures = TRUE)

# save data
save("Jacobs_principal_neurons.rda")
save("Jacobs_interneurons.rda")

