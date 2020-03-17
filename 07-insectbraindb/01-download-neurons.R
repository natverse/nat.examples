## This script assumed that you have run the file "07-insectbraindb/00-setup.R"

## Functions to talk to insectbraindb.org, a site primarily curated by Prof. Stanley Heinze, are wrapped up into the package neuromorphr

## So let's use that

## What neurons does the insectbraindb.org host?
available.neurons = insectbraindb_neuron_info()

## Let's just download all of the neurons in the database to play with,
## there are not very many:
nrow(available.neurons)

## First, we call the read neurons function, with ids set to NULL
insect.neurons = insectbraindb_read_neurons(ids = NULL)

## Hmm, let's see how many neurons we have perspecies
table(insect.neurons[,"common_name"])

## So, it seem the Monarch Butterfly is the clear winner there, 
## maybe let's just have those
butterfly.neurons = subset(insect.neurons, common_name == "Monarch Butterfly")

