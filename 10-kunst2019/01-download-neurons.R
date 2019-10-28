# set the scene
## This script assumed that you have run the file "10-kunst2019/00-setup.R"
# download neuron data
urls=file.path("https://fishatlas.neuro.mpg.de/neurons/download/download_all_neurons")
message("Checking for presence of data ...")
for (url in urls){
  localfile=basename(url)
  if(file.exists(localfile)) next
  message("Downloading data ...")
  t=try(download.file(url, localfile, mode='wb'))
  if(inherits(t,'try-error')) {
    message("unable to download ", url)
    next
  }
}

# unzip
unzip("MPIN-Atlas__Kunst_et_al__neurons_all.zip")

# read neurons
zfish_neurons_right = read.neurons("MPIN-Atlas__Kunst_et_al__neurons_all/Right/")
zfish_neurons_original = read.neurons("MPIN-Atlas__Kunst_et_al__neurons_all/Original/")

# save
save(zfish_neurons_right, file = "zfish_neurons_right.rda")
save(zfish_neurons_original, file = "zfish_neurons_original.rda")
