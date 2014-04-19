# Mouse retinal connectome data from from Helmstaedter et al 2013
## Article
Connectomic reconstruction of the inner plexiform layer in the mouse retina
 
Moritz Helmstaedter, Kevin L. Briggman, Srinivas C. Turaga, Viren Jain, H. Sebastian Seung & Winfried Denk
AffiliationsContributionsCorresponding author

Nature 500, 168–174 (08 August 2013) [doi:10.1038/nature12346](http://dx.doi.org/10.1038/nature12346)

## Warnings
The cell type labelled as "glia" needs verification. Only one unique skeleton 
for each cell is included in the packaged dataset – although 4-6 were typically
traced (each one by a different student).

## Running the code
In the R for Mac GUI (and perhaps others) you can drag and drop source code files
into the console.

```
source("00-setup.R", chdir=TRUE)
source("01-download.R", chdir=TRUE)
```
Then go over `02-sample-plots.R` step by step. You may need to do

```
setwd("/path/to/nat.examples/03-helmstaedter2013")
```

## TODO
* Convert to regular neuron rather than special skel datatype
* include full conversion script in R (currently need to dump single neurons out
  from matlab because R.matlab has trouble reading the original mat file)