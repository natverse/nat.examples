# Mouse retinal connectome data from from Helmstaedter et al 2013
## Article
Connectomic reconstruction of the inner plexiform layer in the mouse retina
 
Moritz Helmstaedter, Kevin L. Briggman, Srinivas C. Turaga, Viren Jain, H. Sebastian Seung & Winfried Denk
AffiliationsContributionsCorresponding author

Nature 500, 168–174 (08 August 2013) [doi:10.1038/nature12346](http://dx.doi.org/10.1038/nature12346)

## Warnings
The cell type labelled as "glia" needs verification

## Running the code

```
source("00-setup.R", chdir=TRUE)
```
etc.

## TODO
* Convert to regular neuron rather than special skel datatype
* include full conversion script in R (currently need to dump single neurons out
  from matlab because R.matlab has trouble reading the original mat file)