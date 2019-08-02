nat.examples
============
<img align="right" width="300px" src="https://raw.githubusercontent.com/natverse/nat.examples/images/hex-natverse_logo.png">

Sample code for the NeuroAnatomy Toolbox ([nat](https://github.com/jefferis/nat)) R package

## Contents
01. Library installation and basic example
02. Grosjean et al 2011 Drosophila Projection Neuron Data
03. Helmstaedter et al 2013 Mouse Retinal Connectome
04. Sumbul et al 2014 Mouse Retinal Ganglion Cells
05. Miyasaka et al 2014 Zebrafish Mitral Cell Projectome
06. Lee et al 2012 Traced Drosophila neurons
07. Heinze et al 2013 Traced Monarch butterfly neurons

## Using these examples

1. Install [R](http://cran.r-project.org/) for your platform
2. Optionally install [RStudio IDE](http://www.rstudio.com/ide/download/)
3. Start R or RStudio
4. Now install (once only) and load nat package from within R

```
install.packages("nat")
library(nat)
```

Each example (except the very basic `01-setup`) assumes that R's current working directory has been set to the
relevant folder. This will be handled by the `00-setup.R` script in each folder if
it is called like this:

```
source('/path/to/nat.examples/02-grosjean20011/00-setup.R')
```
Do not use the `chdir=TRUE` option (even if your IDE wants to do this for you).

You can also set the path yourself manually e.g.

```
setwd('/path/to/nat.examples/02-grosjean20011/')
```

Once you have run the `00-setup.R` script you can start start running bits of
code interactively (e.g. with copy paste or RStudio's "Run selection* menu option).

## Acknowledgements
We insist that you cite the original authors of each study if you make use of
the data that they have released. Please also cite this repository 
(https://github.com/jefferis/nat.examples) if you make use of specific example
code along with the https://github.com/jefferis/nat package.
