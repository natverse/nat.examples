basic.nat
============

Some example code that demonstrate nat's basic functionality, looking at 'lateral horn neurons' of the fly brain. 

With this code, we explore data classes defined by nat, including EM neurons that contain synapse annotations.

It uses a package of fly neuron data called 'lhns'. This package is not part of the *natverse*, and is a little large. If you have an issue installing it (which can happen, for example, if you have a slow internet connections), you can try the following:

`remotes::install_git("git://github.com/jefferislab/lhns.git")`
(should take 5-30m depending on connection)

OR  in the terminal. Choose a directory that you like
`cd  ~/natverse`
`git clone --depth 1 https://github.com/jefferislab/lhns`

And then, in R
`remotes::install_local("~/natverse/lhns")`

The second methods takes serveral minutes to download, but it is much less fragile that the install_github mechanism and you get some indication of progress.

