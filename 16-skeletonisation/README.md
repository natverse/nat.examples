# Skeletonising planarian muscle

## Article

The natverse package *fafbseg* contains an R wrapper for python code built by P. Schlegel, for skeletonising mesh data. We have used it extensively to skeletonise neurons. In this example, we examine a more general case - that of planarian muscle. [Planarians](https://en.wikipedia.org/wiki/Planarian) are small flatworms with remarkable regeneration capabilities.

This example uses a .obj file given to us by [James Cleland](https://twitter.com/jpcleland) from the [Rink lab](https://twitter.com/rinklab). It is from a [TrakEM2 project](https://imagej.net/TrakEM2), where James has reconstructed a muscle fibre in a flatworm.

In this example, we attempt to skeletonise the muscle fibre and discover its 'cable length'.

In order for this code to work, you will have to set yourself up to run python code for skeletonisation, i.e. [the skeletor library](https://github.com/schlegelp/skeletor) by P. Schlegel. Detailed instructions for this can be found [here](http://natverse.org/fafbseg/articles/articles/installing-cloudvolume-meshparty.html). The in those instructions has to do with working with [flywire](https://ngl.flywire.ai/) neurons.

## Running the code

See `Using these examples` in the main `README.md` file
