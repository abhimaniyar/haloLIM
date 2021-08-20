# haloLIM

Code to calculate the spherically averaged power spectrum for [CII]
and CO(2-1) lines

Some warnings/advice:
* As of now, the code calculates the clustering and shot-noise
spherically averaged power spectrum for [CII] and CO2-1 lines
using some recipes (SFR to L_{[CII]} or L_{CO}) available.
The SFR required is calculated using the best fit parameterization
of Maniyar 2021 et al. (https://arxiv.org/pdf/2006.16329.pdf)
which was developed for a CIB halo model.

* In the current version of the code, we do not do a full halo-model
calculation. We just divide the power spectra in the clustering and
shot noise terms. However, it is pretty straightforward to implement
the halo model formalism and get the power spectra as 1-halo, 2-halo,
and shot noise terms.

* It has to be noted that currently the code is written such that it
calculates P(k, z) for several redshifts at the same time. Thus, if you
give only one redshift value as an input, the clustering power spectrum
will come out with the shape (len(k), 1). On the other hand, shot noise
is just a number for every redshift and not an array.

* If you want to run the code several times, or in a MCMC run, or just
in general want to speed up the calculations, pre-calculate and store
the halo mass function for a given mass and redshift range and then
just interpolate for the required mass and redshift ranges. Same
argument applicable to halo bias and especially the Fourier transform
of the NFW profile if you want to use the full halo model formalism.

* Requires numpy, scipy, and astropy for calculations and matplotlib
for plotting purposes. scipy and astropy dependancies are redundant
and code can actually be written fully using numpy if required.

You can just clone the repository, and run
driver.py file:
```
python driver.py
```

If you have any comments, suggestions or questions, please do not hesitate
to contact me: abhishek.maniyar@nyu.edu
