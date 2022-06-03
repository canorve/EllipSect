.. contents::
   :depth: 3
..

**EllipSect: Full Output**
==========================

|DOI|

EllipSect creates surface brightness profiles and extracts other
photometric data from the GALFIT output peng et al. (2002).

**Script example**
------------------

If you want to use EllipSect inside your own python script, you can call
it like the following example:

::

       from ellipsect import ArgParsing 
       from ellipsect import SectorsGalfit

       #put all the argument parsing in a list:
       args=['galfit.01','--logx', '--phot','--noplot']


       parser_args = ArgParsing(args)

       photapi = SectorsGalfit(parser_args)

       print("Akaike Criterion: ",photapi.AICrit)
       print("Bulge to Total: ",photapi.BulgeToTotal)

In the previous example, the option “–phot” is necessary to produce the
output variables such as “photapi.AICirt”.

**Variables of the output class**
---------------------------------

Using the name photapi as an example (it can be any name), the full
output variables are explained below:

**variables that stores the names of the files**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| photapi.galfile: stores the name of the GALFIT file
| photapi.sboutput: name of the output surface brightness files
| photapi.output: name of the output photometry file
| photapi.inputmodel: name of the fits model file
| photapi.objname: name of the object or galaxy for search in NED,
  e.g. m51
| photapi.namefile: prefix name for the generation of output names
| photapi.namepng: name of the plot image.
| photapi.namesec: name of the sectors image. This is created by the
  sectors_photometry function
| photapi.namemul: name of the multiplot image.
| photapi.namemod: name of the sectors image for the model. This is
  created by the sectors_photometry function
| photapi.namesub: name of the individual components fits file
| photapi.namesig: name of the sigma image. This is created by GALFIT
| photapi.namesnr: name of the Signal-to-Noise ratio image
| photapi.namened: name of the xml NED file
| photapi.namecheck: name of the ellipse fits image where the photometry
  was done.
| photapi.namering: name of the ring file where the sky was computed.

**Variables that stores information of the GALFIT file**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| photapi.xc and photapi.yc: center of the galaxy photapi.q: galaxy axis
  ratio
| photapi.ang: position angle (same as GALFIT)
| photapi.skylevel: sky background
| photapi.scale: Plate scale arcsec/pixel
| photapi.inputimage: Name of the input image
| photapi.outimage: Name of the GALFIT output file
| photapi.maskimage: Name of the Mask image
| photapi.mgzpt: photometric zeropoint
| photapi.exptime: exposition time
| photapi.tempmask: temporary mask
| photapi.xmin, xmax, ymin, ymax : coordinates of the fitting region
  (same as GALFIT)
| photapi.band: photometric band. used for NED file

*The sky results given in the following variables are not involved in
the photometric output:*

| photapi.gradskymean: sky mean computed with the gradient method
| photapi.gradskystd: sky standard deviation with the gradient method
| photapi.gradskymed: sky median computed with the gradient method

| photapi.randskymean: sky mean computed with the random box method
| photapi.randskystd: sky standard deviation with the random box method
| photapi.randskymed: sky median computed with the random box method

**Variables that stores photometry of each GALFIT model components**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The variables below are arrays, because it stores the individual
components

| photapi.Comps: array boolean flags for components
| photapi.N: number of components
| photapi.NameComp: Name of the components
| photapi.PosX and PosY: X,Y Position of each individual components
| photapi.Mag: Magnitude of each individual compoment
| photapi.Rad: Radius of each components. If it is Sersic, this means
  the effective radius
| photapi.Exp: exponent of the component. If it is Sersic, this means
  the Sersic Index
| photapi.Exp2: Second exponent if the component have it, for example
  the moffat model
| photapi.Exp3: Third exponent
| photapi.AxRat: Axis ratio the individual component
| photapi.PosAng: Position angle of the individual component
| photapi.skip: parameter “z” of the GALFIT file refering to skip this
  model in the output image
| photapi.freepar: Number of free parameters of each individual
  components

| photapi.Rad50: Radius containing 50% of the component light. In case
  of Sersic, this is the same as effective radius
| photapi.SerInd: Sersic index
| photapi.Rad50kpc: Radius containing 50% of light in kpc
| photapi.Rad50sec: Radius containing 50% of light in arc sec
| photapi.Rad90: Radius containing 90% of light
| photapi.AbsMagComp: Absolute Magnitude of the component
| photapi.LumComp: Luminosity of the component
| photapi.Flux: Flux of the component
| photapi.PerLight: Percentage of the total light
| photapi.me: surface brightness at effective radius
| photapi.mme: mean surface brightness at effective radius
| photapi.kser: parameter coupled to n to allow to have surface
  brightness at Re
| photapi.KronRad: Kron radius
| photapi.PetRad: Petrosian radius computed with nu = 0.2

Take into account that there are some components that can not extracted
some of the variables. For instance, Petrosian, Effective radius can not
be computed for the Nuker model. In general all the extra variables is
computed for the Sersic models following the recipe of Graham et.
al. (2005)

**Variables that stores general photometry of the Galaxy**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some the quantities below are extracted from NED. Other quantities are
computed within the ellipse (to see this ellipse check the file ended in
\*-check.fits.

| photapi.aell: semi-major axis of the ellipse
| photapi.bell: semi-minor axis of the ellipse
| photapi.GalExt: Galactic extinction
| photapi.DistMod: Distance modulus
| photapi.DistMod2: Distance modulus (redshift independent)
| photapi.Scalekpc: arc sec to kilo parsec ratio
| photapi.SbDim: Surface brightness dimming
| photapi.magalaper: magnitude computed using ellipse’s aperture
| photapi.magmodaper: magnitude computed using ellipse’s aperture for
  the model
| photapi.totFlux: total flux of the galaxy
| photapi.totMag:total magnitude of the galaxy (sum of all components)
| photapi.BulgeToTotal: Bulge to total luminosity ratio
| photapi.tidal: Tidal value (check Tal et. al. (2009)
| photapi.objchinu: Chinu computed within ellipse (GALFIT computes this
  value for whole fitting region
| photapi.bump: bumpiness value see Blakeslee (2006)
| photapi.snr: mean of the Signal to noise ratio
| photapi.stdsnr: standard deviation of the signal to noise ratio
| photapi.totsnr: sum of the signal to noise ratio
| photapi.rss: residual sum of squares
| photapi.ndof: number of free paramters
| photapi.AbsMag: Absolute magnitude computed with DistMod
| photapi.AbsMag2: Absolute magnitude computed with DistMod2
| photapi.Lum: Luminosity of the galaxy (in solar luminosities)
| photapi.AICrit: Akaike Information Criterion
| photapi.BICrit: Bayesian Information Criterion
| photapi.BICres: Bayesian Information Criterion (using PSF area instead
  of number of pixels)

.. |DOI| image:: https://zenodo.org/badge/282223217.svg
   :target: https://zenodo.org/badge/latestdoi/282223217
