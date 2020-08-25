---
title: 'EllipSect: A Python tool for surface brightness analysis'
tags:
  - Python
  - astronomy
  - galaxies
  - surface brightness
  - photometry
authors:
  - name: Christopher Añorve
    orcid: 0000-0002-3721-8869
    affiliation: 1
  - name: Omar U. Reyes-Amador
    affiliation: 2
  - name: Emmanuel Ríos-López
    orcid: 0000-0002-4436-221X
    affiliation: "1, 3"

affiliations:
 - name: Facultad de Ciencias de la Tierra y el Espacio, Universidad Autónoma de Sinaloa, Blvd. de la Americas y Av. Universitarios S/N, Ciudad Universitaria, C.P. 80010 Culiacán, Sinaloa, México
   index: 1
 - name:  Instituto de Radioastronomía y Astrofísica, UNAM, Campus Morelia, AP 3-72, CP 58089, México
   index: 2
 - name: Instituto Nacional de Astrofísica Óptica y Electrónica (INAOE), Apartado Postal 51 y 216, 72000 Puebla, Mexico    
   index: 3
date: 21 August 2020
bibliography: paper.bib
---

# Summary

Galaxies are the building blocks of the large scale structure of the Universe. 
The larger ones contain billions of stars that form stellar components such 
as bulges, bars, disks and rings. Consequently, these components make a diverse variety in galaxy morphology. The quantification of the galaxy images is a fundamental step to understand their structure and composition [emmanuel comentario]. Using models that fit their light distribution is one way to do it. Such models are mathematical functions of surface brightness for the different components of the galaxies. A suitable model that reliably represents the physical properties requires a detailed inspection of the fitted models.

A well-known program for modeling the surface brightness of astronomical sources is GALFIT [[@peng02] 1752 cites at the moment of writing this paper]. It allows to use a wide variety of standard functions such as Sérsic [@sersic68], de Vaucouleurs [@devau48], Nuker, gaussian, among others. GALFIT provides the fitted model parameters, errors, and a FITS (Flexible Image Transport System) cube image to check if the galaxy model is the appropriate one.
The FITS file contains the galaxy, model and residual images. Typically, a visual check on those images can be difficult since it relies on image contrast. Nevertheless, users construct surface brightness profiles from the GALFIT's output data to compare the galaxy and model images.

GALFIT's users have been using plots of surface brightness vs. radius to guide the eye for deviations from the galaxy and the model. To do this, they have been using IRAF's(Image Reduction and Analysis Facility) task *ellipse* [@jed87] which is another well-know program to extract surface brightness profiles through ellipse fitting of the galaxy isophotes (regions of the galaxy where the surface brightness is constant). This process requires the data format translation from GALFIT to ellipse. This take time if the user needs to test various models to select the appropriate one for the galaxy. An additional issue is that the development and maintenance of IRAF is discontinued since 2013. Nowadays, IRAF is actually supported by the astronomy community. 

Hence, we introduce ``EllipSect`` is a Python tool to make surface brightness profiles and extract complementary photometry from the GALFIT output. The program aids the users to select, remove or change model components. The goal is to provide the most information to select the best model. ``EllipSect`` outputs include graphs of the surface brightness profiles for the galaxy and model. For multiple galaxy fits, it takes into account the surface brightness of nearby galaxies [emmanuel comentario]. This is unfeasible to do with IRAF's task ellipse. It can also includes the individual model components for a detailed analysis. Furthermore, ``EllipSect`` complements the GALFIT photometry by adding other data besides the ones extracted from the model's parameters, such as the total magnitude, luminosity, component to total luminosity ratio, among others photometric variables (see section below). 

Various scripts for GALFIT have been used before [@haussler13], however they cover other needs. For instance, their code run GALFIT to fit thousands of objects without user interaction on large images [emmanuel comentario DGCG].

We designed ``EllipSect`` to be easy to use for any researcher from the 
astronomy community. It omits any direct interaction with the code or translation of GALFIT's data format. ``EllipSect`` has been used to analyze galaxy images from 2MASS (Two Micron All Sky Survey) and LINERs galaxies in which estimations of  morphological and structural parameters have been obtained through photometric decompositions using GALFIT [emmanuel comentario].  


# ``ELLIPSECT``

``EllipSect`` only requires the GALFIT output file. In this simple mode, ``EllipSect`` makes two graphs: one contains the surface brightness average along major axis, and the other contains the surface brightness for different angles displayed in multiple plots. See figure 1 for a fit of 7 gaussian models [emmanuel comentario] for an elliptical galaxy [emmanuel comentario].  

![Example of EllipSect output for an elliptical galaxy and its model that was fitted with 7 gaussian components. In both panels red color represents the galaxy and blue the GALFIT model. Left panel: Surface brightness average vs. radius along the major axis. Model also has error bars since it is the average of individual model components. Right panel: multiple plots of surface brightness of galaxy and model at different angles from major axis (major axis is the one with $0\deg$). The percentage error is shown at the right side of the multi plot. ](Fig1.png)


## Different Modes

``EllipSect`` has different input options to modify the original plots or 
to compute other photometric variables. Additionally the program creates 
files of the surface brightness data used in the graphs.

Below is shown a summary of the different features for ``EllipSect``:

- **Components**: If multiple models compose a surface brightness model, users can enable ``ElliSect`` to include the surface brightness of each sub-component. This allows the user to check if these are fitted as desired.

- **Sky**: The program uses the GALFIT's sky value from the input file, but alternatively, users can enter their own sky value. Furthermore, just for comparison, ``EllipSect`` can calculate the sky background in the sector region where the gradient turns positive. 

- **Complementary photometric variables.**:  ``EllipSect`` can compute photometric variables that are indirectly extracted from the model parameters. For instance: total magnitude, flux, mean surface brightness at effective radius, radius at 90% of total light, bulge to total ratio, and component to total light ratio.
  
- **Photometric aperture parameters.**: The program uses _sectors\_photometry_ from MGEfit library [@cappellari02], which divides an ellipse around the galaxy into sectors to compute the counts in each subregion. Using this ellipse, ``EllipSect`` calculates the following parameters: Tidal [@tal09], Bumpiness [@blakeslee06], Signal to Noise Ratio, $\chi^2_{\nu}$, Akaike information criterion [@akaike74], Bayesian information criterion [@schwarz78].

- **NED**: ``EllipSect`` connects to NED (the NASA/IPAC Extragalactic Database (NED) is funded by the National Aeronautics and Space Administration and operated by the California Institute of Technology.) to download data for the galaxy to estimate absolute magnitude, luminosity, galactic extinction, distance modulus, cosmology corrected scale and surface brightness dimming.  
 
 

# Command line execution

The program is easily executed via the command line. It only requires 
the latest GALFIT output file. Example: 

``` 
./EllipSect galfit.01
``` 

This will generate the two graphs like the ones shown in figure 1. The _-help_ option will show additional features of the program.

# Future

This is part of a larger project where this program will be adapted to analyze 
data images that contains hundreds of galaxies such as galaxy clusters. 

# Acknowledgements

We acknowledge Chien Peng for their invaluable help thorough his Facebook's GALFIT page, and the MGE fitting method and software by Cappellari (2002)

# References
