---
title: 'EllipSect: A Python tool for surface brightness analysis for GALFIT'
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
    orcid: 0000-0001-7707-7389
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
The larger ones contain billions of stars that can form one or various stellar components such as bulges, bars, disks and rings. Consequently, these make a wide variety in galaxy morphology. The quantification of the galaxy components is a fundamental step to understand their structure and composition. For instance, they can be compared with predictions of simulations of galaxy formation and see their evolution over time. Using models that fit their light distribution is one way to study those components. Such models are mathematical functions of surface brightness for the different components of the galaxies. A suitable model that reliably represents the physical properties requires a detailed inspection of the fitted models.

A well-known program for modeling the surface brightness of astronomical sources is GALFIT [[@peng02] 1753 cites at the moment of writing this document]. It allows to use a wide variety of standard functions such as Sérsic [@sersic68], de Vaucouleurs [@devau48], Nuker, gaussian, among others. GALFIT provides the fitted model parameters, errors, and a FITS (Flexible Image Transport System) cube image to check if the galaxy model is the appropriate one. The FITS file contains the galaxy, model and residual images. Typically, a visual check on the residual image can be difficult to interpret without additional information in addition to the $\chi^2_{\nu}$. Nevertheless, users construct surface brightness profiles from the GALFIT's output data to compare the galaxy and model images.

GALFIT users have been using plots of surface brightness vs. radius to guide the eye for deviations from the galaxy and the model. To do this, they have been using IRAF's (Image Reduction and Analysis Facility) task *ellipse* [@jed87], which is another well-known program to extract surface brightness profiles through ellipse fitting of the galaxy isophotes (regions of the galaxy where the surface brightness is constant). This process requires the data format translation from GALFIT to *ellipse*. This takes time if the user needs to test various models to select the appropriate one for the galaxy. An additional issue is that the development and maintenance of IRAF is discontinued since 2013. Nowadays, IRAF is actually supported by the astronomy community. 

Hence, we introduce ``EllipSect``, which is a Python tool to make surface brightness profiles and extract complementary photometry from the GALFIT's output. The program aids the users to select, remove or change model components. The goal is to provide as much information as possible to select the best model. ``EllipSect``'s outputs include graphs of the surface brightness profiles for the galaxy and the model. In case of multiple simultaneous galaxy fitting, it takes into account the surface brightness of neighbor galaxies to the galaxy of interest. This is unfeasible to do with IRAF's task *ellipse* since it only takes one at a time. It can also include the individual model components for a detailed analysis. Furthermore, ``EllipSect`` complements the GALFIT photometry by adding other data besides the ones extracted from the model's parameters, such as the total magnitude, luminosity, component to total luminosity ratio, among others photometric variables (see section below). 

Various scripts for GALFIT have been used before [@haussler13; @barden12; @anorve12; @vikram10], however, they cover other needs. For instance, their codes run GALFIT to automatically fit thousands of objects without user interaction on images containing multiple galaxies.

We designed ``EllipSect`` to be easy to use for any researcher from the 
Astronomy community. It omits any direct interaction with the code or translation of GALFIT's data format. ``EllipSect`` has been used to analyze galaxy images from *2MASS* (Two Micron All Sky Survey) in which estimations of  morphological and structural parameters have been obtained through photometric decompositions using GALFIT. Moreover,  the program is being used in a study of dust nuclear structures of a sample of Active Galactic Nuclei (AGN) in the local Universe through Hubble Space Telescope (HST) images, surface brightness models (with GALFIT) and radiative transfer simulations using SKIRT (advanced radiative transfer in dusty systems). 

# ``EllipSect``

The program only requires the GALFIT's output file. In this simple mode, ``EllipSect`` makes two graphs: one contains the surface brightness average along major axis, and the other one contains the surface brightness for different angles displayed in multiple plots.

The surface brightness is averaged through the perimeter of concentric ellipses along the major axis of the galaxy. The multiple plots aids the users to visualize where a fitted model fails to match with the galaxy. See figure 1 for an example using 7 Gaussian components for an elliptical galaxy.

![Example of EllipSect output for an elliptical galaxy and its model that was fitted with 7 gaussian components. In both panels the color red represents the galaxy and the blue one the GALFIT model. The color for each component is shown in the (top right) box. Left panel: Surface brightness average vs. radius along the major axis. The model also has error bars since it is the average of individual model components. Right panel: Multiple plots of surface brightness of galaxy and model at different angles from major axis (major axis is the one with $0\deg$). The error percentage is shown at the right side of the multi plot. ](Fig1.png)


## Different Modes

``EllipSect`` has different input options to modify the original plots or 
to compute other photometric variables. Additionally, the program creates 
files of the surface brightness data used in the graphs.

Below is shown a summary of the different features for ``EllipSect``:

- **Components**: If multiple models compose a surface brightness model, users can enable ``EllipSect`` to include the surface brightness of each sub-component. This allows the user to check if these are fitted as desired.

- **Sky**: The program uses the GALFIT's sky value from the input file, but alternatively, users can enter their own sky value. Furthermore, just for comparison, ``EllipSect`` can calculate the sky background using the slope of the image counts, and it returns the sky value in the region where the gradient turns positive. 

- **Complementary photometric variables**:  ``EllipSect`` can compute photometric variables that are indirectly extracted from the model parameters. For instance: total magnitude, flux, mean surface brightness at effective radius, radius at 90% of total light, bulge to total ratio, and component to total light ratio.
  
- **Photometric aperture parameters**: The program uses MGEfit library [@cappellari02] to divide an ellipse around the galaxy into sectors to compute the counts in each subregion. Within this ellipse, ``EllipSect`` calculates the following parameters: Tidal [@tal09], Bumpiness [@blakeslee06], Signal to Noise Ratio, $\chi^2_{\nu}$, Akaike information criterion [@akaike74], and Bayesian information criterion [@schwarz78].

- **NED**: ``EllipSect`` connects to NED (NASA/IPAC Extragalactic Database) to download data for the galaxy to estimate absolute magnitude, luminosity, galactic extinction, distance modulus, cosmology corrected scale and surface brightness dimming.  
 
 

# Command line execution

The program is easily executed via the command line. It only requires 
the latest GALFIT's output file. Example: 

``` 
./EllipSect galfit.01
``` 

This will generate the two graphs like the ones shown in figure 1. The _-help_ option will show additional features of the program.

# Future

This is part of a larger project where this program will be adapted to analyze 
data images that contain hundreds of galaxies such as galaxy clusters. 

# Acknowledgements

We acknowledge Chien Peng for his invaluable help through his GALFIT's Facebook page, as well as the MGE fitting method and software by Cappellari (2002).

# References
