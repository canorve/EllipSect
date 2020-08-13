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
    orcid: 
    affiliation: 2


- name: Emmanuel Ríos-López
    orcid: 
    affiliation: "1, 3"



affiliations:
 - name: Facultad de Ciencias de la Tierra y el Espacio, Universidad Autónoma de Sinaloa, Blvd. de la Americas y Av. Universitarios S/N, Ciudad Universitaria, C.P. 80010 Culiacán, Sinaloa, México
   index: 1

 - name:  Instituto de Radioastronomía y Astrofísica, UNAM, Campus Morelia, AP 3-72, CP 58089, México
   index: 2

 - name: Instituto Nacional de Astrofísica Óptica y Electrónica (INAOE), Apartado Postal 51 y 216, 72000 Puebla, Mexico    
  index: 3


date: 1 August 2020
bibliography: paper.bib
---

# Summary

Galaxies are the building blocks of the large scale structure of the Universe. 
Those contain billions of stars that make up a diverse variety of galaxy morphology. As a consequence, they have various stellar components within galaxies such as bulges, bars, disks or rings among others. Hence, the image analysis for quantifying galaxy data is a fundamental step to understanding the formation, structure, and composition of the galaxies. One way to analyze them is through model fitting of the light distribution of a galaxy. Such models are mathematical functions of surface brightness that vary for the different components of the galaxies. A suitable model that reliably represents the physical properties requires a detailed inspection of the fitted models.

A well-known program to fit stellar surface brightness models is GALFIT ([@peng02] 1752 cites at the moment of writing this paper). It provides a wide variety of standards  functions such as Sersic [@sersic68], de Vaucouleurs [@devau48], Nuker, gaussian etc. To check if the galaxy model is a good fit, GALFIT provides a file with the fitted model parameters, errors, and a FITS[^](Flexible Image Transport System) cube image. The latter contains the galaxy, model and residual images. Typically, a visual check of those images can be tricky since it relies on image contrast. The provided data can be 
used to carry out a deeper analysis of the model and its residuals. 

GALFIT's users have been using graphs of surface brightness vs. radius to guide the eye for deviations from the galaxy and the model. They have been using IRAF's[^](Image Reduction and Analysis Facility) task ellipse [@jed87] which is another well-know program to extract surface brightness profiles. This process needs the format translation from GALFIT to Ellipse. This take time if the user needs to test various models to select the appropriate one for the galaxy. An additional issue is, unfortunately, that the development and maintenance of IRAF is discontinued since 2013. Nowadays, IRAF is actually supported by the astronomy community. 

``EllipSect`` is a python tool for making surface brightness profiles and extracting complementary photometry from the GALFIT output. The goal is to provide the greatest amount of information to select the best model. Its outputs are graphs that include surface brightness profiles of the galaxy and model. It also includes Surface brightness profiles of individual model components for a careful analysis. 
Furthermore, ``EllipSect`` complements the GALFIT photometry by adding other data that besides the ones extracted from the model's parameters, such as the total magnitud, luminosity, component to total luminosity ratio among others photometric variables. 

``EllipSect`` was designed to be easy to use for any researcher from the 
astronomy community. Users can make use of the program to decide to add, remove or change model components. It omits any direct interaction with the code or translation of GALFIT's data format. ``EllipSect`` has been used to analyse galaxy images from 2MASS and LINERs galaxies.  

Using wrapping scripts for GALFIT has been used before [@haussler13], but they cover other needs. For instance, their code run GALFIT to fit thousands of objects without user interaction on large images.


# ``ELLIPSECT``

``EllipSect`` only requires GALFIT output file. In this simple mode, ``EllipSect`` makes two graphs: one with the average of surface brightness along major axis, and the other contains multiple plots showing the surface brightness at different angles. See figure 1 for a fit of 7 gaussian models for an elliptical galaxy.  

![EllipSect output sample for an elliptical galaxy that was fitted with 7 gaussian models. In both panels red color represents the galaxy and blue the GALFIT model. Left panel: Surface brightness average vs. radius of both galaxy and model along the major axis. Model also has error bars since it is the average of individual model components. Right panel: multi plot of surface brightness of galaxy and model at different angles from major axis (major axis is the one with $0\deg$). At the right of the surface brightness plot is the one showing the percentage error at that angle. ](Fig1.png)


## Different Modes

``EllipSect`` has different input options to modify the original plots or 
compute other photometric variables. The user can manipulate the plot, and insert 
the surface brightness of the model sub-components. The program creates output files of the data of surface brightness of the graphs.

Below is shown a summary of the different features of ``EllipSect``:

- **Components** If a surface brightness model is composed from multiple components, users can enable ``ElliSect`` to plot the surface brightness of each sub-component. This allows the user to check if sub-component if fitted as desired.

- **Sky** The program uses the GALFIT's sky value from the input file, but alternatively, the user can enter their own sky value. Furthermore, just for comparison, ``EllipSect`` can calculate the sky background in the sector region where the gradient turns positive. 


- **Extra photometric variables.**  ``EllipSect`` can compute photometric variables that are indirectly extracted from the model parameters. Examples of such variables are: total magnitude, flux, mean surface brightness at effective radius, radius at 90% of total light, bulge to total ratio, and component to total light ratio.
  
- **Photometric aperture parameters.** The program uses _sectors\_photometry_ [@cappellari02] which divides an ellipse into sectors around to compute the counts in each sector. ``EllipSect`` computes, within this same ellipse, the following parameters: Tidal [@tal09], Bumpiness [@blakeslee06], SNR, Chinu, Akaike information criterion [@akaike74], Bayesian information criterion [schwarz78].

- **NED** ``EllipSect`` connects to NED[^](The NASA/IPAC Extragalactic Database (NED) is funded by the National Aeronautics and Space Administration and operated by the California Institute of Technology.) to download data of the galaxy to estimate absolute magnitude, luminosity, galactic extinction, distance modulus, cosmology corrected scale and surface brightness dimming.  
 
 

# Command line execution

The program is easily executed via the command line. It only requires 
the latest GALFIT output file (galfit.XX where XX is the largest number if there are various galfit.XX files). Example: 

``` 
./EllipSect galfit.01
``` 

This will generate the two graphs like the ones shown in figure 1. The _-help_ option will show additional features of the program.

# Future

This is part of a larger project where this program will be adaptated to analyze 
data images that includes hundreds of galaxies such as galaxy clusters. 

# Acknowledgements

We acknowledge Chien Peng for their invaluable help thorough his Facebook's GALFIT page, and the MGE fitting method and software by Cappellari (2002)

# References
