---
title: 'EllipSect: A Python tool for surface brightness model analysis for GALFIT'
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

 - name:  Instituto de Radioastronomia y Astrofísica, UNAM, Campus Morelia, AP 3-72, CP 58089, México
   index: 2

 - name: Instituto Nacional de Astrofísica Óptica y Electrónica (INAOE), Apartado Postal 51 y 216, 72000 Puebla, Mexico    
  index: 3


date: 1 August 2020
bibliography: paper.bib
---

# Summary

Galaxies are the building blocks of the large scale structure of the Universe. 
These contain billions of stars that make up a diverse variety of stellar components within galaxies such as bulges, bars, disks or rings. Hence, the study of those is a fundamental step to understanding the formation, structure, and composition of the galaxies. One way to analyze these components is fitting the components of a galaxy image using a stellar model. Such models are mathematical functions of surface brightness that vary for the different components of the galaxies. A good stellar model that truly represents the physical components requires a careful analysis of the fitted models of the galaxy.

A well-known program to fit stellar surface brightness model is GALFIT ([@peng02] 1752 cites at the moment of writing this paper). It provides a wide variety of very well-known  surface brightness functions such as Sersic, de Vaucouleurs, Nuker, moffat, gaussian etc. In order for the user to check whether galaxy model is a good fit, it provides a file with the fitted model parameters and a FITS (Flexible Image Transport System) cube image. This contains the galaxy image, model image and residual image. Sometimes, a visual check of the model and residual image can be tricky since it relies on image contrast. This is not enough for model analysis since it requires a plot analysis of the surface brightness profile vs. radius to compare galaxy and model, and a residual diagram.  

On the other hand, another well-know common used package to extract surface brightness  
from galaxy's isophotes is IRAF's ellipse [@jed87]. This method is one-dimensional and have its caveats. Development and maintenance of IRAF is discontinued since 2013.Nowadays, IRAF is actually supported by the astronomy community. 

Most GALFIT's users have been using IRAF's Ellipse tool to create surface brightness profiles of their models to compare with the one of the galaxies. This process need the translation of the output data from GALFIT to Ellipse. Besides, this needs a new configuration for Ellipse parameter's entry. All these can take time if the astronomer needs to fit several models to careful selects which one is the appropriate for the galaxy.


``EllipSect`` is a python tool to create surface brightness profiles and additional photometry from GALFIT output. The objective is to provide the most information to aid the user to select the best model. Its outputs are graphs that include surface brightness profiles of the galaxy and model. It also includes Surface brightness profiles of individual model components for a careful analysis of those. 
In addition, ``EllipSect`` complements the GALFIT photometry adding other that can not be extracted directly from the models parameters such as the total magnitud, luminosity, component to Total luminosity ratio among others photometric variables. 
(adds the galapagos software)

``EllipSect`` was designed to be used by GALFIT's users or any researcher from the 
astronomy community. GALFIT was designed to provide the user a quick decision over
the surface brigthness model and ``EllipSect`` is also designed to follow this objective. It provides the user the most useful information of the model in order to decide to add, remove or change model components. It does not require any direct interaction with the code or translation from GALFIT's data format. It uses directly 
the GALFIT output file which is helpful to analyze large amounts of galaxies.  ``EllipSect`` has been used to analyse models of 2MASS and LINERs galaxies.  


# ``ELLIPSECT``

The program was designed to be easy to use. It only requires 
GALFIT output file. In this simple mode, ``EllipSect`` makes two graphs: one with the average of surface brightness along major axis, and the other with multiple plots showing the surface brightness at different angles measured from major axis (i.e. major axis is the one with $0\deg$). See figure 1 for a fit of 7 gaussian models for an elliptical galaxy.  

Graphs of SB profiles includes the one along the major axis and a multi plot of SB at different angles to analyse the model with detail.  


![EllipSect output sample for an elliptical galaxy that was fitted with 7 gaussian models. In both panels red color represents the galaxy and blue the GALFIT model. Left panel: Surface brightness average vs. radius of both galaxy and model along the major axis. Model also has error bars since it is the average of individual model components. Right panel: multi plot of surface brightness of galaxy and model at different angles from major axis (major axis is the one with $0\deg$). At the right of the surface brightness plot is the one showing the percentage error at that angle. ](Fig1.png)


## Different Modes

``EllipSect`` has different input options which allows to modify the original plots or 
compute other photometric variables besides the ones computed by fitted surface brightness  models. Using different input options the user can: alter the X/Y range axis, include the surface brightness individual components to the plot, put a grid on the plot, change top X-axis to pixels. If desired, the program makes files of the surface brightness vs. radius values in such a way that the user can make its own 
plots in another graphic tool. 

Below is shown a summary of the different features that can be used with ``EllipSect``:

- **Components** If a surface brigthness model is composed from multiple models, users can enable ``ElliSect`` to plot the surface brightness of each model. This is important because users can check if an individual model has been property fitted to the desired galaxy component. 

- **Sky** The program already use the GALFIT sky value that was used in the fit, but alternatively,the user can introduce its own sky value. In addition, ``EllipSect`` can compute the sky background in the region area where the gradient becomes positive. The latter is not used for photometric computation (see below) and it is only shown for comparison. 

- **Additional photometric model parameters.**  ``EllipSect`` can compute photometric variables that are not directly extracted from the fitted model parameters. Examples of such variables are: total magnitude, flux, mean surface brigthness at effective radius, radius at 90% of total light, bulge to total ratio, and component light fraction of each component.
  
- **Photometric aperture.** The program uses _sectors\_photometry_ [@cappellari02] which divides into sectors an ellipse around the galaxy to compute the counts in each sectors. Taking advantage of this, ``EllipSect`` compute photometric variables within it. This means that the program computes: Tidal [@tal09], Bumpiness [@blakeslee06], SNR, Chinu (recomputed inside this ellipse), Akaike information criterion [@akaike74], Bayesian information criterion [schwarz78].

- **NED** ``EllipSect`` connects to NED[^](The NASA/IPAC Extragalactic Database (NED) is funded by the National Aeronautics and Space Administration and operated by the California Institute of Technology.) to extract information of the galaxy to compute absolute magnitude, luminosity, galactic extinction, distance modulus, cosmology corrected scale and surface brightness dimming.  
 
 

# Command line execution

Once downloaded, the program is easily executed via the command line. It only requires 
the latest output file that was generated with GALFIT (galfit.XX where XX is the largest number if there are several galfit.XX files). Example: 

``` 
./EllipSect galfit.01
``` 

This will generate the two graphs like the ones shown in figure 1. The _-help_ option will show additional features that can be used with the program.

# Future

This is part of a larger project where this program will be incorporated to analyze 
large data images which includes hundreds of galaxies such as galaxy clusters. 

# Acknowledgements

We acknowledge Chien Peng for their invaluable help thorough his Facebook's GALFIT page and the MGE fitting method and software by Cappellari (2002)

# References
