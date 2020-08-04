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
They are composed of stellar structures like bulges, bars, disks or rings. The 
study of those structures is a fundamental step to understanding the formation, history
and composition of the galaxies. A critical step to study them is by 
the means of appropriate model selection of surface brightness. Such step
requires a careful analysis of the models that can be used to fit 
the galaxy properly. 

A popular program to fit such stellar surface brightness model is GALFIT (10000 cites 
at the moment of writing this paper). It can provide a long variety of very well know 2D surface brightness such as Sersic, de Vaucouleurs, etc. It uses the Levenbergh-Marquart algorithm to find the optimal chinu that fits the model to 
the galaxy. In order to check whether galaxy model is a good fit, it provides 
a cube image which it contains the galaxy image, model image and residual image. The 
visual check of the model and residual image can be tricky since it relies on image 
contrast. If it is a good or bad model can be arbitrary to the user. This is not enough 
to model analysis.

GALFIT's users have been using IRAF's Ellipse tool to create surface brightness models
to compare galaxy and model and check if the model was a good fit. Nowadays IRAF is no
longer maintained and actually supported by the astronomy community.  Moreover, 
data translation between GALFIT output to ellipse can take time if the astronomer
needs to fit several model and careful model selection can take longer. 

``EllipSect`` is a python tool to analyse GALFIT output in order to select the 
best model. It reads the x,y center, axis ratio and angular position of the model to  divide the sectors of an ellipse and compute the surface brightness in those individual sectors. Its outputs are graphs that include surface brightness profiles of the 
galaxy and model. It also includes Surface brightness profiles of individual 
model components for a careful analysis of those. SB profiles includes 
the one along the major axis and a multi graph at different angles to analyse 
with detail inspection of the model. In addition, ``ElliSect`` complements 
the GALFIT photometry adding the total magnitud, luminosity, component to Total luminosity ratio, Akaike information criteriion, among other photometric variables. 

``EllipSect`` was designed to be used by GALFIT's users researchers of the 
astronomy community. GALFIT was designed to provide the user a quick decision over
the surface brigthness model and ``EllipSect`` is also designed to follow this. It
can provide the user to extract the most useful information of the model in order 
to assist the user to decide to add, remove or change model components. It does not
require to any direct interaction with the code. ``EllipSect`` 
has been used to analyse models of 2MASS and LINERs galaxies.  


# ``ELLIPSECT``

The program was designed simple to use. It only requires GALFIT output file (galfit.XX).
It reads the model In this simple mode, ``EllipSect`` makes two graphs: one with the average of surface brightness along major axis, and the other with multiple plots showing the surface brightness at different angles measured from major axis (i.e. major axis is the one with $0\deg$). See figure 1 for a fit of 7 gaussian models for an Elliptical galaxy.  

![EllipSect output sample for an elliptical galaxy that was fitted with 7 gaussian models. In both panels red color represents the galaxy and blue the GALFIT model. Left panel: Surface brightness average vs. radius of both galaxy and model along the major axis. Model also has error bars since it is the average of individual model components. Right panel: multi plot of surface brightness of galaxy and model at different angles from major axis (major axis is the one with $0\deg$). At the right of the surface brightness plot is the one showing the percentual error at that angle. ](Fig1.png)


## Different Modes

``EllipSect`` has different input option which allows to modify the original plots or 
compute other photometric variables besides the ones computed by the GALFIT models.
By default the program plots the graphs shown in Figure 1. Using different input options the user can: alter the X/Y range axis, include the surface brightness individual components to the plot, put a grid on the plot, change top X-axis to pixels. The program can also makes surface brightness to output to a file in such a way that 
the user can make his/her own plots in another graphic tool. 

Below is shown a summary of the different features that can be used with ``EllipSect``:

- **Sky** The programs reads GALFIT sky that was used in the fit, but alternatively user can introduce their own sky value. In addition, ``EllipSect`` can compute the sky background in the region area where the gradient becomes positive. This is not used in the photometry output (see below) and it is only shown for comparison. 

- **Additional photometric model parameters.**  If the user enables it, ``EllipSect`` can compute photometric variables that are not directly extracted from the fitted model parameters. Examples of such variables are: total magnitude, flux, mean surface brigthness at effective radius, radius at 90% of total light, bulge to total ratio, and light fraction per component. 

- **Photometric aperture.** The program uses _sectors\_photometry_ (cite) which divides into sectors an ellipse around the galaxy to compute the counts in those individual sectors. Taking advantage of this already created ellipse, ``EllipSect`` compute photometric variables within it. Those include: Tidal, Bumpiness, SNR, Chinu (recomputed inside this ellipse), Akaike information criterion, Bayesian information criterion (cites).

- **NED** ``EllipSect`` connects to NED (cite) to extract information of the galaxy needed to compute absolute magnitude, luminosity, galactic extinction, distance modulus, cosmology corrected scale, surface brightness dimming.  

# Command line execution

Once downloaded, the program is easily executed via the command line. It only requires 
the latest output file that was generated with GALFIT (galfit.XX where XX is the largest number if multiple runs have been using with GALFIT). Example: 

./EllipSect galfit.01

This will generate the two graphs shown above. The _-help_ option will show additional
features that can be used with the program.

# Future

This is part of a larger project where this program will be incorporated to analyze 
large data images which includes hundreds of galaxies such as galaxy clusters. 

# Acknowledgements

We acknowledge Chien Peng for their invaluable help thorough his Facebook's GALFIT page and the MGE fitting method and software by Cappellari (2002)
# References
