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
best model. Its output are graphs that include surface brightness profiles of the 
galaxy and model. It also includes Surface brightness profiles of individual 
model components for a careful analysis of those. SB profiles includes 
the one along the major axis and a multi graph at different angles to analyse 
with detail inspection of the model. In addition, ``ElliSect`` complements 
the GALFIT photometry adding the total magnitud, luminosity, component to Total luminosity ratio, Akaike information criteriion, among other photometric variables. 

``EllipSect`` was designed to be used by GALFIT's users researchers of the 
astronomy community. GALFIT was designed to provide the user a quick decision over
the surface brigthness model and ``EllipSect`` is also designed to follow this. It
can provide the user to extract the most useful information of the model in order 
to assist the user to decide to add, remove or change model components. ``EllipSect`` 
has been used to analyse models of 2MASS and LINERs galaxies.  


[//]: <> (below is the sample paper: )

The forces on stars, galaxies, and dark matter under external gravitational
fields lead to the dynamical evolution of structures in the universe. The orbits
of these bodies are therefore key to understanding the formation, history, and
future state of galaxies. The field of "galactic dynamics," which aims to model
the gravitating components of galaxies to study their structure and evolution,
is now well-established, commonly taught, and frequently used in astronomy.
Aside from toy problems and demonstrations, the majority of problems require
efficient numerical tools, many of which require the same base code (e.g., for
performing numerical orbit integration).

``Gala`` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for ``Gala`` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. ``Gala`` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the ``Astropy`` package [@astropy] (``astropy.units`` and
``astropy.coordinates``).

``Gala`` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in ``Gala`` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike. The source code for ``Gala`` has been
archived to Zenodo with the linked DOI: [@zenodo]

# Acknowledgements

We acknowledge Chien Peng for their invaluable help thorough his Facebook's GALFIT page and the MGE fitting method and software by Cappellari (2002)
# References
