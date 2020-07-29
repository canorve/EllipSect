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
They are composed of several stellar structures like bulge, bar, disk rings. The 
study of those structures is a fundamental step to understanding the formation, history
and composition of the galaxies. A critical step to understand them is by 
the means of appropriate model selection of surface brightness. Such step
requires a careful analysis of the models that can be used to fit the galaxy properly. 


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

We acknowledge Chien Peng and  the MGE fitting method and software by Cappellari (2002)
# References
