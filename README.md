# metal_cat_site_reqs

This repository contains a set of Python codes used to calculate the excluded areas of adsorbates on metal surfaces based on their molecular structures and volumes determined from van der Waals radii. Two models are included:
(1) Direct projection of the adsorbate's van der Waals volume onto the xy-plane of the geometry file (corresponding to the catalyst surface plane)
(2) Area exclude by a circular probe tracing the van der Waals area of the adsorbate. One version considers an H-adatom as the probe, and another considers a circle with an area equivalent to the adsorbate considered.

It currently supports adsorbates composed of H, C, and O. 

The code requires the Atomic Simulation Environment (https://wiki.fysik.dtu.dk/ase/index.html) package to read geometry structure.

Codes accompanying manuscript: "Determining Site Requirements for Multi-Site Catalysis on Metal Surfaces Using Molecular Areas"

Version history:
