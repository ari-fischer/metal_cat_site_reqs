# metal_cat_site_reqs

This repository contains (A) a set of Python codes used to calculate the excluded areas of adsorbates on metal surfaces based on their molecular structures and volumes determined from van der Waals radii, and (B) MATLAB code for regressing kinetic parameters to fit turnover frequency data with multi-site kinetic model for benzene hydrogenation.

These codes accompany the manuscript: "Determining Site Requirements for Multi-Site Catalysis on Metal Surfaces Using Molecular Areas" published in the Journal of Catalysis (https://doi.org/10.1016/j.jcat.2025.116179)

(A) Two models are included:

(1) Direct projection of the adsorbate's van der Waals volume onto the xy-plane of the geometry file (corresponding to the catalyst surface plane)
(2) Area exclude by a circular probe tracing the van der Waals area of the adsorbate. One version considers an H-adatom as the probe, and another considers a circle with an area equivalent to the adsorbate considered.

It currently supports adsorbates composed of H, C, and O. 

The code requires the Atomic Simulation Environment (https://wiki.fysik.dtu.dk/ase/index.html) package to read geometry structure.


(B) MATLAB code that runs lsqcurvfit function to regress a temperature-dependent kinetic model to a set of turnover frequency measurements at different temperatures and benzene pressures. 

Version history:

Accepted manuscript (Journal of Catalysis) 28 April 2025
	Added the MATLAB codes (B)
