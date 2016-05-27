# README #

This code simulates a simple explicit polymer model with various energy terms. These polymers can be stiff and sticky, and can be grafted onto a surface or isolated.

### Contents ###

* A class representing a collection of chains, with methods for perturbing the structure and checking its energy.
* A script that runs a Monte Carlo simulation of such an object.
* A VMD script for visualizing the resulting trajectory.

### Things that might be added ###

* Analysis code for evaluating properties like form factors, gyration radii, extensions, number of sticky bonds, etc, for finished runs.
* Simulated annealing for sampling the conformations of very sticky coils, which have rough energy landscapes.

### Documentation ###

* To install, download the code and run "chainSimulation.py". 
* Run "chainSimulation.py -help" for command-line usage.
* The underlying model is described in detail in the NXUS report for Statens Serum Institut (www.nxus.dk).