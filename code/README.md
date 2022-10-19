# README

This code simulates a simple explicit polymer model with various energy terms. These polymers can be stiff and sticky, and can be grafted onto a surface or isolated.

## Model
The physical model is one where $N$ hard spheres are connected along a chain. The distance between adjacent beads is $d$, and angle between three adjacent beads is free up to a maximum angle $\theta_\mathrm{max}$. Non-adjacent beads interact via a square-well potential of depth $\eps$, and a width such that beads are considered bound when they are closer together than $1.2d$. In the sampling of the model, a Metropolis condition based on $e^{-\Delta n \cdot \epsilon /RT}$ is applied for each step where $\Delta n$ bonds are broken.

## Contents

* A class representing a collection of chains, with methods for perturbing the structure and checking its energy.
* A script that runs a Monte Carlo simulation of such an object, with possible simulated annealing.
* A VMD script for visualizing the resulting trajectory.

## Things that aren't included but could have been

* Analysis code for evaluating properties like form factors, gyration radii, extensions, number of sticky bonds, etc, for finished runs.
* Parallelize the damn thing.
* Add a Monte Carlo routine for randomly placing chains densely on a surface.

## Documentation

* To install, download the code and run "chainSimulation.py". 
* Run "chainSimulation.py -help" for command-line usage.
* The underlying model is described in detail in the NXUS report for Statens Serum Institut (www.nxus.dk).