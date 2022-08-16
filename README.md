
# Modelling the demographic history of _A. thaliana_
## Project for Osmond Lab

This project is intended to provide a test case for Matt Osmond's [_sparg_](https://github.com/mmosmond/sparg.git) software, which estimates dispersal rates and locates genetic ancestors from tree sequence data using genome-wide genealogies.

First, I created an individual-based forward simulation roughly approximating the demographic history of _A. thaliana_, as per [Fulgione and Hancock (2018)](https://doi.org/10.1111/nph.15244), coded up using the R-package _slendr_ and the forward simulation software SLiM. The genomes of the present day populations were then sampled throughout the population ranges and analyzed on Python, using _sparg_, to infer the dispersal rate and origin of the individuals. This was then compared to the actual dispersal rate and origin used to model the populations, giving us an idea of how well _sparg_ performs on genomes sampled from populations with complex histories.

#### I. Setting up the model on _slendr_
The R script used to set up the basic geographical boundaries and population dynamics can be found in the home directory. The files in the model directory generated by this R script can be found in model/. Note that there are two SLiM scripts provided here: script_original.slim is the default backend skeleton script that _slendr_ uses to run the forward simulations on SLiM, and script.slim is the edited SLiM script (see below) that can be run from the command line (without _slendr_) using a dedicated bash script.

#### II. Bash scripts that edit and run the SLiM script
The SLiM script generated by this code does not capture the selfing behavior of _A. thaliana_, so we next wrote a bash script that edits the SLiM script to get the behavior we want. 

Since _slendr_ is not available on CRAN yet, installing it on a cluster is somewhat tedious. Thus, we wanted a simple bash script that uses the model directory and the edited SLiM script to do what the ```slim()``` function would do in _slendr_ - given some parameter values, run the simulation and produce the sampled tree sequence outputs at the end.

These bash scripts can be found in bash/.

#### III. Using _sparg_ to analyze the tree sequences
To be completed.

#### IV. Interpreting the results
To be completed.

#### Citations
1. Fulgione, A. and Hancock, A.M. (2018), Archaic lineages broaden our view on the history of Arabidopsis thaliana. New Phytol, 219: 1194-1198. https://doi.org/10.1111/nph.15244
