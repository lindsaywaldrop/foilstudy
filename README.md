# foilstudy -- An example of performance space creation and analysis

This repository accompanies the publication Hebdon, He, Battista, and Waldrop (in prep) 
and contains the necessary code to reproduce the figures therein. 

## Required Software

The following software is required for operating the code in this repository. Specific versions
are listed although other versions may also work. 

 - R version 4.4.1 (2024-06-14)
 - RStudio Version 2024.09.0+375
 - 

Note that the repository contains the results from the MATLAB simulations, so MATLAB is not required 
to generate surrogates or performance-space analyses. However, if you wish to reproduce the simulations 
themselves, MATLAB and XFOIL for MATLAB are required. 

 - MATLAB 2024a. 
 - XFOIL for MATLAB: https://www.mathworks.com/matlabcentral/fileexchange/50070-xfoil-for-matlab
Download and install by adding the installation folder to your MATLAB path. 

## Instructions for Reproducing Results

Please follow the steps below to reproduce the analyses:

 1. Clone the github repository or download the release and unarchive the code.
 2. Open the foilstudy.Rproj in RStudio. 

To reproduce the simulations (MATLAB required): 

 3. Open the `Creating_cambers.Rmd` file in `doc/` and run all lines. 
 4. Open the `xfoil_step.mlx` file in `doc/` MATLAB and run all sections. Note that this will take several hours to complete.
 5. Open the `Visualizing_results.Rmd` in `doc/` and run all lines. This will provide information on the raw results of the simulations.
 
To reproduce the surrogate construction and analysis: 



