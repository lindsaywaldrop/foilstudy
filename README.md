# foilstudy -- An example of performance space creation and analysis

This repository accompanies the publication Hebdon, He, Battista, and Waldrop (in prep) 
and contains the necessary code to reproduce the figures therein. 

## Required Software

The following software is required for operating the code in this repository. Specific versions
are listed although other versions may also work. 

 - R version 4.4.1 (2024-06-14)
 - RStudio Version 2024.09.0+375
 - R Packages (available on CRAN): `pracma`, `ggplot2`, `patchwork`, `spacefillr`

Note that the repository contains the results from the MATLAB simulations, so MATLAB is not required 
to generate surrogates or performance-space analyses. However, if you wish to reproduce the simulations 
themselves, MATLAB and XFOIL for MATLAB are required. 

 - MATLAB 2024a. 
 - XFOIL for MATLAB: https://www.mathworks.com/matlabcentral/fileexchange/50070-xfoil-for-matlab
Download and install by adding the installation folder to your MATLAB path. 

## Instructions for Reproducing Results

Please follow the steps below to reproduce the analyses. Note that several steps can be performed with either R or MATLAB, depending on your preference.

 1. Clone the github repository or download the release and unarchive the code. If running in RStudio, open the `foilstudy.Rproj` in RStudio. 
 2. Generate the parameter spaces using grid sampling, gPC sampling, and neural network sampling.
    * MATLAB: open and run all lines `./doc/Generating_points_MATLAB.mlx` in MATLAB. 
    * R: open and run all lines `./doc/Generating_points.Rmd` in RStudio. (Note: NN sampling is not yet implemented in R!)
 3. Generate the airfoil files required to run in XFOIL. 
    * MATLAB: Open the `./doc/Creating_cambers_MATLAB.mlx` in MATLAB and run all lines.
    * R: Open the `./doc/Creating_cambers.Rmd` file in RStudio and run all lines, including the last code chunk which requires that a line be uncommented.  
 4. (MATLAB Only) Open the `./doc/xfoil_step.mlx` file MATLAB and run all sections. Note that this will take several hours to complete. Results are included in the project folder `./results/`, so it is not necessary to complete this step if only using R.
 5. Open the `Visualizing_results.Rmd` in `doc/` and run all lines. This will provide information on the raw results of the simulations.
 
To reproduce the surrogate construction and analysis: 



