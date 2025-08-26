# foilstudy -- An example of performance space creation and analysis

This repository accompanies the publication Hebdon, He, Battista, and Waldrop (in prep) 
and contains the necessary code to reproduce the performance-space generation figure therein. 

## Required Software

The following software is required for operating the code in this repository. Specific versions
are listed although other versions may also work. 

 - R version 4.5.1 (2025-06-13)
 - RStudio Version 2025.05.0+496 
 - R Packages (available on CRAN): `pracma`, `ggplot2`, `patchwork`, `spacefillr`

Note that the repository contains the results from the MATLAB simulations, so MATLAB is not required 
to generate surrogates or performance-space analyses. However, if you wish to reproduce the simulations 
themselves, MATLAB and XFOIL for MATLAB are required. 

 - MATLAB 2024a. 
 - XFOIL for MATLAB: https://www.mathworks.com/matlabcentral/fileexchange/50070-xfoil-for-matlab
Download and install by adding the installation folder to your MATLAB path. 

## Instructions for Reproducing Results

Please follow the steps below to reproduce the analyses. Note that several steps can be performed with either R or MATLAB, depending on your preference. Start by cloning the github repository or download the release and unarchive the code. If running in RStudio, open the `foilstudy.Rproj` in RStudio. You may run the `install_R_packages.R` script in `./src/r-scripts/` folder to install the required R packages. 

Next, run the following scripts in order: 

 1. Generate the parameter spaces using grid sampling, gPC sampling, and neural-network sampling.
    * MATLAB: open and run all lines `./doc/1-Generating_points_MATLAB.mlx` in MATLAB. 
    * R: open and run all lines `./doc/1-Generating_points.Rmd` in RStudio. (Note: NN sampling is not yet implemented in R!)
 2. Generate the airfoil files required to run in XFOIL. 
    * MATLAB: Open the `./doc/2-Creating_cambers_MATLAB.mlx` in MATLAB and run all lines.
    * R: Open the `./doc/2-Creating_cambers.Rmd` file in RStudio and run all lines, including the last code chunk which requires that a line be uncommented.  
 3. (MATLAB Only) Open the `./doc/3-xfoil_step.mlx` file MATLAB and run all sections. Note that this will take several hours to complete. Results are included in the project folder `./results/`, so it is not necessary to complete this step if only using R.
 4. (R Only) Open the `4-Visualizing_results.Rmd` in `doc/` and run all lines. This will provide information on the raw results of the simulations.
 5. (R Only) Open the `5-creating _surrogates.Rmd` in `doc/` and run all lines. This will provide space to generate the gPC and NN surrogates, including training the neural network model. 
 
__Note:__ each RMD should have the Knit directory to `project` to successfully knit the document and run code. 

Rerunning the code will not reproduce the results exactly because of some randomness in the process of point generation and XFOIL functions. These small differences will not affect the overall information presented in the final figure.
 
## More Information on Surrogate-construction Methods
 
There are code tutorials available for generalized polynomial chaos (R: `gPC_primer.Rmd`) and neural networks (R: `nn_primer.Rmd`). These use general examples to create surrogates and can be used as a code basis for your analysis. These tutorials, and the Creating Surrogates files, also include notes about using each method but should not be considered exhaustive on either method. Specialized expertise is required to use these methods to produce meaningful results for your system.  

