# RS-HRTF Computation MATLAB Software

Collection of MATLAB scripts and functions for computing an RS-HRTF and reproducing the results in the publication:  

R. Sridhar and E. Y. Choueiri, "Optimal Series Truncation of the Rigid-Sphere Head-Related Transfer Function for Accurate Binaural Perception," J. Audio Eng. Soc., 2021 (to appear)

**Note:** This software has only been tested on a full version of MATLAB R2018a. It uses multiple pre-packaged MATLAB toolboxes such as the Signal Processing and Curve Fitting toolboxes. Furthermore, some of the scripts and functions require additional, third-party toolboxes/software to work. For more information on the latter, see the README file in the "dependencies" folder.

## Usage Guide

**Note:** Please go through the README file in the "dependencies" folder before proceeding.

###### Starting the Software

In MATLAB, navigate to the directory containing "rshrtf_start.m" and run the command rshrtf_start from the Command Window. This should launch any dependencies and add required paths to MATLAB's default search path.

###### Running Scripts

In MATLAB, navigate to the "scripts" folder. First read the descriptions and usage instructions at the top of each script (above the boilerplate). **Do not run an entire script all at once**. Instead, run each cell in the script one at a time in top-down order, following any additional instructions given in each cell.