# Software Dependencies

Some scripts and functions in the RS-HRTF Computation software require additional toolboxes, APIs and software to work. If these toolboxes (see list below) are not already added to the MATLAB search path, perform **any one** of the following tasks:

1. If the toolboxes have already been downloaded manually, either place them here in the "dependencies" folder or anywhere in the default MATLAB folder (type "userpath" in MATLAB's Command Window to get the full path to this folder on your device). Alternatively, if stored elsewhere, manually add the location to the MATLAB search path. Finally, in MATLAB, navigate to the RS-HRTF_Computation > MATLAB folder and call the function "rshrtf_start.m" from MATLAB's Command Window.

2. Connect to the internet, navigate to RS-HRTF_Computation > MATLAB from within MATLAB itself and call the function "rshrtf_start.m" from MATLAB's Command Window. If the function does not find the required toolboxes in the MATLAB search path, the default MATLAB folder, or here in the "dependencies" folder, an attempt will be made to download them from the web and store them here in the "dependencies" folder.

## List of Dependencies (Toolboxes, APIs, and Software)

1. 3D3A MATLAB Toolbox - https://github.com/PrincetonUniversity/3D3A-MATLAB-Toolbox/archive/v0.1.0.zip
2. Hatchfill2 - https://www.mathworks.com/matlabcentral/fileexchange/53593-hatchfill2; already included here in "dependencies" folder.
3. 3D3A Research Lite Toolbox - direct download not available; already included here in "dependencies" folder.

Note: The 3D3A MATLAB Toolbox requires additional dependencies. For more information, refer to the README file in that toolbox.