function rshrtf_start()
%RSHRTF_START Start RS-HRTF Computation software
%   RSHRTF_START adds to the MATLAB search path, folders and additional 
%   toolboxes required to run scripts in the RS-HRTF Computation software
%   package.

%   =======================================================================
%   This file is part of the RS-HRTF Computation software.
%   
%   Rahulram Sridhar <rahulram@princeton.edu>
%   3D Audio and Applied Acoustics (3D3A) Laboratory
%   Princeton University, Princeton, New Jersey 08544, USA
%   
%   MIT License
%   
%   Copyright (c) 2021 Rahulram Sridhar
%   
%   Permission is hereby granted, free of charge, to any person obtaining a
%   copy of this software and associated documentation files (the 
%   "Software"), to deal in the Software without restriction, including 
%   without limitation the rights to use, copy, modify, merge, publish, 
%   distribute, sublicense, and/or sell copies of the Software, and to 
%   permit persons to whom the Software is furnished to do so, subject to 
%   the following conditions:
%   
%   The above copyright notice and this permission notice shall be included
%   in all copies or substantial portions of the Software.
%   
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
%   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
%   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
%   CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
%   TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
%   SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%   =======================================================================

% Check input count
narginchk(0,0)

disp('Starting RS-HRTF Computation software...')
[parentDir,~,~] = fileparts(which('rshrtf_start'));
depDir = fullfile(parentDir,'dependencies');

% Try to find the 3D3A MATLAB Toolbox
disp('Looking for the 3D3A MATLAB Toolbox...')
found3D3A = false;
if exist('start3D3AMATLABToolbox','file') == 2
    found3D3A = true;
else
    dir3D3A1 = dir(fullfile(userpath,'**','start3D3AMATLABToolbox.m'));
    dir3D3A2 = dir(fullfile(depDir,'**','start3D3AMATLABToolbox.m'));
    if ~isempty(dir3D3A1)
        addpath(dir3D3A1(1).folder)
        found3D3A = true;
    elseif ~isempty(dir3D3A2)
        addpath(dir3D3A2(1).folder)
        found3D3A = true;
    else
        disp('Could not find local copy of the 3D3A MATLAB Toolbox.')
        disp('Attempting to download from web...')
        url = ['https://github.com/PrincetonUniversity/',...
            '3D3A-MATLAB-Toolbox/archive/v0.1.0.zip'];
        try
            unzip(url,fullfile(depDir,'3D3A-MATLAB-Toolbox'));
            found3D3A = true;
        catch
            disp('Unable to download/unzip the 3D3A MATLAB Toolbox.')
        end
        
        if found3D3A
            dir3D3A = dir(fullfile(depDir,'**',...
                'start3D3AMATLABToolbox.m'));
            if ~isempty(dir3D3A)
                addpath(dir3D3A(1).folder)
                disp('3D3A MATLAB Toolbox downloaded successfully.')
            end
        end
    end
end

if found3D3A
    disp('3D3A MATLAB Toolbox found.')
    start3D3AMATLABToolbox
else
    warning(['Could not find the 3D3A MATLAB Toolbox. Some parts of',...
        ' the software will not work.'])
end

addpath(genpath(fullfile(parentDir)))
if exist(fullfile(parentDir,'dependencies'),'dir') == 7
    rmpath(fullfile(parentDir,'dependencies'));
end

disp('RS-HRTF Computation software initialization complete.')

end
