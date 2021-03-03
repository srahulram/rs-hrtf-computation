function exportFigure(FIG,FPATH)
%EXPORTFIGURE Export a figure.
%   EXPORTFIGURE(FIG,FPATH) exports the figure specified by the handle FIG
%   to the path specified by FPATH. The extension to the file name
%   determines the export format. Valid extensions are '.pdf' and '.eps'. 
%   If a valid extension is not provided, a '.eps' is exported by default.

%   =======================================================================
%   This file is part of the 3D3A Research Lite Toolbox.
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

narginchk(2,2)

validateattributes(FIG,{'matlab.ui.Figure'},{},'exportFigure','FIG',1);
validateattributes(FPATH,{'char'},{'scalartext'},'exportFigure','FPATH',2);
[FDIR,FNAME,FEXT] = fileparts(FPATH);

if exist(FDIR,'dir') ~= 7
    mkdir(FDIR);
end

switch lower(FEXT)
    case '.pdf'
        exportTag = '-dpdf';
        renameFlag = false;
    case '.eps'
        exportTag = '-dpsc'; % To prevent tight bounding box
        renameFlag = true; 
    otherwise
        warning(['Specified file name does not have a valid extension.',...
            ' Exporting to eps.'])
        FEXT = '.eps';
        exportTag = '-dpsc';
        renameFlag = true;
end

% Set additional figure properties for exporting
posVec = FIG.Position;
totalPlotWidth = posVec(3);
totalPlotHeight = posVec(4);
FIG.PaperUnits = 'centimeters';
FIG.PaperPositionMode = 'auto';
FIG.PaperSize = [totalPlotWidth+0.1,totalPlotHeight+0.1];

% Export figure
if renameFlag
    tempExportPath = fullfile(FDIR,[FNAME,'.ps']);
    print(FIG,tempExportPath,exportTag);
    
    exportPath = fullfile(FDIR,[FNAME,FEXT]);
    [status,~,~] = movefile(tempExportPath,exportPath,'f');
else
    exportPath = fullfile(FDIR,[FNAME,FEXT]);
    print(FIG,exportPath,exportTag);
    
    status = true;
end

if status
    fprintf('Figure exported to file:\n')
    fprintf('%s\n',rel2abs(exportPath));
end

end
