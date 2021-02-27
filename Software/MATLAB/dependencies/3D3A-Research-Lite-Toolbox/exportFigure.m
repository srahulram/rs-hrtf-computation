function exportFigure(FIG,FPATH)
%EXPORTFIGURE Export a figure to a PDF.
%   EXPORTFIGURE(FIG,FPATH) exports the figure specified by the handle FIG
%   to a PDF file with path specified by FPATH.

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

if ~strcmpi(FEXT,'.pdf')
    warning(['Specified file name does not have .pdf extension.',...
        ' Exporting to pdf.'])
    FEXT = '.pdf';
end

% Set additional figure properties for exporting
posVec = FIG.Position;
totalPlotWidth = posVec(3);
totalPlotHeight = posVec(4);
FIG.PaperUnits = 'centimeters';
FIG.PaperPositionMode = 'auto';
FIG.PaperSize = [totalPlotWidth+0.1,totalPlotHeight+0.1];

% Export figure to pdf
print(FIG,fullfile(FDIR,[FNAME,FEXT]),'-dpdf');

fprintf('Figure exported to file:\n')
fprintf('%s\n',rel2abs(fullfile(FDIR,[FNAME,FEXT])));

end
