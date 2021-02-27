function varargout = formatLinePlot(F,P,N,varargin)
%FORMATLINEPLOT Apply formatting to line plot.
%   F = FORMATLINEPLOT(F,P,N) takes a figure object, F, consisting of a 2D 
%   line plot (generated using the 'plot' command), P, with N unique data 
%   vectors (i.e., N dependent variables being plotted) and formats it for
%   publication. The formatted figure is returned in F. Note that P must be
%   a plot object generated using the 'plot' command.
%
%   F = FORMATLINEPLOT(___,Name,Value) allows specification of one or 
%   more *optional* Name,Value pairs that control various aspects of plot
%   formatting. The following pairs may be specified:
%       'BXLabel', BXLAB - BXLAB must be a character vector and corresponds
%       to the bottom x-axis label. LaTeX allowed.
%
%       'LYLabel', LYLAB - LYLAB must be a character vector and corresponds
%       to the left y-axis label. LaTeX allowed.
%
%       'TXLabel', TXLAB - TXLAB must be a character vector and corresponds
%       to the top x-axis label. LaTeX allowed.
%
%       'RYLabel', RYLAB - RYLAB must be a character vector and corresponds
%       to the right y-axis label. LaTeX allowed.
%
%       'BXLim', XLIMSB - XLIMSB must be a two-element vector specifying 
%       the plot limits for the bottom x-axis.
%
%       'LYLim', YLIMSL - YLIMSL must be a two-element vector specifying 
%       the plot limits for the left y-axis.
%
%       'TXLim', XLIMST - XLIMST must be a two-element vector specifying 
%       the plot limits for the top x-axis. By default, this is the same
%       as XLIMSB.
%
%       'RYLim', YLIMSR - YLIMSR must be a two-element vector specifying 
%       the plot limits for the right y-axis. By default, this is the same
%       as YLIMSL.
%
%       'BXScale', XSCALEB - XSCALEB can be 'linear' (default) or 'log'.
%
%       'LYScale', YSCALEL - YSCALEL can be 'linear' (default) or 'log'.
%
%       'TXScale', XSCALET - XSCALET can be 'linear' (default) or 'log'.
%
%       'RYScale', YSCALER - YSCALER can be 'linear' (default) or 'log'.
%
%       'BXTick', XTICKVECB - XTICKVECB must be a vector specifying the
%       values on the bottom x-axis where tick marks should be drawn. By 
%       default, 5 equally spaced ticks are generated. To prevent ticks
%       from being drawn, specify XTICKVECB as an empty vector.
%
%       'LYTick', YTICKVECL - YTICKVECL must be a vector specifying the
%       values on the left y-axis where tick marks should be drawn. By 
%       default, 5 equally spaced ticks are generated. To prevent ticks
%       from being drawn, specify YTICKVECL as an empty vector.
%
%       'TXTick', XTICKVECT - XTICKVECT must be a vector specifying the
%       values on the top x-axis where tick marks should be drawn. By 
%       default, 5 equally spaced ticks are generated. To prevent ticks
%       from being drawn, specify XTICKVECT as an empty vector.
%
%       'RYTick', YTICKVECR - YTICKVECR must be a vector specifying the
%       values on the right y-axis where tick marks should be drawn. By 
%       default, 5 equally spaced ticks are generated. To prevent ticks
%       from being drawn, specify YTICKVECR as an empty vector.
%
%       'BXTickLabel', XTICKLBLSB - XTICKLBLSB must be a cell array of the
%       same length as XTICKVECB where each element must be a character
%       vector that is the tick label corresponding to the tick value in
%       XTICKVECB. To prevent labels from being generated specify 
%       XTICKLBLSB as an empty cell array.
%
%       'LYTickLabel', YTICKLBLSL - YTICKLBLSL must be a cell array of the
%       same length as YTICKVECL where each element must be a character
%       vector that is the tick label corresponding to the tick value in
%       YTICKVECL. To prevent labels from being generated specify 
%       YTICKLBLSL as an empty cell array.
%
%       'TXTickLabel', XTICKLBLST - XTICKLBLST must be a cell array of the
%       same length as XTICKVECT where each element must be a character
%       vector that is the tick label corresponding to the tick value in
%       XTICKVECT. To prevent labels from being generated specify 
%       XTICKLBLST as an empty cell array.
%
%       'RYTickLabel', YTICKLBLSR - YTICKLBLSR must be a cell array of the
%       same length as YTICKVECR where each element must be a character
%       vector that is the tick label corresponding to the tick value in
%       YTICKVECR. To prevent labels from being generated specify 
%       YTICKLBLSR as an empty cell array.
%
%       'BXMinorTick', xMinTickBFlag - xMinTickBFlag can be 'on' or 'off'
%       (default) and specifies whether or not to include minor ticks in
%       the bottom x axis.
%
%       'BXMinorTickVal', XMINTICKVECB - XMINTICKVECB must be a vector 
%       specifying the values on the bottom x-axis where minor tick marks 
%       should be drawn. To prevent minor ticks from being drawn, specify 
%       xMinTickBFlag (above) as 'off'.
%
%       'TXMinorTick', xMinTickTFlag - xMinTickTFlag can be 'on' or 'off'
%       (default) and specifies whether or not to include minor ticks in
%       the top x axis.
%
%       'TXMinorTickVal', XMINTICKVECT - XMINTICKVECT must be a vector 
%       specifying the values on the top x-axis where minor tick marks 
%       should be drawn. To prevent minor ticks from being drawn, specify 
%       xMinTickTFlag (above) as 'off'.
%
%       'LYMinorTick', yMinTickLFlag - yMinTickLFlag can be 'on' or 'off'
%       (default) and specifies whether or not to include minor ticks in
%       the left y axis.
%
%       'LYMinorTickVal', YMINTICKVECL - YMINTICKVECL must be a vector 
%       specifying the values on the left y-axis where minor tick marks 
%       should be drawn. To prevent minor ticks from being drawn, specify 
%       yMinTickLFlag (above) as 'off'.
%
%       'RYMinorTick', yMinTickRFlag - yMinTickRFlag can be 'on' or 'off'
%       (default) and specifies whether or not to include minor ticks in
%       the right y axis.
%
%       'RYMinorTickVal', YMINTICKVECR - YMINTICKVECR must be a vector 
%       specifying the values on the right y-axis where minor tick marks 
%       should be drawn. To prevent minor ticks from being drawn, specify 
%       yMinTickRFlag (above) as 'off'.
%
%       'XGrid', xGridFlag - xGridFlag can be 'on' (default) or 'off' and 
%       specifies whether or not to include grid lines for the bottom x-
%       axis.
%
%       'XMinorGrid', xMinorGridFlag - xMinorGridFlag can be 'on' or 'off' 
%       (default) and specifies whether or not to include minor grid lines 
%       for the bottom x-axis.
%
%       'YGrid', yGridFlag - yGridFlag can be 'on' (default) or 'off' and 
%       specifies whether or not to include grid lines for the left y-
%       axis.
%
%       'YMinorGrid', yMinorGridFlag - yMinorGridFlag can be 'on' or 'off' 
%       (default) and specifies whether or not to include minor grid lines 
%       for the left y-axis.
%
%       'includeColor', colorFlag - colorFlag can be true (default) or
%       false and specifies whether the plot should be colored or
%       grayscale. Colors are automatically chosen with a different color
%       being applied to each line.
%
%       'colorList', lineColorMat - lineColorMat should be a matrix 
%       specifying line colors in RGB format.
%
%       'lineStyleList', lineCell - lineCell should be a cell array
%       specifying the line styles in the order they should be plotted. If
%       specified, the number of elements in the input cell array should
%       correspond to the number of individual data vectors being plotted. 
%       To prevent a line being drawn for a given set of data, specify the 
%       line style corresponding to the data as 'none'. By default, all 
%       data are plotted as solid lines.
%
%       'markerStyleList', markerCell - markerCell should be a cell array
%       specifying the marker styles in the order they should be plotted.
%       If specified, the number of elements in the input cell array should
%       correspond to the number of individual data vectors being plotted. 
%       To prevent a marker being drawn for a given set of data, specify 
%       the marker style corresponding to the data as 'none'. By default, 
%       all data are plotted without markers.
%
%       'markerFillFlag', markerFill - markerFill can be true (to generate 
%       filled markers) or false (to generate empty markers). The default
%       value is false.
%
%       'plotAspect', AR - AR should be the aspect ratio of the plot. The
%       default value is 4/3.
%
%   [F,P,A1,A2] = FORMATLINEPLOT(___) optionally returns the plot object, 
%   P, and axes objects, A1 (corresponding to primary axes, namely the 
%   bottom x-axis and left y-axis), and A2 (corresponding to secondary 
%   axes, namely the top x-axis and right y-axis).
%
%   See also PLOT.

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

% Check input count
narginchk(3,inf);

% Validate required inputs
validateattributes(F,{'matlab.ui.Figure'},{},'formatLinePlot','F',1);
validateattributes(P,{'matlab.graphics.chart.primitive.Line',...
    'matlab.graphics.chart.primitive.Stair',...
    'matlab.graphics.chart.primitive.ErrorBar'},{},'formatLinePlot','P',2);
validateattributes(N,{'numeric'},{'scalar','finite','positive'},...
    'formatLinePlot','N',3);

% Get current axes
A1 = get(F,'CurrentAxes');

% Validate optional inputs
indx = find(strcmpi(varargin,'BXLabel'),1);
if isempty(indx)
    BXLAB = '';
else
    BXLAB = varargin{indx+1};
    validateattributes(BXLAB,{'char'},{'scalartext'},'formatLinePlot',...
        'BXLAB');
end

indx = find(strcmpi(varargin,'LYLabel'),1);
if isempty(indx)
    LYLAB = '';
else
    LYLAB = varargin{indx+1};
    validateattributes(LYLAB,{'char'},{'scalartext'},'formatLinePlot',...
        'LYLAB');
end

indx = find(strcmpi(varargin,'TXLabel'),1);
if isempty(indx)
    TXLAB = '';
else
    TXLAB = varargin{indx+1};
    validateattributes(TXLAB,{'char'},{'scalartext'},'formatLinePlot',...
        'TXLAB');
end

indx = find(strcmpi(varargin,'RYLabel'),1);
if isempty(indx)
    RYLAB = '';
else
    RYLAB = varargin{indx+1};
    validateattributes(RYLAB,{'char'},{'scalartext'},'formatLinePlot',...
        'RYLAB');
end

indx = find(strcmpi(varargin,'BXLim'),1);
if isempty(indx)
    XLIMSB = A1.XLim;
else
    XLIMSB = sort(varargin{indx+1});
    validateattributes(XLIMSB,{'numeric'},{'vector','finite','real'},...
        'formatLinePlot','XLIMSB');
end

indx = find(strcmpi(varargin,'LYLim'),1);
if isempty(indx)
    YLIMSL = A1.YLim;
else
    YLIMSL = sort(varargin{indx+1});
    validateattributes(YLIMSL,{'numeric'},{'vector','finite','real'},...
        'formatLinePlot','YLIMSL');
end

indx = find(strcmpi(varargin,'TXLim'),1);
if isempty(indx)
    XLIMST = XLIMSB;
else
    XLIMST = sort(varargin{indx+1});
    validateattributes(XLIMST,{'numeric'},{'vector','finite','real'},...
        'formatLinePlot','XLIMST');
end

indx = find(strcmpi(varargin,'RYLim'),1);
if isempty(indx)
    YLIMSR = YLIMSL;
else
    YLIMSR = sort(varargin{indx+1});
    validateattributes(YLIMSR,{'numeric'},{'vector','finite','real'},...
        'formatLinePlot','YLIMSR');
end

indx = find(strcmpi(varargin,'BXScale'),1);
if isempty(indx)
    XSCALEB = 'linear';
else
    XSCALEB = varargin{indx+1};
    validateattributes(XSCALEB,{'char'},{'scalartext'},'formatLinePlot',...
        'XSCALEB');
end

indx = find(strcmpi(varargin,'LYScale'),1);
if isempty(indx)
    YSCALEL = 'linear';
else
    YSCALEL = varargin{indx+1};
    validateattributes(YSCALEL,{'char'},{'scalartext'},'formatLinePlot',...
        'YSCALEL');
end

indx = find(strcmpi(varargin,'TXScale'),1);
if isempty(indx)
    XSCALET = 'linear';
else
    XSCALET = varargin{indx+1};
    validateattributes(XSCALET,{'char'},{'scalartext'},'formatLinePlot',...
        'XSCALET');
end

indx = find(strcmpi(varargin,'RYScale'),1);
if isempty(indx)
    YSCALER = 'linear';
else
    YSCALER = varargin{indx+1};
    validateattributes(YSCALER,{'char'},{'scalartext'},'formatLinePlot',...
        'YSCALER');
end

indx = find(strcmpi(varargin,'BXTick'),1);
if isempty(indx)
    XTICKVECB = A1.XTick;
    BXTickFlag = false;
else
    XTICKVECB = varargin{indx+1};
    validateattributes(XTICKVECB,{'numeric'},{'2d','finite','real'},...
        'formatLinePlot','XTICKVECB');
    BXTickFlag = true;
end

indx = find(strcmpi(varargin,'LYTick'),1);
if isempty(indx)
    YTICKVECL = A1.YTick;
    LYTickFlag = false;
else
    YTICKVECL = varargin{indx+1};
    validateattributes(YTICKVECL,{'numeric'},{'2d','finite','real'},...
        'formatLinePlot','YTICKVECL');
    LYTickFlag = true;
end

indx = find(strcmpi(varargin,'TXTick'),1);
if isempty(indx)
    XTICKVECT = XTICKVECB;
    TXTickFlag = false;
else
    XTICKVECT = varargin{indx+1};
    validateattributes(XTICKVECT,{'numeric'},{'2d','finite','real'},...
        'formatLinePlot','XTICKVECT');
    TXTickFlag = true;
end

indx = find(strcmpi(varargin,'RYTick'),1);
if isempty(indx)
    YTICKVECR = YTICKVECL;
    RYTickFlag = false;
else
    YTICKVECR = varargin{indx+1};
    validateattributes(YTICKVECR,{'numeric'},{'2d','finite','real'},...
        'formatLinePlot','YTICKVECR');
    RYTickFlag = true;
end

indx = find(strcmpi(varargin,'BXTickLabel'),1);
if isempty(indx)
    numXTicks = length(XTICKVECB);
    if numXTicks > 0
        if BXTickFlag
            XTICKLBLSB = cell(1,numXTicks);
            for ii = 1:numXTicks
                XTICKLBLSB{1,ii} = num2str(XTICKVECB(ii));
            end
        else
            XTICKLBLSB = A1.XTickLabel;
        end
    else
        XTICKLBLSB = {};
    end
else
    XTICKLBLSB = varargin{indx+1};
    validateattributes(XTICKLBLSB,{'cell'},{'2d'},'formatLinePlot',...
        'XTICKLBLSB');
    if ~isempty(XTICKLBLSB) && (length(XTICKLBLSB) ~= length(XTICKVECB))
        error('XTICKLBLSB must be empty or the same length as XTICKVECB')
    end
end

indx = find(strcmpi(varargin,'LYTickLabel'),1);
if isempty(indx)
    numYTicks = length(YTICKVECL);
    if numYTicks > 0
        if LYTickFlag
            YTICKLBLSL = cell(1,numYTicks);
            for ii = 1:numYTicks
                YTICKLBLSL{1,ii} = num2str(YTICKVECL(ii));
            end
        else
            YTICKLBLSL = A1.YTickLabel;
        end
    else
        YTICKLBLSL = {};
    end
else
    YTICKLBLSL = varargin{indx+1};
    validateattributes(YTICKLBLSL,{'cell'},{'2d'},'formatLinePlot',...
        'YTICKLBLSL');
    if ~isempty(YTICKLBLSL) && (length(YTICKLBLSL) ~= length(YTICKVECL))
        error('YTICKLBLSL must be empty or the same length as YTICKVECL')
    end
end

indx = find(strcmpi(varargin,'TXTickLabel'),1);
if isempty(indx)
    numXTicks = length(XTICKVECT);
    if numXTicks > 0
        if TXTickFlag
            XTICKLBLST = cell(1,numXTicks);
            for ii = 1:numXTicks
                XTICKLBLST{1,ii} = num2str(XTICKVECT(ii));
            end
        else
            XTICKLBLST = {};
        end
    else
        XTICKLBLST = {};
    end
else
    XTICKLBLST = varargin{indx+1};
    validateattributes(XTICKLBLST,{'cell'},{'2d'},'formatLinePlot',...
        'XTICKLBLST');
    if ~isempty(XTICKLBLST) && (length(XTICKLBLST) ~= length(XTICKVECT))
        error('XTICKLBLST must be empty or the same length as XTICKVECT')
    end
end

indx = find(strcmpi(varargin,'RYTickLabel'),1);
if isempty(indx)
    numYTicks = length(YTICKVECR);
    if numYTicks > 0
        if RYTickFlag
            YTICKLBLSR = cell(1,numYTicks);
            for ii = 1:numYTicks
                YTICKLBLSR{1,ii} = num2str(YTICKVECR(ii));
            end
        else
            YTICKLBLSR = {};
        end
    else
        YTICKLBLSR = {};
    end
else
    YTICKLBLSR = varargin{indx+1};
    validateattributes(YTICKLBLSR,{'cell'},{'2d'},'formatLinePlot',...
        'YTICKLBLSR');
    if ~isempty(YTICKLBLSR) && (length(YTICKLBLSR) ~= length(YTICKVECR))
        error('YTICKLBLSR must be empty or the same length as YTICKVECR')
    end
end

indx = find(strcmpi(varargin,'BXMinorTick'),1);
if isempty(indx)
    xMinTickBFlag = 'off';
else
    xMinTickBFlag = varargin{indx+1};
    validateattributes(xMinTickBFlag,{'char'},{'scalartext'},...
        'formatLinePlot','xMinTickBFlag');
end

indx = find(strcmpi(varargin,'BXMinorTickVal'),1);
if isempty(indx)
    BXMinorTickValFlag = false;
else
    XMINTICKVECB = varargin{indx+1};
    validateattributes(XMINTICKVECB,{'numeric'},{'2d','finite','real'},...
        'formatLinePlot','XMINTICKVECB');
    BXMinorTickValFlag = true;
end

indx = find(strcmpi(varargin,'TXMinorTick'),1);
if isempty(indx)
    xMinTickTFlag = 'off';
else
    xMinTickTFlag = varargin{indx+1};
    validateattributes(xMinTickTFlag,{'char'},{'scalartext'},...
        'formatLinePlot','xMinTickTFlag');
end

indx = find(strcmpi(varargin,'TXMinorTickVal'),1);
if isempty(indx)
    TXMinorTickValFlag = false;
else
    XMINTICKVECT = varargin{indx+1};
    validateattributes(XMINTICKVECT,{'numeric'},{'2d','finite','real'},...
        'formatLinePlot','XMINTICKVECT');
    TXMinorTickValFlag = true;
end

indx = find(strcmpi(varargin,'LYMinorTick'),1);
if isempty(indx)
    yMinTickLFlag = 'off';
else
    yMinTickLFlag = varargin{indx+1};
    validateattributes(yMinTickLFlag,{'char'},{'scalartext'},...
        'formatLinePlot','yMinTickLFlag');
end

indx = find(strcmpi(varargin,'LYMinorTickVal'),1);
if isempty(indx)
    LYMinorTickValFlag = false;
else
    YMINTICKVECL = varargin{indx+1};
    validateattributes(YMINTICKVECL,{'numeric'},{'2d','finite','real'},...
        'formatLinePlot','YMINTICKVECL');
    LYMinorTickValFlag = true;
end

indx = find(strcmpi(varargin,'RYMinorTick'),1);
if isempty(indx)
    yMinTickRFlag = 'off';
else
    yMinTickRFlag = varargin{indx+1};
    validateattributes(yMinTickRFlag,{'char'},{'scalartext'},...
        'formatLinePlot','yMinTickRFlag');
end

indx = find(strcmpi(varargin,'RYMinorTickVal'),1);
if isempty(indx)
    RYMinorTickValFlag = false;
else
    YMINTICKVECR = varargin{indx+1};
    validateattributes(YMINTICKVECR,{'numeric'},{'2d','finite','real'},...
        'formatLinePlot','YMINTICKVECR');
    RYMinorTickValFlag = true;
end

indx = find(strcmpi(varargin,'XGrid'),1);
if isempty(indx)
    xGridFlag = 'on';
else
    xGridFlag = varargin{indx+1};
    validateattributes(xGridFlag,{'char'},{'scalartext'},...
        'formatLinePlot','xGridFlag');
end

indx = find(strcmpi(varargin,'XMinorGrid'),1);
if isempty(indx)
    xMinorGridFlag = 'off';
else
    xMinorGridFlag = varargin{indx+1};
    validateattributes(xMinorGridFlag,{'char'},{'scalartext'},...
        'formatLinePlot','xMinorGridFlag');
end

indx = find(strcmpi(varargin,'YGrid'),1);
if isempty(indx)
    yGridFlag = 'on';
else
    yGridFlag = varargin{indx+1};
    validateattributes(yGridFlag,{'char'},{'scalartext'},...
        'formatLinePlot','yGridFlag');
end

indx = find(strcmpi(varargin,'YMinorGrid'),1);
if isempty(indx)
    yMinorGridFlag = 'off';
else
    yMinorGridFlag = varargin{indx+1};
    validateattributes(yMinorGridFlag,{'char'},{'scalartext'},...
        'formatLinePlot','yMinorGridFlag');
end

indx = find(strcmpi(varargin,'includeColor'),1);
if isempty(indx)
    colorFlag = true;
else
    colorFlag = varargin{indx+1};
    validateattributes(colorFlag,{'logical'},{'scalar'},...
        'formatLinePlot','colorFlag');
end

indx = find(strcmpi(varargin,'colorList'),1);
if isempty(indx)
    if colorFlag
        lineColorMat = (1/255)*[0,0,0; % black
            255,0,0; % red
            0,128,0; % green
            0,0,255; % blue
            241,152,46; % orange
            185,102,245; % purple
            4,210,201; % turqoise blue
            0,255,0]; % bright green
    else
        if N == 1
            lineColorMat = [0,0,0];
        else
            lineColorMat = repmat(linspace(0,0.6,N).',1,3);
        end
    end
else
    lineColorMat = varargin{indx+1};
    validateattributes(lineColorMat,{'numeric'},{'2d','nonnegative',...
        'nonempty','real','<=',1,'size',[NaN,3]},'formatLinePlot',...
        'lineColorMat');
end

indx = find(strcmpi(varargin,'lineStyleList'),1);
if isempty(indx)
    lineCell = {'-'};
else
    lineCell = varargin{indx+1};
    validateattributes(lineCell,{'cell'},{'vector'},'formatLinePlot',...
        'lineStyleList');
end

indx = find(strcmpi(varargin,'markerStyleList'),1);
if isempty(indx)
    markerCell = {'none'};
else
    markerCell = varargin{indx+1};
    validateattributes(markerCell,{'cell'},{'vector'},'formatLinePlot',...
        'markerStyleList');
end

indx = find(strcmpi(varargin,'markerFillFlag'),1);
if isempty(indx)
    markerFill = false;
else
    markerFill = varargin{indx+1};
    validateattributes(markerFill,{'logical'},{'scalar'},...
        'formatLinePlot','markerFill');
end

indx = find(strcmpi(varargin,'plotAspect'),1);
if isempty(indx)
    AR = 4/3;
else
    AR = varargin{indx+1};
    validateattributes(AR,{'numeric'},{'scalar','finite','positive'},...
        'formatLinePlot','AR');
end

% Desired plot box dimension (excluding tick and axes labels), in cm
plotBoxWidthCm = 12;
% Uniform padding for tick and axes labels, in cm
plotBoxPadCm = 3;
% Other formatting parameters
fontSize = 22;
fontName = 'Times';
plotBoxLineWidth = 0.75;
plotLineWidth = 1.5;

% General Properties
set(0,'DefaultTextInterpreter','latex')

% Line Properties
numColors = size(lineColorMat,1);
numLineStyles = numel(lineCell);
for ii = 1:N
    P(ii).LineWidth = plotLineWidth;
    colorIndx = mod(ii,numColors);
    colorIndx(colorIndx == 0) = numColors;
    P(ii).Color = lineColorMat(colorIndx,:);
    lineStyleIndx = mod(ii,numLineStyles);
    lineStyleIndx(lineStyleIndx == 0) = numLineStyles;
    P(ii).LineStyle = lineCell{1,lineStyleIndx};
end

% Marker Properties
numMarkerStyles = numel(markerCell);
for ii = 1:N
    P(ii).LineWidth = plotLineWidth;
    colorIndx = mod(ii,numColors);
    colorIndx(colorIndx == 0) = numColors;
    P(ii).Color = lineColorMat(colorIndx,:);
    markerStyleIndx = mod(ii,numMarkerStyles);
    markerStyleIndx(markerStyleIndx == 0) = numMarkerStyles;
    P(ii).Marker = markerCell{1,markerStyleIndx};
    P(ii).MarkerSize = 8;
    P(ii).MarkerEdgeColor = lineColorMat(colorIndx,:);
    if markerFill
        P(ii).MarkerFaceColor = lineColorMat(colorIndx,:);
    else
        P(ii).MarkerFaceColor = 'none';
    end
end

% Axes Properties
A1.Units = 'centimeters';
A1.FontSize = fontSize;
A1.FontName = fontName;
A1.TickLabelInterpreter = 'latex';

% Axes Ruler Properties
A1.XLim = XLIMSB;
A1.YLim = YLIMSL;
A1.XColor = 'black';
A1.YColor = 'black';
A1.XScale = XSCALEB;
A1.YScale = YSCALEL;

% Apply X Ticks and Labels
A1.XTick = XTICKVECB;
A1.XTickLabel = XTICKLBLSB;
A1.XMinorTick = xMinTickBFlag;
if strcmpi(xMinTickBFlag,'on')
    if BXMinorTickValFlag
        A1.XAxis.MinorTickValues = XMINTICKVECB;
    end
end

% Apply Y Ticks and Labels
A1.YTick = YTICKVECL;
A1.YTickLabel = YTICKLBLSL;
A1.YMinorTick = yMinTickLFlag;
if strcmpi(yMinTickLFlag,'on')
    if LYMinorTickValFlag
        A1.YAxis.MinorTickValues = YMINTICKVECL;
    end
end

% Grid Properties
A1.XGrid = xGridFlag;
A1.YGrid = yGridFlag;
A1.XMinorGrid = xMinorGridFlag;
A1.YMinorGrid = yMinorGridFlag;

% Axes and Title Labels
A1.XLabel.String = BXLAB;
A1.XLabel.FontSize = fontSize;
A1.XLabel.FontName = fontName;
A1.YLabel.String = LYLAB;
A1.YLabel.FontSize = fontSize;
A1.YLabel.FontName = fontName;

% Plot Box Properties
A1.Color = 'none';
A1.LineWidth = plotBoxLineWidth;
A1.XRuler.Axle.LineWidth = plotBoxLineWidth;
A1.YRuler.Axle.LineWidth = plotBoxLineWidth;
A1.TickDir = 'in';
A1.TickLength = [0.014 0.035];
A1.Box = 'off';

% Apply formatting to secondary axes
A1_pos = A1.Position;
A2 = axes('Units','centimeters','Position',A1_pos,'XAxisLocation','top',...
    'YAxisLocation','right','Color','none');

% Common Axes Properties
A2.FontSize = fontSize;
A2.FontName = fontName;
A2.TickLabelInterpreter = 'latex';

% Axes Ruler Properties
A2.XLim = XLIMST;
A2.YLim = YLIMSR;
A2.XColor = 'black';
A2.YColor = 'black';
A2.XScale = XSCALET;
A2.YScale = YSCALER;

% Apply X Ticks and Labels
A2.XTick = XTICKVECT;
A2.XTickLabel = XTICKLBLST;
A2.XMinorTick = xMinTickTFlag;
if strcmpi(xMinTickTFlag,'on')
    if TXMinorTickValFlag
        A2.XAxis.MinorTickValues = XMINTICKVECT;
    end
end

% Apply Y Ticks and Labels
A2.YTick = YTICKVECR;
A2.YTickLabel = YTICKLBLSR;
A2.YMinorTick = yMinTickRFlag;
if strcmpi(yMinTickRFlag,'on')
    if RYMinorTickValFlag
        A2.YAxis.MinorTickValues = YMINTICKVECR;
    end
end

% Axes Labels
A2.XLabel.String = TXLAB;
A2.XLabel.FontSize = fontSize;
A2.XLabel.FontName = fontName;
A2.YLabel.String = RYLAB;
A2.YLabel.FontSize = fontSize;
A2.YLabel.FontName = fontName;

% Plot Box Properties
A2.Color = 'none';
A2.LineWidth = plotBoxLineWidth;
A2.XRuler.Axle.LineWidth = plotBoxLineWidth;
A2.YRuler.Axle.LineWidth = plotBoxLineWidth;
A2.TickDir = 'in';
A2.TickLength = [0.014 0.035];
A2.Box = 'off';

% Calculate plot resize parameters
axTightInset = A1.TightInset;
ax2TightInset = A2.TightInset;
maxTightInset = max([axTightInset;ax2TightInset]);
if any(maxTightInset > plotBoxPadCm)
    warning(['Margins may be too small to accommodate all tick and ',...
        'axes labels.'])
end

% Set aspect ratio
plotBoxHeightCm = plotBoxWidthCm*(1/AR);
% Resize outer box
totalPlotWidth = plotBoxWidthCm+(2*plotBoxPadCm);
totalPlotHeight = plotBoxHeightCm+(2*plotBoxPadCm);
A1.OuterPosition = [0,0,totalPlotWidth,totalPlotHeight];
% Resize plot box
A1.Position = [plotBoxPadCm,plotBoxPadCm,plotBoxWidthCm,plotBoxHeightCm];
A2.Position = [plotBoxPadCm,plotBoxPadCm,plotBoxWidthCm,plotBoxHeightCm];

% Resize figure
F.Units = 'centimeters';
figPos = F.Position;
F.Position = [figPos(1),figPos(2),totalPlotWidth,totalPlotHeight];

% Return outputs
switch nargout
    case 1
        varargout = {F};
    case 4
        varargout = {F,P,A1,A2};
    otherwise
        error('formatLinePlot expects request for 1 or 4 outputs.')
end

% Reset General Properties
set(0,'DefaultTextInterpreter','tex')

end
