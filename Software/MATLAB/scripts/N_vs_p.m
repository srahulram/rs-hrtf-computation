%% Truncation order versus numerical precision.
%
% In this script, we see how truncation order, N, varies with numerical
% precision, p.
%
% Note: Calculations done for left ear HRTF only.

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

%% Run initialization script

rshrtf_init

%% Define parameters

% Head parameters
head.a = 0.09;
head.eL = [90,0];

% Source parameters
rho = inf;
source.r = rho*(head.a);
source.azVec = 90:5:270;
source.elVec = 0;

% Computation parameters
muVec = [0.1;5;10;20;30];
numMus = length(muVec);
dsp.fVec = muVec*getSoundSpeed()/(2*pi*head.a);
pVec = (0:20).';
numPs = length(pVec);

%% Compute N for various p using the method by Sridhar and Choueiri (2020)

fprintf('Computing using ''sridharchoueiri2020'' method...\n');

% Compute N for various p values
N_sc = cell(numPs,1);
for ii = 1:numPs
    dsp.method = {'sridharchoueiri2020',pVec(ii),4};
    fprintf('Computing for p = %d...\n',pVec(ii));
    [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
    N_sc{ii,1} = dsp.NMatL;
end

clearvars ii

%% Prepare data for plotting

N_sc_samp = zeros(numMus,numPs);
N_sc_mean = zeros(numMus,numPs);
N_sc_max = zeros(numMus,numPs);
N_sc_min = zeros(numMus,numPs);
for ii = 1:numPs
    N_sc_samp(:,ii) = N_sc{ii,1}(:,35); % CHANGE INDEX AS NEEDED
    N_sc_mean(:,ii) = mean(N_sc{ii,1},2);
    N_sc_max(:,ii) = max(N_sc{ii,1},[],2);
    N_sc_min(:,ii) = min(N_sc{ii,1},[],2);
end

clearvars ii

%% Plot N vs p for a few different frequencies and for a chosen source
% distance

% Get existing variables defined in workspace
varsbefore = who;

fig1 = figure();
p1 = plot(pVec,N_sc_samp.');
[fig1,p1,a1,a2] = formatLinePlot(fig1,p1,numMus,...
    'BXLabel','Numerical precision, $p$',...
    'LYLabel','Truncation order, $N$',...
    'BXLim',[-0.7,20.7],...
    'LYLim',[-5,80],...
    'BXTick',0:4:20,...
    'LYTick',0:15:75,...
    'TXTick',0:4:20,...
    'RYTick',0:15:75,...
    'TXTickLabel',{},...
    'RYTickLabel',{},...
    'BXMinorTick','on',...
    'BXMinorTickVal',setdiff(0:20,0:4:20,'stable'),...
    'TXMinorTick','on',...
    'TXMinorTickVal',setdiff(0:20,0:4:20,'stable'),...
    'LYMinorTick','on',...
    'RYMinorTick','on',...
    'includeColor',false,...
    'lineStyleList',{'none'},...
    'markerStyleList',{'o','s','d','^','v'});

% Annotate with text (CHANGE AS NEEDED)
textFlag = true;
if textFlag
    % Positions correspond to actual data values
    t = text(a1,[6.5,10.9,12.7,13.5,13.3],[71,57,42,30,14],...
        {['$\mu = ',num2str(round(muVec(5),1)),'~(f = 18.2$ kHz)'],...
        [num2str(round(muVec(4),1)),'~(12.1)'],...
        [num2str(round(muVec(3),1)),'~(6)'],...
        [num2str(round(muVec(2),1)),'~(3)'],...
        [num2str(round(muVec(1),1)),'~(0.06)']});
    for ii = 1:numMus
        t(ii).Interpreter = 'latex';
        t(ii).FontSize = 22;
        t(ii).FontName = 'Times';
    end
end

% Clear variables generated only by this cell
varsafter = who;
newvars = setdiff(varsafter,varsbefore);
varstokeep = {'fig1','p1','a1','a2'};
varstoremove = setdiff(newvars,varstokeep);
clear(varstoremove{:});
clearvars varstokeep varstoremove varsafter newvars

%% (Optional) Export plot

% Specify plot file name
plotFileName = 'N_vs_p_r-inf_grayscale'; % CHANGE AS NEEDED

exportFigure(fig1,fullfile(plotsPath,'Computation',[plotFileName,'.eps']))
