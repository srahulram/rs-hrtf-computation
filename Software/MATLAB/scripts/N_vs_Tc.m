%% Truncation order versus convergence threshold.
%
% In this script, we see how truncation order, N, varies with convergence
% threshold, Tc.
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
thVec = logspace(-7,-3,200);
numTH = length(thVec);

%% Compute N for various T_c using the method by Cooper and Bauck (1980)

fprintf('Computing using ''cooperbauck1980'' method...\n');

% Compute N for various Tc values
N_cb = cell(numTH,1);
for ii = 1:numTH
    dsp.method = {'cooperbauck1980',thVec(ii),2};
    fprintf('Computing for T_c = 10^(%2.1f)...\n',...
        round(log10(thVec(ii)),1));
    [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
    N_cb{ii,1} = dsp.NMatL;
end

clearvars ii

%% Prepare data for plotting

N_cb_samp = zeros(numMus,numTH);
N_cb_mean = zeros(numMus,numTH);
N_cb_max = zeros(numMus,numTH);
N_cb_min = zeros(numMus,numTH);
for ii = 1:numTH
    N_cb_samp(:,ii) = N_cb{ii,1}(:,35); % CHANGE INDEX AS NEEDED
    N_cb_mean(:,ii) = mean(N_cb{ii,1},2);
    N_cb_max(:,ii) = max(N_cb{ii,1},[],2);
    N_cb_min(:,ii) = min(N_cb{ii,1},[],2);
end

clearvars ii

%% Plot N as a function of mu, rho, and Tc.

% Get existing variables defined in workspace
varsbefore = who;

finalPlotData = N_cb_samp.';
numPlots = size(finalPlotData,2);
fig1 = figure();
p1 = plot(thVec,finalPlotData);
[fig1,p1,a1,a2] = formatLinePlot(fig1,p1,numPlots,...
    'BXLabel','Convergence threshold, $T_c$',...
    'LYLabel','Truncation order, $N$',...
    'BXScale','log',...
    'TXScale','log',...
    'BXLim',[0.85*10^(-7),1.2*10^(-3)],...
    'LYLim',[-3,63],...
    'TXLim',[0.85*10^(-7),1.2*10^(-3)],...
    'RYLim',[-3,63],...
    'BXTick',10.^(-7:-3),...
    'LYTick',0:15:60,...
    'TXTick',10.^(-7:-3),...
    'RYTick',0:15:60,...
    'BXTickLabel',{'$10^{-7}$','$10^{-6}$','$10^{-5}$','$10^{-4}$',...
        '$10^{-3}$'},...
    'TXTickLabel',{},...
    'RYTickLabel',{},...
    'BXMinorTick','on',...
    'TXMinorTick','on',...
    'LYMinorTick','on',...
    'RYMinorTick','on',...
    'includeColor',false,...
    'lineStyleList',{'-'},...
    'markerStyleList',{'none'});

% Annotate with text (CHANGE AS NEEDED)
textFlag = true;
if textFlag
    % Positions correspond to actual data values
    t = text(a1,[4e-7,6.5e-6,1.5e-5,4e-5,8e-5],[7,14.5,23.5,34.5,...
        47.5],{['$\mu = ',num2str(round(muVec(1),1)),...
        '~(f = 0.06\,\mathrm{kHz}$)'],...
        [num2str(round(muVec(2),1)),' (3)'],...
        [num2str(round(muVec(3),1)),' (6.1)']...
        [num2str(round(muVec(4),1)),' (12.1)'],...
        [num2str(round(muVec(5),1)),' (18.2)']});
    for ii = 1:numPlots
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
plotFileName = 'N_vs_Tc_r-inf_grayscale'; % CHANGE AS NEEDED

exportFigure(fig1,fullfile(plotsPath,'Computation',[plotFileName,'.eps']))

%% Observations
%
% - For any given value of rho and theta, we see that N increases with
% decreasing Tc, as expected, provided Tc is decreased by a sufficiently
% large amount that is a complicated function of mu, rho, and theta. 
%
% - For any given value of rho and mu, N doesn't vary significantly with
% theta.
%
% - For any given value of mu and theta, N increases with decreasing rho,
% with the increase in N being significant for near-field sources (i.e., 
% rho less than approx. 10).
