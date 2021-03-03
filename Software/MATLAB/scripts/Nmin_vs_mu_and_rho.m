%% Optimal truncation order versus source frequency and distance
%
% In this script, we see how optimal truncation order, Nmin, varies with 
% non-dimensional frequency, mu, and non-dimensional source distance, rho, 
% for a given angle of incidence, theta. Dependence of Nmin on theta is not 
% considered since we find it to be negligible.
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

%% Compute RS-HRTF using the iterative approach 
%
% Perform computation for the frontal position, a fixed value of 
% precision, p, and various values of rho.

% Head parameters
head.a = 0.09;
head.eL = [90,0];
% Source parameters
rhoVec = [2;3;10]; % Use for plotting
rVec = (head.a)*rhoVec;
rVecLen = length(rVec);
source.azVec = 150;
source.elVec = 0;
% Computation parameters
dsp.fS = 40000;
dsp.T = 0.01;

% Compute N_min for various combinations of parameter values
Nmin_opt = cell(1,rVecLen);
Nmin_form = cell(1,rVecLen);
for jj = 1:rVecLen
    fprintf('Computing Nmin iteratively for rho = %2.1f...\n',...
        round(rhoVec(jj),1));
    dsp.method = {'sridharchoueiri2020',4,4};
    source.r = rVec(jj);
    [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
    Nmin_opt{1,jj} = dsp.NMatL;
    
    fprintf('Computing formula-based Nmin for rho = %2.1f...\n',...
        round(rhoVec(jj),1));
    dsp.method = {'formulanminsimp'};
    source.r = rVec(jj);
    [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
    Nmin_form{1,jj} = dsp.NMatL;
end

clearvars jj

%% Generate plot data

Nmin_opt = cell2mat(Nmin_opt);
Nmin_form = cell2mat(Nmin_form);

aoi = round(mod(source.azVec-head.eL(1),360)); % Angle of incidence.
numPlotPts = 20; % Approx. number of data points to plot. Change as needed.
fVec = getFreqVec(dsp.fS,dsp.irLen);
logFVec = getLogFreqVec(dsp.fS,dsp.irLen);
muVec = 2*pi*fVec*head.a/getSoundSpeed;
nyqIndx = ceil((dsp.irLen+1)/2);

% Extract indices of data points to plot
plotFVec = downsample(getLogFreqVec(dsp.fS/2,nyqIndx),...
    floor(nyqIndx/numPlotPts));
numPlotPts = length(plotFVec);
plotPtIndxs = zeros(numPlotPts,1);
for ii = 1:numPlotPts
    [~,plotPtIndxs(ii)] = findNearest(fVec,plotFVec(ii));
end
plotPtIndxs = unique(plotPtIndxs);
numPlotPts = length(plotPtIndxs);

clearvars ii

%% Plot dependency of Nmin for various mu and rho and for fixed p and Theta

% Get existing variables defined in workspace
varsbefore = who;

txTicks = [10,20,50,100,200,500,1000,2000,5000,10000,20000,50000];
txTickLabels = cell(1,length(txTicks));
for ii = 1:length(txTicks)
    txTickLabels{1,ii} = num2str(txTicks(ii)/1000);
end

numPlots1 = size(Nmin_opt,2);
numPlots2 = size(Nmin_form,2);
fig1 = figure();
p1 = plot(muVec(plotPtIndxs),Nmin_opt(plotPtIndxs,:));
hold on
p2 = stairs(muVec(2:nyqIndx),Nmin_form(2:end,:));
hold off
[fig1,p1,a1,a2] = formatLinePlot(fig1,p1,numPlots1,...
    'TXTick',[],...
    'RYTick',[],...
    'includeColor',false,...
    'lineStyleList',{'none'},...
    'markerStyleList',{'o','s','d'},...
    'markerFillFlag',true);
set(fig1,'CurrentAxes',a1);
[fig1,p2,a1,a2] = formatLinePlot(fig1,p2,numPlots2,...
    'BXLabel','Non-dimensional frequency, $\mu$',...
    'LYLabel','Truncation order, $N_{\mathrm{min}}$',...
    'TXLabel','Frequency, $f$ (kHz)',...
    'BXScale','log',...
    'TXScale','log',...
    'LYScale','log',...
    'RYScale','log',...
    'BXLim',[logFVec(1)-15,(dsp.fS/2)+3000]*(2*pi*head.a)/getSoundSpeed,...
    'LYLim',[1.65,70],...
    'TXLim',[logFVec(1)-15,(dsp.fS/2)+3000],...
    'RYLim',[1.65,70],...
    'BXTick',[0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50],...
    'LYTick',[2,3,5,10,20,40,60],...
    'TXTick',txTicks,...
    'RYTick',[2,3,5,10,20,40,60],...
    'TXTickLabel',txTickLabels,...
    'RYTickLabel',{},...
    'BXMinorTick','on',...
    'TXMinorTick','on',...
    'LYMinorTick','on',...
    'RYMinorTick','on',...
    'XMinorGrid','on',...
    'includeColor',false,...
    'lineStyleList',{'-'},...
    'markerStyleList',{'none'});

% Annotate with text
textFlag = true; % Change as needed
if textFlag
    t = text(a1,[0.18,0.23,0.20],[13.5,6,3.5],...
        {'$\rho = 2$','$3$','$10$'});
    for ii = 1:numPlots2
        t(ii).Interpreter = 'latex';
        t(ii).FontSize = 22;
        t(ii).FontName = 'Times';
    end
end

% Clear variables generated only by this cell
varsafter = who;
newvars = setdiff(varsafter,varsbefore);
varstokeep = {'fig1','p1','p2','a1','a2'};
varstoremove = setdiff(newvars,varstokeep);
clear(varstoremove{:});
clearvars varstokeep varstoremove varsafter newvars

%% (Optional) Export plot

% Specify plot file name
plotFileName = 'N-min_vs_mu-and-rho_grayscale'; % CHANGE AS NEEDED

exportFigure(fig1,fullfile(plotsPath,'Computation',[plotFileName,'.eps']))
