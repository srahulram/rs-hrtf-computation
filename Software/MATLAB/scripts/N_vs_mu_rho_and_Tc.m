%% Truncation order versus Tc, source frequency and distance.
%
% In this script, we see how truncation order, N, varies with Tc, non-
% dimensional frequency, mu, and non-dimensional source distance, rho, for 
% a given angle of incidence, theta. Dependence of N on theta is not 
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

%% Compute N as a function of mu and rho

% Head parameters
head.a = 0.09;
head.eL = [90,0];

% Source parameters
rhoVec = [inf;3;2;1.5];
rVec = (head.a)*rhoVec;
rVecLen = length(rVec);
source.azVec = 260;
source.elVec = 0;

% Computation parameters
thVec = 1e-6;
numTH = length(thVec);
muVec = logspace(-2,log10(40),400);
dsp.fVec = muVec*getSoundSpeed()/(2*pi*(head.a));

% Compute N
N_Cell = cell(1,rVecLen*numTH);
indx = 1;
for ii = 1:rVecLen
    for jj = 1:numTH
        dsp.method = {'dudamartens1998',thVec(jj),2};
        fprintf('Computing for rho = %3.1f and Tc = %5.4f...\n',...
            round(rhoVec(ii),1),thVec(jj));
        source.r = rVec(ii);
        [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
        N_Cell{1,indx} = dsp.NMatL;
        indx = indx + 1;
    end
end

N_mat = cell2mat(N_Cell);

clearvars ii jj indx

%% Prepare data for plotting

N_gd = round(exp(1)*muVec.');
finalPlotData = [N_gd,N_mat];

%% Plot N as a function of mu, rho, and Tc.

% Get existing variables defined in workspace
varsbefore = who;

txTicks = [10,100,500,1000,2000,5000,10000,20000,40000];
txTickLabels = cell(1,length(txTicks));
for ii = 1:length(txTicks)
    txTickLabels{1,ii} = num2str(txTicks(ii)/1000);
end

clearvars ii

numPlots = size(finalPlotData,2);
fig1 = figure();
p1 = plot(muVec.',finalPlotData);
[fig1,p1,a1,a2] = formatLinePlot(fig1,p1,numPlots,...
    'BXLabel','Non-dimensional frequency, $\mu$',...
    'LYLabel','Truncation order, $N$',...
    'TXLabel','Frequency, $f$ (kHz)',...
    'BXScale','log',...
    'TXScale','log',...
    'BXLim',[0.008,45],...
    'LYLim',[-3,63],...
    'TXLim',[0.008,45]*getSoundSpeed()/(2*pi*(head.a)),...
    'BXTick',[0.01,0.1,0.5,1,2,5,10,40],...
    'LYTick',0:15:60,...
    'TXTick',txTicks,...
    'RYTick',0:15:60,...
    'TXTickLabel',txTickLabels,...
    'RYTickLabel',{},...
    'BXMinorTick','on',...
    'TXMinorTick','on',...
    'LYMinorTick','on',...
    'RYMinorTick','on',...
    'includeColor',false,...
    'lineStyleList',{'--','-','-','-','-'},...
    'markerStyleList',{'none'});

% Annotate with text (CHANGE AS NEEDED)
textFlag = true; % Change as needed
if textFlag
    t = text(a1,[0.01,0.04,0.07,0.04,2],[34,23,14,7,2],...
        {'$\rho = 3/2$','$2$','$3$','$\rho \to \infty$',...
        'Eq. (13)'});
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
plotFileName = 'N_vs_mu_rho_and_Tc_grayscale'; % CHANGE AS NEEDED

exportFigure(fig1,fullfile(plotsPath,'Computation',[plotFileName,'.eps']))

%% (Optional) Export plot data
%
% This is for generating plots using another piece of software (e.g.,
% Mathematica)

% Specify file name for exporting plot data
plotDataFileName = 'N_vs_mu_rho_and_Tc_plotdata'; % CHANGE AS NEEDED

xData = muVec.';
yData = finalPlotData;
a = head.a;
c = getSoundSpeed();
save(fullfile(dataPath,'Plot Data',[plotDataFileName,'.mat']),'xData',...
    'yData','a','c');

clearvars xData yData a c

%% Observations
%
% - For a fixed Tc, N increases with increasing mu and decreasing rho.
%
% - The formula by Gumerov and Duraiswami (dashed line), which is supposed
% to approximate N for rho -> infty, underestimates N for mu < 7 and
% overestimates it for mu > 7.
