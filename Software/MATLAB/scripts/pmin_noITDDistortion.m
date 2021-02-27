%% Optimization of modified iterative approach based on IPD distortion
%
% Recommended usage: Run each independent cell as required instead of
% running the entire script all at once.
%
% Determine the smallest value of numerical precision, p, to use in the 
% modified iterative approach described by Sridhar and Choueiri [1] in 
% order to compute the RS-HRTF such that the maximum absolute IPD 
% distortion (relative to benchmarks) for mu < 2 and for sources within 
% +/- 30 deg. of the median plane is less than 10 microseconds. Assume 
% antipodal ears.
%
% Ref:
%   [1]. R. Sridhar and E. Y. Choueiri, "Optimal Series Truncation of the 
%   Rigid-Sphere Head-Related Transfer Function for Accurate Binaural 
%   Perception," J. Audio Eng. Soc., 2021 (to appear)

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

%% (Optional) Load pre-computed data
%
% Useful for regenerating plots in the paper by Sridhar and Choueiri [1].

loadFlag = true;

if loadFlag
    load(fullfile(dataPath,'Computation','pmin_noIPDDistortion_data.mat'))
end

%% Compute reference IPD spectra

if ~loadFlag
    % Specify parameters for RS-HRTF computation
    % See 'batchComputeSphereHRTFs' function for parameter descriptions
    
    % Head parameters
    head.a = 0.09;
    head.eL = [90,0];
    head.eR = [270,0];
    
    % Source parameters
    rhoVec = [logspace(log10(1.5),log10(20),20).';inf];
    rVec = (head.a)*rhoVec;
    rVecLen = length(rVec);
    source.azVec = mod(-30:1:30,360);
    source.elVec = 0;
    
    % Computation parameters
    dsp.fS = 48000;
    dsp.T = 0.01;
    dsp.causalflag = true;
    dsp.method = {'exact',inf};
    
    % Main calculation
    refIPD = cell(rVecLen,1);
    for ii = 1:rVecLen
        fprintf('Computing for rho = %3.1f...\n',round(rhoVec(ii),1));
        source.r = rVec(ii);
        [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
        refIPD{ii,1} = estimateIPD(head.hrirL,head.hrirR,dsp.fS);
    end
    
    clearvars ii
end

%% Iteratively determine lowest value of p for desired accuracy.

if ~loadFlag
    muC = 2; % Approximate upper cut-off value of mu for computing IPD
    pVec = 1:5; % Determined apriori to save computation time
    numPs = length(pVec);
    
    % Cut-off value in Hz corresponding to muC (rounded to the nearest 100)
    fc = floor(getSoundSpeed()*muC/((head.a)*2*pi*100))*100;
    fVec = getFreqVec(dsp.fS,dsp.irLen);
    [~,fcIndx] = min(abs(fVec-fc));
    maxDistMat = zeros(numPs,rVecLen); 
    meanDistMat = zeros(numPs,rVecLen); 
    minDistMat = zeros(numPs,rVecLen); 
    maxDistVec = zeros(numPs,1);
    meanDistVec = zeros(numPs,1);
    minDistVec = zeros(numPs,1);
    for ii = 1:numPs
        fprintf('Computing for iteration %d (p = %d)\n',ii,pVec(ii));
        dsp.method = {'sridharchoueiri2020',pVec(ii),4};
        for jj = 1:rVecLen
            fprintf('Computing for rho = %3.1f\n',round(rhoVec(jj),1));
            source.r = rVec(jj);
            
            [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
            testIPD = estimateIPD(head.hrirL,head.hrirR,dsp.fS);
            diffIPD = refIPD{jj,1}(1:fcIndx,:)-testIPD(1:fcIndx,:);
            maxDistMat(ii,jj) = max(max(abs(diffIPD)));
            meanDistMat(ii,jj) = mean(mean(abs(diffIPD)));
            minDistMat(ii,jj) = min(min(abs(diffIPD)));
            fprintf('Current value of maxDist (in mu-s): %d\n',...
                round(maxDistMat(ii,jj)*10^6));
        end
        maxDistVec(ii,1) = max(maxDistMat(ii,:));
        meanDistVec(ii,1) = mean(meanDistMat(ii,:));
        minDistVec(ii,1) = min(minDistMat(ii,:));
    end
    
    clearvars ii jj fcIndx testIPD diffIPD
end

%% (Optional) Export computed data

exportFlag = false;

if ~loadFlag && exportFlag
    exportFileName = 'pmin_noIPDDistortion_data'; % CHANGE AS NEEDED
    
    [~,exportFileName,~] = fileparts(exportFileName); % Removes extension
    head = rmfield(head,{'hrirL','hrirR'});
    save(fullfile(dataPath,'Computation',[exportFileName,'.mat']),...
        '-regexp',['^(?!(dataPath|plotsPath|projectPath|',...
        'exportFileName|exportFlag|loadFlag)$).'],'-v7.3');
end

%% Plot max. distortion as a function of p

% Get existing variables defined in workspace
varsbefore = who;

fig1 = figure();
e1 = errorbar(pVec.',meanDistVec*10^6,[],(maxDistVec-meanDistVec)*10^6);
hold on
p1 = plot(linspace(0,6,50).',10*ones(50,1));
hold off
[fig1,e1,a1,a2] = formatLinePlot(fig1,e1,1,...
    'TXTick',[],...
    'RYTick',[],...
    'includeColor',false,...
    'lineStyleList',{'none'},...
    'markerStyleList',{'o'});
set(fig1,'CurrentAxes',a1);
[fig1,p1,a1,a2] = formatLinePlot(fig1,p1,1,...
    'BXLabel','Numerical precision, $p$',...
    'LYLabel','Absolute ITD error ($\mu$s)',...
    'BXScale','linear',...
    'TXScale','linear',...
    'LYScale','log',...
    'RYScale','linear',...
    'BXLim',[0.8,5.2],...
    'LYLim',[0.33e-3,3e3],...
    'TXLim',[0.8,5.2],...
    'RYLim',log10([0.33e-3,3e3]),...
    'BXTick',1:5,...
    'LYTick',1.*10.^((-3:1:3).'),...
    'RYTick',log10(1.*10.^((-3:1:3).')),...
    'LYTickLabel',{'$10^{-3}$','','$10^{-1}$','','$10$','','$10^{3}$'},...
    'RYTickLabel',{},...
    'includeColor',false,...
    'lineStyleList',{'--'});

% Draw Patch
ptch = patch(a2,[0.8,5.2,5.2,0.8],[1,1,log10(3e3),log10(3e3)],...
    [176,176,176]/255);
ptch.FaceColor = 'none';
ptch.EdgeColor = 'none';
ptch.LineStyle = 'none';
hatchfill2(ptch,'single','HatchAngle',45,'HatchDensity',20,...
    'HatchColor',[176,176,176]/255,'HatchLineWidth',1.5);

% Clear variables generated only by this cell
varsafter = who;
newvars = setdiff(varsafter,varsbefore);
varstokeep = {'fig1','e1','p1','a1','a2'};
varstoremove = setdiff(newvars,varstokeep);
clear(varstoremove{:});
clearvars varstokeep varstoremove varsafter newvars

%% (Optional) Export plot

% Specify plot file name
plotFileName = 'abs_ITD_err_vs_p_grayscale'; % CHANGE AS NEEDED

exportFigure(fig1,fullfile(plotsPath,'Computation',[plotFileName,'.pdf']))
