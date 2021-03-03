%% Optimization of modified iterative approach based on ILD distortion
%
% Recommended usage: Run each independent cell as required instead of
% running the entire script all at once.
%
% Determine the smallest value of numerical precision, p, to use in the 
% modified iterative approach described by Sridhar and Choueiri [1] in 
% order to compute the RS-HRTF such that the maximum absolute ILD 
% distortion (relative to benchmarks) for mu < 8 and for sources within 
% +/- 90 deg. of the median plane is less than 1 dB. Assume antipodal ears.
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
    load(fullfile(dataPath,'Computation','pmin_noILDDistortion_data.mat'))
end

%% Compute reference ITD spectra

if ~loadFlag
    % Common parameters for RS-HRTF computation
    % See 'batchComputeSphereHRTFs' function for parameter descriptions
    
    % Head parameters
    head.a = 0.09;
    head.eL = [90,0];
    head.eR = [270,0];
    
    % Source parameters
    rhoVec = [logspace(log10(1.5),log10(20),20).';inf];
    rVec = (head.a)*rhoVec;
    rVecLen = length(rVec);
    source.azVec = mod(-90:3:90,360);
    source.elVec = 0;
    
    % Computation parameters
    dsp.fS = 48000;
    dsp.T = 0.01;
    dsp.causalflag = true;
    dsp.method = {'exact',inf};
    
    % Main calculation
    refILD = cell(rVecLen,1);
    for ii = 1:rVecLen
        fprintf('Computing for rho = %3.1f...\n',round(rhoVec(ii),1));
        source.r = rVec(ii);
        [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
        refILD{ii,1} = estimateILD(head.hrirL,head.hrirR);
    end
    
    clearvars ii
end

%% Iteratively determine lowest value of p for desired accuracy.

if ~loadFlag
    muC = 8; % Approximate upper cut-off value of mu for computing IPD
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
            testILD = estimateILD(head.hrirL,head.hrirR);
            diffILD = refILD{jj,1}(1:fcIndx,:)-testILD(1:fcIndx,:);
            maxDistMat(ii,jj) = max(max(abs(diffILD)));
            meanDistMat(ii,jj) = mean(mean(abs(diffILD)));
            minDistMat(ii,jj) = min(min(abs(diffILD)));
            fprintf('Current value of maxDist (in dB): %d\n',...
                round(maxDistMat(ii,jj)));
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
    exportFileName = 'pmin_noILDDistortion_data'; % CHANGE AS NEEDED
    
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
e1 = errorbar(pVec.',meanDistVec,[],(maxDistVec-meanDistVec));
hold on
p1 = plot(linspace(0,6,50).',ones(50,1));
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
    'LYLabel','Absolute ILD error (dB)',...
    'BXScale','linear',...
    'TXScale','linear',...
    'LYScale','log',...
    'RYScale','linear',...
    'BXLim',[0.8,5.2],...
    'LYLim',[0.33e-6,3e2],...
    'TXLim',[0.8,5.2],...
    'RYLim',log10([0.33e-6,3e2]),...
    'BXTick',1:5,...
    'LYTick',1.*10.^((-6:1:2).'),...
    'RYTick',log10(1.*10.^((-6:1:2).')),...
    'LYTickLabel',{'$10^{-6}$','','$10^{-4}$','','$10^{-2}$','','$1$',...
        '','$10^{2}$'},...
    'RYTickLabel',{},...
    'includeColor',false,...
    'lineStyleList',{'--'});

% Draw Patch
ptch = patch(a2,[0.8,5.2,5.2,0.8],[0,0,log10(3e2),log10(3e2)],...
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
plotFileName = 'abs_ILD_err_vs_p_grayscale'; % CHANGE AS NEEDED

exportFigure(fig1,fullfile(plotsPath,'Computation',[plotFileName,'.eps']))
