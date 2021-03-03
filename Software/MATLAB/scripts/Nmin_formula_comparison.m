%% Comparison of formulas for RS-HRTF truncation order
%
% Recommended usage: Run each independent cell as required instead of
% running the entire script all at once.
%
% In this script, we compare the formula for Nmin derived using the script
% Nmin_formula_derivation.m and the rule-of-thumb formula derived by 
% Gumerov and Duraiswami (see the paper by Sridhar and Choueiri [1] for 
% more information). This script generates the plots shown in Fig. 6 in the 
% paper by Sridhar and Choueiri [1].
%
% Note: Calculations done for left ear HRTF only.
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

loadFlag = true; % Set to false to recompute data

if loadFlag
    load(fullfile(dataPath,'Computation','N_formula_comparison_data.mat'))
end

%% Compute IPD spectra

if ~loadFlag
    % Specify parameters for RS-HRTF computation
    % See 'batchComputeSphereHRTFs' function for parameter descriptions
    
    % Head parameters
    head.a = 0.09;
    head.eL = [90,0];
    head.eR = [270,0];
    
    % Source parameters
    rhoVec = [logspace(log10(2),1,9).';inf];
    rVec = (head.a)*rhoVec;
    rVecLen = length(rVec);
    source.azVec = mod(-30:1:30,360);
    source.elVec = 0;
    
    % Computation parameters
    dsp.fS = 48000;
    dsp.T = 0.01;
    dsp.causalflag = true;
    
    % Main calculation    
%     refIPD = cell(rVecLen,1);
    optIPD = cell(rVecLen,1);
    scIPD = cell(rVecLen,1);
    gdIPD = cell(rVecLen,1);
%     N_ref = cell(1,rVecLen);
    N_opt = cell(1,rVecLen);
    N_sc = cell(1,rVecLen);
    N_gd = cell(1,rVecLen);
%     T_ref = zeros(rVecLen,1);
    T_opt = zeros(rVecLen,1);
    T_sc = zeros(rVecLen,1);
    T_gd = zeros(rVecLen,1);
    for ii = 1:rVecLen
%         % Compute reference IPD
%         fprintf('Computing reference for rho = %3.1f...\n',...
%             round(rhoVec(ii),1));
%         source.r = rVec(ii);
%         dsp.method = {'exact',inf};
%         tStart = tic;
%         [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
%         T_ref(ii,1) = toc(tStart);
%         N_ref{1,ii} = median(dsp.NMatL,2);
%         refIPD{ii,1} = estimateIPD(head.hrirL,head.hrirR,dsp.fS);
        
        % Compute optimal IPD (using iterative approach)
        fprintf('Computing optimal for rho = %3.1f...\n',...
            round(rhoVec(ii),1));
        source.r = rVec(ii);
        dsp.method = {'sridharchoueiri2020',4,4};
        tStart = tic;
        [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
        T_opt(ii,1) = toc(tStart);
        N_opt{1,ii} = round(mean(dsp.NMatL,2));
        optIPD{ii,1} = estimateIPD(head.hrirL,head.hrirR,dsp.fS);
        
        % Compute IPD when using Sridhar and Choueiri's formula for Nmin
        fprintf('Computing SC for rho = %3.1f...\n',round(rhoVec(ii),1));
        source.r = rVec(ii);
        dsp.method = {'formulanmin'};
        tStart = tic;
        [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
        T_sc(ii,1) = toc(tStart);
        N_sc{1,ii} = round(mean(dsp.NMatL,2));
        scIPD{ii,1} = estimateIPD(head.hrirL,head.hrirR,dsp.fS);
        
        % Compute IPD when using Gumerov and Duraiswami's formula for N
        fprintf('Computing GD for rho = %3.1f...\n',round(rhoVec(ii),1));
        source.r = rVec(ii);
        dsp.method = {'formulagd'};
        tStart = tic;
        [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
        T_gd(ii,1) = toc(tStart);
        N_gd{1,ii} = round(mean(dsp.NMatL,2));
        gdIPD{ii,1} = estimateIPD(head.hrirL,head.hrirR,dsp.fS);
    end
    
    clearvars ii tStart
end

%% Compute ILD spectra

if ~loadFlag
    % Specify parameters for RS-HRTF computation
    % See 'batchComputeSphereHRTFs' function for parameter descriptions
    
    % Head parameters
    head.a = 0.09;
    head.eL = [90,0];
    head.eR = [270,0];
    
    % Source parameters
    rhoVec = [logspace(log10(2),1,9).';inf];
    rVec = (head.a)*rhoVec;
    rVecLen = length(rVec);
    source.azVec = mod(-90:3:90,360);
    source.elVec = 0;
    
    % Computation parameters
    dsp.fS = 48000;
    dsp.T = 0.01;
    dsp.causalflag = true;
    
    % Main calculation    
%     refILD = cell(rVecLen,1);
    optILD = cell(rVecLen,1);
    scILD = cell(rVecLen,1);
    gdILD = cell(rVecLen,1);
    for ii = 1:rVecLen
%         % Compute reference ILD
%         fprintf('Computing reference for rho = %3.1f...\n',...
%             round(rhoVec(ii),1));
%         source.r = rVec(ii);
%         dsp.method = {'exact',inf};
%         [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
%         refILD{ii,1} = estimateILD(head.hrirL,head.hrirR);
        
        % Compute optimal ILD (using iterative approach)
        fprintf('Computing optimal for rho = %3.1f...\n',...
            round(rhoVec(ii),1));
        source.r = rVec(ii);
        dsp.method = {'sridharchoueiri2020',4,4};
        [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
        optILD{ii,1} = estimateILD(head.hrirL,head.hrirR);
        
        % Compute ILD when using Sridhar and Choueiri's formula for Nmin
        fprintf('Computing SC for rho = %3.1f...\n',round(rhoVec(ii),1));
        source.r = rVec(ii);
        dsp.method = {'formulanmin'};
        [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
        scILD{ii,1} = estimateILD(head.hrirL,head.hrirR);
        
        % Compute ILD when using Gumerov and Duraiswami's formula for N
        fprintf('Computing GD for rho = %3.1f...\n',round(rhoVec(ii),1));
        source.r = rVec(ii);
        dsp.method = {'formulagd'};
        [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
        gdILD{ii,1} = estimateILD(head.hrirL,head.hrirR);
    end
    
    clearvars ii
end

%% Compute IPD and ILD spectra for simplified version of the formula 
% proposed by Sridhar and Choueiri [1]. 
%
% This simplified formula is the final, published version of the formula 
% and the calculations performed here are used to generate the plots in
% Fig. 6 in the paper by Sridhar and Choueiri [1].

if ~loadFlag
    dsp.method = {'formulanminsimp'};
    scIPD_simp = cell(rVecLen,1);
    scILD_simp = cell(rVecLen,1);
    N_sc_simp = cell(1,rVecLen);
    T_sc_simp = zeros(rVecLen,1);
    for ii = 1:rVecLen
        source.r = rVec(ii);
        
        % Compute IPD when using Sridhar and Choueiri's formula for Nmin
        fprintf('Computing IPD for rho = %3.1f...\n',round(rhoVec(ii),1));
        source.azVec = mod(-30:1:30,360);
        tStart = tic;
        [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
        T_sc_simp(ii,1) = toc(tStart);
        N_sc_simp{1,ii} = round(mean(dsp.NMatL,2));
        scIPD_simp{ii,1} = estimateIPD(head.hrirL,head.hrirR,dsp.fS);
        
        % Compute ILD when using Sridhar and Choueiri's formula for Nmin
        fprintf('Computing ILD for rho = %3.1f...\n',round(rhoVec(ii),1));
        source.azVec = mod(-90:3:90,360);
        [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp);
        scILD_simp{ii,1} = estimateILD(head.hrirL,head.hrirR);
    end
    
    clearvars ii tStart
end

%% (Optional) Export computed data

exportFlag = false;

if exportFlag && ~loadFlag
    exportFileName = 'N_formula_comparison_data'; % CHANGE AS NEEDED
    
    [~,exportFileName,~] = fileparts(exportFileName); % Removes extension
    save(fullfile(dataPath,'Computation',[exportFileName,'.mat']),...
        '-regexp',['^(?!(dataPath|plotsPath|projectPath|',...
        'exportFileName|exportFlag|loadFlag)$).'],'-v7.3');
end

%% Compute error metrics and prepare data for analysis
%
% Note: The values for relative increase in computation speed computed here 
% are not the values published in the paper by Sridhar and Choueiri [1] 
% because the timing computation here is not accurate. Please see the 
% script 'sc2020formula_vs_iter_vs_gd2002_timing.m' for timing 
% calculations.

% Compute relative increase in computation speed
relSpeedGain = zeros(rVecLen,4);
relSpeedGain(:,1) = ((T_opt-T_gd)./T_gd)*100;
relSpeedGain(:,2) = ((T_opt-T_sc)./T_sc)*100;
relSpeedGain(:,3) = ((T_opt-T_sc_simp)./T_sc_simp)*100;
relSpeedGain(:,4) = ((T_gd-T_sc_simp)./T_sc_simp)*100;

% Compute errors in truncation orders
N_err_gd = cell2mat(N_opt)-cell2mat(N_gd);
N_err_sc = cell2mat(N_opt)-cell2mat(N_sc);
N_err_sc_simp = cell2mat(N_opt)-cell2mat(N_sc_simp);

% Compute IPD and ILD errors
% Cut-off value in Hz corresponding to muC (rounded to the nearest 100)
muC = 2; % Approximate upper cut-off value of mu for computing IPD
fc = floor(getSoundSpeed()*muC/((head.a)*2*pi*100))*100;
fVec = getFreqVec(dsp.fS,dsp.irLen);
[~,fcIndx1] = min(abs(fVec-fc));

% Cut-off value in Hz corresponding to muC (rounded to the nearest 100)
muC = 8; % Approximate upper cut-off value of mu for computing ILD
fc = floor(getSoundSpeed()*muC/((head.a)*2*pi*100))*100;
fVec = getFreqVec(dsp.fS,dsp.irLen);
[~,fcIndx2] = min(abs(fVec-fc));

% maxIPDDistVec = zeros(rVecLen,1);
% meanIPDDistVec = zeros(rVecLen,1);
maxIPDDistMat = zeros(rVecLen,3);
meanIPDDistMat = zeros(rVecLen,3);
% maxILDDistVec = zeros(rVecLen,1);
% meanILDDistVec = zeros(rVecLen,1);
maxILDDistMat = zeros(rVecLen,3);
meanILDDistMat = zeros(rVecLen,3);
for ii = 1:rVecLen
    % Compute IPD distortions
%     diffIPD = refIPD{ii,1}(1:fcIndx1,:)-optIPD{ii,1}(1:fcIndx1,:);
%     maxIPDDistVec(ii,1) = max(max(abs(diffIPD)))*10^6;
%     meanIPDDistVec(ii,1) = mean(mean(abs(diffIPD)))*10^6;
    
    diffIPD = optIPD{ii,1}(1:fcIndx1,:)-gdIPD{ii,1}(1:fcIndx1,:);
    maxIPDDistMat(ii,1) = max(max(abs(diffIPD)))*10^6;
    meanIPDDistMat(ii,1) = mean(mean(abs(diffIPD)))*10^6;
    
    diffIPD = optIPD{ii,1}(1:fcIndx1,:)-scIPD{ii,1}(1:fcIndx1,:);
    maxIPDDistMat(ii,2) = max(max(abs(diffIPD)))*10^6;
    meanIPDDistMat(ii,2) = mean(mean(abs(diffIPD)))*10^6;
    
    diffIPD = optIPD{ii,1}(1:fcIndx1,:)-scIPD_simp{ii,1}(1:fcIndx1,:);
    maxIPDDistMat(ii,3) = max(max(abs(diffIPD)))*10^6;
    meanIPDDistMat(ii,3) = mean(mean(abs(diffIPD)))*10^6;
    
    % Compute ILD distortions
%     diffILD = refILD{ii,1}(1:fcIndx2,:)-optILD{ii,1}(1:fcIndx2,:);
%     maxILDDistVec(ii,1) = max(max(abs(diffILD)));
%     meanILDDistVec(ii,1) = mean(mean(abs(diffILD)));
    
    diffILD = optILD{ii,1}(1:fcIndx2,:)-gdILD{ii,1}(1:fcIndx2,:);
    maxILDDistMat(ii,1) = max(max(abs(diffILD)));
    meanILDDistMat(ii,1) = mean(mean(abs(diffILD)));
    
    diffILD = optILD{ii,1}(1:fcIndx2,:)-scILD{ii,1}(1:fcIndx2,:);
    maxILDDistMat(ii,2) = max(max(abs(diffILD)));
    meanILDDistMat(ii,2) = mean(mean(abs(diffILD)));
    
    diffILD = optILD{ii,1}(1:fcIndx2,:)-scILD_simp{ii,1}(1:fcIndx2,:);
    maxILDDistMat(ii,3) = max(max(abs(diffILD)));
    meanILDDistMat(ii,3) = mean(mean(abs(diffILD)));
end

clearvars ii diffIPD diffILD muC fc

%% Plot absolute ITD error as a function of rho

% Get existing variables defined in workspace
varsbefore = who;

bxTicks = [2,3,5,10,20];
numBXTicks = length(bxTicks);
bxTickLabels = cell(1,numBXTicks);
for ii = 1:numBXTicks-1
    bxTickLabels{1,ii} = num2str(round(bxTicks(ii),1));
end
bxTickLabels{1,numBXTicks} = '$\to \infty$';

fig1 = figure();
plotRhoVec = [rhoVec(1:rVecLen-1);20];
e1 = errorbar(repmat(plotRhoVec,1,2),meanIPDDistMat(:,[1,3]),[],...
    maxIPDDistMat(:,[1,3]));
hold on
p1 = plot(linspace(0,25,100).',10*ones(100,1));
hold off
[fig1,e1,a1,a2] = formatLinePlot(fig1,e1,2,...
    'TXTick',[],...
    'RYTick',[],...
    'includeColor',false,...
    'lineStyleList',{'none'},...
    'markerStyleList',{'o','s'});
set(fig1,'CurrentAxes',a1);
[fig1,p1,a1,a2] = formatLinePlot(fig1,p1,1,...
    'BXLabel','Non-dimensional source distance, $\rho$',...
    'LYLabel','Absolute ITD error ($\mu$s)',...
    'BXScale','log',...
    'TXScale','log',...
    'LYScale','log',...
    'RYScale','log',...
    'BXLim',[1.8,25],...
    'LYLim',[1e-2,1e3],...
    'TXLim',[1.8,25],...
    'RYLim',[1e-2,1e3],...
    'BXTick',bxTicks,...
    'BXTickLabel',bxTickLabels,...
    'LYTick',1.*10.^((-2:1:3).'),...
    'RYTick',1.*10.^((-2:1:3).'),...
    'LYTickLabel',{'$10^{-2}$','$10^{-1}$','$1$','$10$','$10^{2}$',...
        '$10^{3}$'},...
    'RYTickLabel',{},...
    'includeColor',false,...
    'lineStyleList',{'--'});

% Annotate with text (CHANGE AS NEEDED)
textFlag = true;
if textFlag
    % Positions correspond to actual data values
    t = text(a1,[10.3,10.3],[2*10^2,1.7],{'Eq. (13)','Eq. (19)'});
    t(1).Interpreter = 'latex';
    t(1).FontSize = 22;
    t(1).FontName = 'Times';
    t(2).Interpreter = 'latex';
    t(2).FontSize = 22;
    t(2).FontName = 'Times';
end

% Clear variables generated only by this cell
varsafter = who;
newvars = setdiff(varsafter,varsbefore);
varstokeep = {'fig1','e1','p1','a1','a2'};
varstoremove = setdiff(newvars,varstokeep);
clear(varstoremove{:});
clearvars varstokeep varstoremove varsafter newvars

%% (Optional) Export plot

% Specify plot file name
plotFileName = 'abs_ITD_err_vs_rho_grayscale'; % CHANGE AS NEEDED

exportFigure(fig1,fullfile(plotsPath,'Computation',[plotFileName,'.eps']))

%% Plot absolute ILD error as a function of rho

% Get existing variables defined in workspace
varsbefore = who;

bxTicks = [2,3,5,10,20];
numBXTicks = length(bxTicks);
bxTickLabels = cell(1,numBXTicks);
for ii = 1:numBXTicks-1
    bxTickLabels{1,ii} = num2str(round(bxTicks(ii),1));
end
bxTickLabels{1,numBXTicks} = '$\to \infty$';

fig2 = figure();
plotRhoVec = [rhoVec(1:rVecLen-1);20];
e1 = errorbar(repmat(plotRhoVec,1,2),meanILDDistMat(:,[1,3]),[],...
    maxILDDistMat(:,[1,3]));
hold on
p1 = plot(linspace(0,25,100).',ones(100,1));
hold off
[fig2,e1,a1,a2] = formatLinePlot(fig2,e1,2,...
    'TXTick',[],...
    'RYTick',[],...
    'includeColor',false,...
    'lineStyleList',{'none'},...
    'markerStyleList',{'o','s'});
set(fig2,'CurrentAxes',a1);
[fig2,p1,a1,a2] = formatLinePlot(fig2,p1,1,...
    'BXLabel','Non-dimensional source distance, $\rho$',...
    'LYLabel','Absolute ILD error (dB)',...
    'BXScale','log',...
    'TXScale','log',...
    'LYScale','log',...
    'RYScale','log',...
    'BXLim',[1.8,25],...
    'LYLim',[1e-4,1e2],...
    'TXLim',[1.8,25],...
    'RYLim',[1e-4,1e2],...
    'BXTick',bxTicks,...
    'BXTickLabel',bxTickLabels,...
    'LYTick',1.*10.^((-4:1:2).'),...
    'RYTick',1.*10.^((-4:1:2).'),...
    'LYTickLabel',{'$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$',...
        '$1$','$10^{1}$','$10^{2}$'},...
    'RYTickLabel',{},...
    'includeColor',false,...
    'lineStyleList',{'--'});

% Annotate with text (CHANGE AS NEEDED)
textFlag = true;
if textFlag
    % Positions correspond to actual data values
    t = text(a1,[10.4,10.4],[5,0.2*10^-3],{'Eq. (13)','Eq. (19)'});
    t(1).Interpreter = 'latex';
    t(1).FontSize = 22;
    t(1).FontName = 'Times';
    t(2).Interpreter = 'latex';
    t(2).FontSize = 22;
    t(2).FontName = 'Times';
end

% Clear variables generated only by this cell
varsafter = who;
newvars = setdiff(varsafter,varsbefore);
varstokeep = {'fig2','e1','p1','a1','a2'};
varstoremove = setdiff(newvars,varstokeep);
clear(varstoremove{:});
clearvars varstokeep varstoremove varsafter newvars

%% (Optional) Export plot

% Specify plot file name
plotFileName = 'abs_ILD_err_vs_rho_grayscale'; % CHANGE AS NEEDED

exportFigure(fig2,fullfile(plotsPath,'Computation',[plotFileName,'.eps']))
