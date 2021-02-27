%% Empirical derivation of formula for Nmin
%
% Recommended usage: Run each independent cell as required instead of
% running the entire script all at once.
%
% In this script, we empirically-derive a formula for Nmin as a function of
% mu, rho, and numerical precision, p. This is the formula proposed by
% Sridhar and Choueiri [1] in their paper.
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

%% Load pre-computed data

loadFlag = true; % Set to false to recompute data

if loadFlag
    load(fullfile(dataPath,'Computation','N_formula_derivation_data.mat'))
end

%% Compute RS-HRTF for many mu, rho, and p values

if ~loadFlag
    % Head parameters
    head.a = 0.09;
    head.eL = [90,0];
    
    % Source parameters
    rhoVec = logspace(log10(2),1,10).';
    rVec = (head.a)*rhoVec;
    rVecLen = length(rVec);
    source.azVec = 90:5:180;
    source.elVec = 0;
    
    % Computation parameters
    dsp.fS = 40000;
    dsp.T = 0.01;
    pVec = 2:15;
    numPs = length(pVec);
    
    N_optimal = cell(numPs,rVecLen);
    for ii = 1:numPs
        dsp.method = {'sridharchoueiri2020',pVec(ii),4};
        for jj = 1:rVecLen
            fprintf('Computing for p = %d and rho = %3.1f...\n',...
                pVec(ii),round(rhoVec(jj),1));
            source.r = rVec(jj);
            [~,~,dsp] = batchComputeSphereHRTFs(source,head,dsp);
            N_optimal{ii,jj} = round(mean(dsp.NMatL,2));
        end
    end
    
    clearvars ii jj
end

%% (Optional) Export computed data

exportFlag = false; % Change as needed

if exportFlag && ~loadFlag
    exportFileName = 'N_formula_derivation_data';
    
    [~,exportFileName,~] = fileparts(exportFileName); % Removes extension
    save(fullfile(dataPath,'Computation',[exportFileName,'.mat']),...
        '-regexp',['^(?!(dataPath|plotsPath|projectPath|',...
        'exportFileName|exportFlag|loadFlag)$).'],'-v7.3');
end

%% Compute variables required for curve fitting

fVec = getFreqVec(dsp.fS,dsp.irLen);
nyqIndx = ceil((dsp.irLen+1)/2);
muVec = 2*pi*fVec(1:nyqIndx)*(head.a)/getSoundSpeed();
N_optimal = N_optimal.';

%% Determine model parameters by curve fitting

a_Mat = zeros(rVecLen,numPs);
b_Mat = zeros(rVecLen,numPs);
c_Mat = zeros(rVecLen,numPs);
rSMat = zeros(rVecLen,numPs);
rmseMat = zeros(rVecLen,numPs);
for ii = 1:numPs
    for jj = 1:rVecLen
        fprintf('Determining fit for p = %d and rho = %3.1f...\n',...
            pVec(ii),round(rhoVec(jj),1));
        
        % Specify fit options
        fitOpts = fitoptions('Method','NonlinearLeastSquares',...
            'Algorithm','Trust-Region',...
            'Lower',[0,-Inf,0],...
            'Upper',[max(N_optimal{jj,ii})+2,Inf,2],...
            'StartPoint',[N_optimal{jj,ii}(2),1,0.5]);
        
        % Define fitting curve
        fitDef = fittype('a+(b*mu^c)','coefficients',{'a','b','c'},...
            'dependent',{'N_fit'},'independent',{'mu'},'options',fitOpts);
        
        % Generate fit
        [N_fit,gof,~] = fit(muVec(2:nyqIndx),...
            N_optimal{jj,ii}(2:nyqIndx),fitDef,fitOpts);
        
        % Get coefficients and fit-quality metric values
        coeffVals = coeffvalues(N_fit);
        a_Mat(jj,ii) = coeffVals(1);
        b_Mat(jj,ii) = coeffVals(2);
        c_Mat(jj,ii) = coeffVals(3);
        rSMat(jj,ii) = gof.rsquare;
        rmseMat(jj,ii) = gof.rmse;
    end
end

clearvars ii jj coeffVals N_fit gof fitDef fitOpts

%% Observations
%
% - Based on analyzing the obtained fits, we ignore the fits for the
% smallest two values of rho, starting instead at rho = 1.3877.

%% Fit curve to previous model coefficients as a function of rho
%
% We find that using a first-order rational function to compute a, b, and c 
% from rho appears to work best.  We denote the numerator coefficients by 
% p1 and p2, and the denominator coefficient by q1. We generate fits only 
% for p >= 2.

% Coefficients for a
p1_a = zeros(numPs,1);
p2_a = zeros(numPs,1);
q1_a = zeros(numPs,1);
err_a = cell(numPs,1);

% Coefficients for b
p1_b = zeros(numPs,1);
p2_b = zeros(numPs,1);
q1_b = zeros(numPs,1);
err_b = cell(numPs,1);

% Coefficients for c
p1_c = zeros(numPs,1);
p2_c = zeros(numPs,1);
q1_c = zeros(numPs,1);
err_c = cell(numPs,1);

for ii = 1:numPs
    fprintf('Determining fit for p = %d...\n',pVec(ii));
    
    % *** Compute coefficients for a ***
    
    % Specify fit options
    fitOpts = fitoptions('Method','NonlinearLeastSquares',...
        'Algorithm','Trust-Region',...
        'Lower',[0,-Inf,-2],...
        'Upper',[Inf,Inf,Inf],...
        'StartPoint',[2,10,-1]);
    
    % Define fitting curve
    fitDef = fittype('rat11');
    
    % Generate fit
    [rhoFit,gof,~] = fit(rhoVec,a_Mat(:,ii),fitDef,fitOpts);
    
    % Get coefficients and fit-quality metric values
    coeffVals = coeffvalues(rhoFit);
    p1_a(ii) = coeffVals(1);
    p2_a(ii) = coeffVals(2);
    q1_a(ii) = coeffVals(3);
    err_a{ii,1} = zeros(2,1);
    err_a{ii,1}(1) = gof.rsquare;
    err_a{ii,1}(2) = gof.rmse;
    
    % *** Compute coefficients for b ***
    
    % Specify fit options
    fitOpts = fitoptions('Method','NonlinearLeastSquares',...
        'Algorithm','Trust-Region',...
        'Lower',[0,-Inf,-2],...
        'Upper',[Inf,Inf,Inf],...
        'StartPoint',[2,2.5,-1]);
    
    % Define fitting curve
    [rhoFit,gof,~] = fit(rhoVec,b_Mat(:,ii),fitDef,fitOpts);
    
    % Get coefficients and fit-quality metric values
    coeffVals = coeffvalues(rhoFit);
    p1_b(ii) = coeffVals(1);
    p2_b(ii) = coeffVals(2);
    q1_b(ii) = coeffVals(3);
    err_b{ii,1} = zeros(2,1);
    err_b{ii,1}(1) = gof.rsquare;
    err_b{ii,1}(2) = gof.rmse;
    
    % *** Compute coefficients for c ***
    
    % Specify fit options
    fitOpts = fitoptions('Method','NonlinearLeastSquares',...
        'Algorithm','Trust-Region',...
        'Lower',[0,-Inf,-2],...
        'Upper',[Inf,Inf,Inf],...
        'StartPoint',[0.8,-1,-1]);
    
    % Generate fit
    [rhoFit,gof,~] = fit(rhoVec,c_Mat(:,ii),fitDef,fitOpts);
    
    % Get coefficients and fit-quality metric values
    coeffVals = coeffvalues(rhoFit);
    p1_c(ii) = coeffVals(1);
    p2_c(ii) = coeffVals(2);
    q1_c(ii) = coeffVals(3);
    err_c{ii,1} = zeros(2,1);
    err_c{ii,1}(1) = gof.rsquare;
    err_c{ii,1}(2) = gof.rmse;
end

clearvars ii fitOpts fitDef rhoFit gof coeffVals

%% Fit curves to previous model coefficients as a function of p
%
% Note: The results in this cell are not published in the paper by Sridhar
% and Choueiri [1].
%
% Use the curve fitting toolbox GUI (type cftool in the Command Window) to
% generate these last set of fits.

% *** Results ***
%
% For a_Mat:
%
% p1_a = 0.3442*pVec - 0.9584;
% p2_a = 3.369*pVec - 3.69;
% q1_a = -(0.7605*pVec + 1.704)./(pVec + 0.7461);
%
% For b_Mat:
%
% p1_b = 1.758*exp(0.1033*pVec);
% p2_b = -(2.043*pVec - 0.1245);
% p3_b = 2.161*pVec - 1.964;
% q1_b = 0.001728*pVec.^2.977 - 3.088;
% q2_b = 2.679 - 0.0008429*pVec.^3.343;
%
% For c_Mat:
%
% p1_c = 1.017-0.1041*pVec.^0.5086;
% p2_c = 0.1267*pVec - 2.859;
% p3_c = -(0.1227*pVec - 2.335);
% q1_c = -(0.02058*pVec + 2.968);
% q2_c = 0.06702*pVec + 2.169;

%% Generate new surface fit
%
% Note: The results in this cell are not published in the paper by Sridhar
% and Choueiri [1].
%
% Generate surface fit for a_Mat, b_Mat and c_Mat as a function of both p
% and rho.

% *** Determine fit for a_Mat ***

% Specify fit options
fitOpts = fitoptions('Method','NonlinearLeastSquares',...
    'Algorithm','Trust-Region');

% Define fitting curve
fitDef = fittype('(((a1*p-a2)*rho)+(a3*p-a4))/(rho-((b1*p+b2)/(p+b3)))',...
    'coefficients',{'a1','a2','a3','a4','b1','b2','b3'},...
    'dependent',{'aMat_fit'},'independent',{'p','rho'},'options',fitOpts);

% Generate fit
[XOut,YOut,ZOut] = prepareSurfaceData(pVec,rhoVec,a_Mat);
[aMat_fit,gof,~] = fit([XOut,YOut],ZOut,fitDef,fitOpts);

% Get coefficients and fit quality metric values
coeffNames_a = coeffnames(aMat_fit);
coeffVals_a = coeffvalues(aMat_fit);
rSMat_a = gof.rsquare;
rmseMat_a = gof.rmse;

% *** Determine fit for b_Mat ***

% Specify fit options
fitOpts = fitoptions('Method','NonlinearLeastSquares',...
    'Algorithm','Trust-Region');

% Define fitting curve
fitDef = fittype(['((a1*exp(a2*p)*rho^2)-((a3*p-a4)*rho)+(a5*p-a6))/',...
    '(rho^2+((b1*p^b2-b3)*rho)+(b4-(b5*p^b6)))'],...
    'coefficients',{'a1','a2','a3','a4','a5','a6','b1','b2','b3','b4',...
    'b5','b6'},'dependent',{'bMat_fit'},'independent',{'p','rho'},...
    'options',fitOpts);

% Generate fit
[XOut,YOut,ZOut] = prepareSurfaceData(pVec,rhoVec,b_Mat);
[bMat_fit,gof,~] = fit([XOut,YOut],ZOut,fitDef,fitOpts);

% Get coefficients and fit quality metric values
coeffNames_b = coeffnames(bMat_fit);
coeffVals_b = coeffvalues(bMat_fit);
rSMat_b = gof.rsquare;
rmseMat_b = gof.rmse;

% *** Determine fit for c_Mat ***

% Specify fit options
fitOpts = fitoptions('Method','NonlinearLeastSquares',...
    'Algorithm','Trust-Region');

% Define fitting curve
fitDef = fittype(['(((a1-(a2*p^a3))*rho^2)+((a4*p-a5)*rho)-(a6*p-a7))/',...
    '(rho^2-((b1*p+b2)*rho)+(b3*p+b4))'],'coefficients',{'a1','a2',...
    'a3','a4','a5','a6','a7','b1','b2','b3','b4'},'dependent',...
    {'cMat_fit'},'independent',{'p','rho'},'options',fitOpts);

% Generate fit
[XOut,YOut,ZOut] = prepareSurfaceData(pVec,rhoVec.',c_Mat);
[cMat_fit,gof,~] = fit([XOut,YOut],ZOut,fitDef,fitOpts);

% Get coefficients and fit-quality metric values
coeffNames_c = coeffnames(cMat_fit);
coeffVals_c = coeffvalues(cMat_fit);
rSMat_c = gof.rsquare;
rmseMat_c = gof.rmse;

%% Export model parameters

mpeFlag = false; % Change as needed

if mpeFlag
    save(fullfile(dataPath,'Computation','N_formula_model_params.mat'),...
        'fVec','nyqIndx','muVec','N_optimal',...
        'a_Mat','b_Mat','c_Mat','rSMat','rmseMat',...
        'pVec','numPs','rhoVec','rVec','rVecLen',...
        'p1_a','p2_a','q1_a','err_a',...
        'p1_b','p2_b','q1_b','err_b',...
        'p1_c','p2_c','q1_c','err_c',...
        '-v7.3');
end
