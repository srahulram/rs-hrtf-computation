%% Compare computation times using stripped-down version of RS-HRTF 
% computation code.
%
% Recommended usage: Run each independent cell as required instead of
% running the entire script all at once.
%
% In this script, we use a stripped-down version of the RS-HRTF computation
% code found in the 3D3A Research Lite Toolbox and 3D3A MATLAB Toolbox to
% measure and compare computation times for various computation approaches.
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
% Useful for verifying results published by Sridhar and Choueiri [1].

loadFlag = true;

if loadFlag
    load(fullfile(dataPath,'Computation',...
        'sc2020formula_vs_iter_vs_gd2002.mat'))
end

%% Define required variables

if ~loadFlag
    r_sph = 0.09; % head radius in m
    rhoVec = [2;3;5;inf];
    numRhos = length(rhoVec);
    thetaVec = 0:45:180;
    numThetas = length(thetaVec);
    irDur = 2; % in ms
    Fs = 40000;
    irLen = round(irDur*Fs/1000);
    fVec = getFreqVec(Fs,irLen);
    nyqIndx = ceil((irLen+1)/2);
    muVec = 2*pi*fVec(2:nyqIndx)*r_sph/getSoundSpeed();
    numMus = length(muVec);
    pVec = 4; % Must be finite
    numPs = length(pVec);
    numChkTerms = 4;
end

%% Compute additional variables

if ~loadFlag
    rhoMat = repmat(rhoVec.',numMus,1);
    muMat = repmat(muVec,1,numRhos);
    muRhoMat = muMat.*rhoMat;
    cosThetaVec = cosd(thetaVec);
end

%% Perform and time iterative approach to RS-HRTF computation
%
% This is based on the modified iterative approach described by Sridhar and 
% Choueiri [1].

if ~loadFlag
    % Get existing variables defined in workspace
    varsbefore = who;
    
    H_iter = cell(numPs,numThetas);
    CT_iter = cell(numPs,numThetas);
    N_iter = cell(numPs,numThetas);
    for pp = 1:numPs
        p = pVec(pp);
        for ii = 1:numThetas
            H_iter{pp,ii} = zeros(numMus,numRhos);
            CT_iter{pp,ii} = zeros(numMus,numRhos);
            N_iter{pp,ii} = zeros(numMus,numRhos);
            cosTheta = cosThetaVec(ii);
            for jj = 1:numRhos
                for kk = 1:numMus
                    rho = rhoMat(kk,jj);
                    mu = muMat(kk,jj);
                    mr = muRhoMat(kk,jj);
                    fprintf(['Computing for p = %d, theta = %d, rho =',...
                        ' %d, and mu = %3.1f...\n'],p,thetaVec(ii),rho,...
                        round(mu,1));
                    
                    if rho == inf
                        % Perform 10 untimed repetitions
                        for qq = 1:10
                            % Initialize hf to store 3 most recent terms
                            hf = zeros(3,1);
                            % Initialize P to store 3 most recent terms
                            P = zeros(3,1);
                            hf(2) = exp(1i*mu)/(1i*mu); % m = 0
                            hf(3) = hf(2)*((1/mu)-1i); % m = 1
                            
                            % Initialize sums
                            psiPS = 0;
                            psiPSVec = zeros(numChkTerms+1,1);
                            
                            m = 0;
                            P(3) = 1; % P value for m = 0
                            dh = computeDH(hf(2:3),m,mu); % m = 0
                            Am = computeAm(m,dh); % m = 0
                            psiPS = psiPS + conj(Am)*P(3);
                            termIndx = 1;
                            psiPSVec(termIndx) = psiPS;
                            lF = true;
                            % Move current value out of P(3)
                            P = circshift(P,-1);
                            % Update current value (i.e., for m = 1)
                            P(3) = cosTheta;
                            while lF
                                m = m + 1;
                                hf = circshift(hf,-1);
                                hf(3) = computeHX(hf(1:2),m+1,mu);
                                dh = computeDH(hf(2:3),m,mu);
                                Am = computeAm(m,dh);
                                psiPS = psiPS + conj(Am)*P(3);
                                termIndx = termIndx + 1;
                                psiPSVec(termIndx) = psiPS;
                                if termIndx == (numChkTerms+1)
                                    metricVec = zeros(numChkTerms,1);
                                    for zz = 1:numChkTerms
                                        metricVec(zz) = psiPSVec(zz+1)-...
                                            psiPSVec(zz);
                                    end
                                    realVec = round(real(metricVec),p);
                                    imagVec = round(imag(metricVec),p);
                                    lF1 = any([realVec;imagVec] ~= 0);
                                    lF2 = ~any(isnan([realVec;imagVec]));
                                    lF = lF1 && lF2;
                                    termIndx = termIndx - 1;
                                    psiPSVec = circshift(psiPSVec,-1);
                                else
                                    lF = true;
                                end
                                P = circshift(P,-1);
                                P(3) = computeP(P(1:2),m+1,cosTheta);
                            end
                            psiPSVec = circshift(psiPSVec,1);
                        end
                        
                        % Perform 1000 timed repetitions
                        tStart = tic;
                        for qq = 1:1000
                            % Initialize hf to store 3 most recent terms
                            hf = zeros(3,1);
                            % Initialize P to store 3 most recent terms
                            P = zeros(3,1);
                            hf(2) = exp(1i*mu)/(1i*mu); % m = 0
                            hf(3) = hf(2)*((1/mu)-1i); % m = 1
                            
                            % Initialize sums
                            psiPS = 0;
                            psiPSVec = zeros(numChkTerms+1,1);
                            
                            m = 0;
                            P(3) = 1; % P value for m = 0
                            dh = computeDH(hf(2:3),m,mu); % m = 0
                            Am = computeAm(m,dh); % m = 0
                            psiPS = psiPS + conj(Am)*P(3);
                            termIndx = 1;
                            psiPSVec(termIndx) = psiPS;
                            lF = true;
                            % Move current value out of P(3)
                            P = circshift(P,-1);
                            % Update current value (i.e., for m = 1)
                            P(3) = cosTheta;
                            while lF
                                m = m + 1;
                                hf = circshift(hf,-1);
                                hf(3) = computeHX(hf(1:2),m+1,mu);
                                dh = computeDH(hf(2:3),m,mu);
                                Am = computeAm(m,dh);
                                psiPS = psiPS + conj(Am)*P(3);
                                termIndx = termIndx + 1;
                                psiPSVec(termIndx) = psiPS;
                                if termIndx == (numChkTerms+1)
                                    metricVec = zeros(numChkTerms,1);
                                    for zz = 1:numChkTerms
                                        metricVec(zz) = psiPSVec(zz+1)-...
                                            psiPSVec(zz);
                                    end
                                    realVec = round(real(metricVec),p);
                                    imagVec = round(imag(metricVec),p);
                                    lF1 = any([realVec;imagVec] ~= 0);
                                    lF2 = ~any(isnan([realVec;imagVec]));
                                    lF = lF1 && lF2;
                                    termIndx = termIndx - 1;
                                    psiPSVec = circshift(psiPSVec,-1);
                                else
                                    lF = true;
                                end
                                P = circshift(P,-1);
                                P(3) = computeP(P(1:2),m+1,cosTheta);
                            end
                            psiPSVec = circshift(psiPSVec,1);
                        end
                        tStop = toc(tStart);
                        H_iter{pp,ii}(kk,jj) = (1/mu^2)*psiPSVec(1);
                        
                        CT_iter{pp,ii}(kk,jj) = tStop/1000;
                    else
                        for qq = 1:10
                            % Initialize hf to store 3 most recent terms
                            hf = zeros(3,1);
                            % Initialize P to store 3 most recent terms
                            P = zeros(3,1);
                            hf(2) = exp(1i*mu)/(1i*mu); % m = 0
                            hf(3) = hf(2)*((1/mu)-1i); % m = 1
                            
                            % Initialize sums
                            psiPS = 0;
                            psiPSVec = zeros(numChkTerms+1,1);
                            
                            % Init. hn to store 3 most recent terms
                            hn = zeros(3,1);
                            % Note: hn is the variable storing the value of
                            % the spherical-Hankel function (of the first 
                            % kind) with mr as the argument.
                            hn(3) = exp(1i*mr)/(1i*mr); % m = 0
                            
                            m = 0;
                            P(3) = 1; % P value for m = 0
                            dh = computeDH(hf(2:3),m,mu); % m = 0
                            Bm = computeBm(m,dh,hn(3)); % m = 0
                            psiPS = psiPS + conj(Bm)*P(3);
                            termIndx = 1;
                            psiPSVec(termIndx) = psiPS;
                            lF = true;
                            P = circshift(P,-1);
                            P(3) = cosTheta;
                            hn = circshift(hn,-1);
                            hn(3) = hn(2)*((1/mr)-1i); % m = 1
                            while lF
                                m = m + 1;
                                hf = circshift(hf,-1);
                                hf(3) = computeHX(hf(1:2),m+1,mu);
                                dh = computeDH(hf(2:3),m,mu);
                                Bm = computeBm(m,dh,hn(3));
                                psiPS = psiPS + conj(Bm)*P(3);
                                termIndx = termIndx + 1;
                                psiPSVec(termIndx) = psiPS;
                                if termIndx == (numChkTerms+1)
                                    metricVec = zeros(numChkTerms,1);
                                    for zz = 1:numChkTerms
                                        metricVec(zz) = psiPSVec(zz+1)-...
                                            psiPSVec(zz);
                                    end
                                    realVec = round(real(metricVec),p);
                                    imagVec = round(imag(metricVec),p);
                                    lF1 = any([realVec;imagVec] ~= 0);
                                    lF2 = ~any(isnan([realVec;imagVec]));
                                    lF = lF1 && lF2;
                                    termIndx = termIndx - 1;
                                    psiPSVec = circshift(psiPSVec,-1);
                                else
                                    lF = true;
                                end
                                P = circshift(P,-1);
                                P(3) = computeP(P(1:2),m+1,cosTheta);
                                hn = circshift(hn,-1);
                                hn(3) = computeHX(hn(1:2),m+1,mr);
                            end
                            psiPSVec = circshift(psiPSVec,1);
                        end
                        
                        tStart = tic;
                        for qq = 1:1000
                            % Initialize hf to store 3 most recent terms
                            hf = zeros(3,1);
                            % Initialize P to store 3 most recent terms
                            P = zeros(3,1);
                            hf(2) = exp(1i*mu)/(1i*mu); % m = 0
                            hf(3) = hf(2)*((1/mu)-1i); % m = 1
                            
                            % Initialize sums
                            psiPS = 0;
                            psiPSVec = zeros(numChkTerms+1,1);
                            
                            % Init. hn to store 3 most recent terms
                            hn = zeros(3,1);
                            % Note: hn is the variable storing the value of
                            % the spherical-Hankel function (of the first 
                            % kind) with mr as the argument.
                            hn(3) = exp(1i*mr)/(1i*mr); % m = 0
                            
                            m = 0;
                            P(3) = 1; % P value for m = 0
                            dh = computeDH(hf(2:3),m,mu); % m = 0
                            Bm = computeBm(m,dh,hn(3)); % m = 0
                            psiPS = psiPS + conj(Bm)*P(3);
                            termIndx = 1;
                            psiPSVec(termIndx) = psiPS;
                            lF = true;
                            P = circshift(P,-1);
                            P(3) = cosTheta;
                            hn = circshift(hn,-1);
                            hn(3) = hn(2)*((1/mr)-1i); % m = 1
                            while lF
                                m = m + 1;
                                hf = circshift(hf,-1);
                                hf(3) = computeHX(hf(1:2),m+1,mu);
                                dh = computeDH(hf(2:3),m,mu);
                                Bm = computeBm(m,dh,hn(3));
                                psiPS = psiPS + conj(Bm)*P(3);
                                termIndx = termIndx + 1;
                                psiPSVec(termIndx) = psiPS;
                                if termIndx == (numChkTerms+1)
                                    metricVec = zeros(numChkTerms,1);
                                    for zz = 1:numChkTerms
                                        metricVec(zz) = psiPSVec(zz+1)-...
                                            psiPSVec(zz);
                                    end
                                    realVec = round(real(metricVec),p);
                                    imagVec = round(imag(metricVec),p);
                                    lF1 = any([realVec;imagVec] ~= 0);
                                    lF2 = ~any(isnan([realVec;imagVec]));
                                    lF = lF1 && lF2;
                                    termIndx = termIndx - 1;
                                    psiPSVec = circshift(psiPSVec,-1);
                                else
                                    lF = true;
                                end
                                P = circshift(P,-1);
                                P(3) = computeP(P(1:2),m+1,cosTheta);
                                hn = circshift(hn,-1);
                                hn(3) = computeHX(hn(1:2),m+1,mr);
                            end
                            psiPSVec = circshift(psiPSVec,1);
                        end
                        tStop = toc(tStart);
                        H_iter{pp,ii}(kk,jj) = -(rho/mu)*exp(1i*mr)*...
                            psiPSVec(1);
                        
                        CT_iter{pp,ii}(kk,jj) = tStop/1000;
                    end
                    N_iter{pp,ii}(kk,jj) = m-termIndx;
                end
            end
        end
    end
    
    % Clear variables generated only by this cell
    varsafter = who;
    newvars = setdiff(varsafter,varsbefore);
    varstokeep = {'H_iter','CT_iter','N_iter'};
    varstoremove = setdiff(newvars,varstokeep);
    clear(varstoremove{:});
    clearvars varstokeep varstoremove varsafter newvars
end

%% Perform and time proposed formula-based RS-HRTF computation
%
% This uses order, Nmin, determined using the formula given by Sridhar and 
% Choueiri [1].

if ~loadFlag
    % Get existing variables defined in workspace
    varsbefore = who;
    
    H_direct = cell(numPs,numThetas);
    CT_direct = cell(numPs,numThetas);
    for pp = 1:numPs
        p = pVec(pp);
        for ii = 1:numThetas
            H_direct{pp,ii} = zeros(numMus,numRhos);
            CT_direct{pp,ii} = zeros(numMus,numRhos);
            cosTheta = cosThetaVec(ii);
            for jj = 1:numRhos
                for kk = 1:numMus
                    rho = rhoMat(kk,jj);
                    mu = muMat(kk,jj);
                    mr = muRhoMat(kk,jj);
                    % Although we don't use p here, the order, N, that we
                    % use depends on p
                    fprintf(['Computing for p = %d, theta = %d, rho =',...
                        ' %d, and mu = %3.1f...\n'],p,thetaVec(ii),rho,...
                        round(mu,1));
                    
                    if rho == inf
                        % Perform 10 untimed repetitions
                        alphaVal = 1.41;
                        betaVal = 2.73;
                        gammaVal = 0.8;
                        N = round(alphaVal+(betaVal*mu.^gammaVal));
                        for qq = 1:10
                            % Initialize hf to store 3 most recent terms
                            hf = zeros(3,1);
                            % Initialize P to store 3 most recent terms
                            P = zeros(3,1);
                            hf(2) = exp(1i*mu)/(1i*mu); % m = 0
                            hf(3) = hf(2)*((1/mu)-1i); % m = 1
                            
                            % Initialize sums
                            psiPS = 0;
                            
                            m = 0;
                            P(3) = 1; % P value for m = 0
                            dh = computeDH(hf(2:3),m,mu); % m = 0
                            Am = computeAm(m,dh); % m = 0
                            psiPS = psiPS + conj(Am)*P(3);
                            P = circshift(P,-1);
                            P(3) = cosTheta;
                            for m = 1:N
                                hf = circshift(hf,-1);
                                hf(3) = computeHX(hf(1:2),m+1,mu);
                                dh = computeDH(hf(2:3),m,mu);
                                Am = computeAm(m,dh);
                                psiPS = psiPS + conj(Am)*P(3);
                                P = circshift(P,-1);
                                P(3) = computeP(P(1:2),m+1,cosTheta);
                            end
                        end
                        
                        % Perform 1000 timed repetitions
                        tStart = tic;
                        alphaVal = 1.41;
                        betaVal = 2.73;
                        gammaVal = 0.8;
                        N = round(alphaVal+(betaVal*mu.^gammaVal));
                        for qq = 1:1000
                            hf = zeros(3,1);
                            P = zeros(3,1);
                            hf(2) = exp(1i*mu)/(1i*mu);
                            hf(3) = hf(2)*((1/mu)-1i);
                            
                            psiPS = 0;
                            
                            m = 0;
                            P(3) = 1;
                            dh = computeDH(hf(2:3),m,mu);
                            Am = computeAm(m,dh);
                            psiPS = psiPS + conj(Am)*P(3);
                            P = circshift(P,-1);
                            P(3) = cosTheta;
                            for m = 1:N
                                hf = circshift(hf,-1);
                                hf(3) = computeHX(hf(1:2),m+1,mu);
                                dh = computeDH(hf(2:3),m,mu);
                                Am = computeAm(m,dh);
                                psiPS = psiPS + conj(Am)*P(3);
                                P = circshift(P,-1);
                                P(3) = computeP(P(1:2),m+1,cosTheta);
                            end
                        end
                        tStop = toc(tStart);
                        H_direct{pp,ii}(kk,jj) = (1/mu^2)*psiPS;
                        
                        CT_direct{pp,ii}(kk,jj) = tStop/1000;
                    else
                        alphaVal = (1.41*rho+3.9)/(rho-1.36);
                        betaVal = (2.73*rho-4.75)/(rho-1.21);
                        gammaVal = (0.8*rho-1.01)/(rho-1.45);
                        N = round(alphaVal+(betaVal*mu.^gammaVal));
                        for qq = 1:10
                            hf = zeros(3,1);
                            P = zeros(3,1);
                            hf(2) = exp(1i*mu)/(1i*mu);
                            hf(3) = hf(2)*((1/mu)-1i);
                            
                            % Initialize sums
                            psiPS = 0;
                            
                            % Init. hn to store 3 most recent terms
                            hn = zeros(3,1);
                            % Note: hn is the variable storing the value of
                            % the spherical-Hankel function (of the first 
                            % kind) with mr as the argument.
                            hn(3) = exp(1i*mr)/(1i*mr); % m = 0
                            
                            m = 0;
                            P(3) = 1; % P for m = 0
                            dh = computeDH(hf(2:3),m,mu); % m = 0
                            Bm = computeBm(m,dh,hn(3)); % m = 0
                            psiPS = psiPS + conj(Bm)*P(3);
                            P = circshift(P,-1);
                            P(3) = cosTheta;
                            hn = circshift(hn,-1);
                            hn(3) = hn(2)*((1/mr)-1i); % m = 1
                            for m = 1:N
                                hf = circshift(hf,-1);
                                hf(3) = computeHX(hf(1:2),m+1,mu);
                                dh = computeDH(hf(2:3),m,mu);
                                Bm = computeBm(m,dh,hn(3));
                                psiPS = psiPS + conj(Bm)*P(3);
                                P = circshift(P,-1);
                                P(3) = computeP(P(1:2),m+1,cosTheta);
                                hn = circshift(hn,-1);
                                hn(3) = computeHX(hn(1:2),m+1,mr);
                            end
                        end
                        
                        tStart = tic;
                        alphaVal = (1.41*rho+3.9)/(rho-1.36);
                        betaVal = (2.73*rho-4.75)/(rho-1.21);
                        gammaVal = (0.8*rho-1.01)/(rho-1.45);
                        N = round(alphaVal+(betaVal*mu.^gammaVal));
                        for qq = 1:1000
                            hf = zeros(3,1);
                            P = zeros(3,1);
                            hf(2) = exp(1i*mu)/(1i*mu);
                            hf(3) = hf(2)*((1/mu)-1i);
                            
                            psiPS = 0;
                            
                            hn = zeros(3,1);
                            hn(3) = exp(1i*mr)/(1i*mr);
                            
                            m = 0;
                            P(3) = 1;
                            dh = computeDH(hf(2:3),m,mu);
                            Bm = computeBm(m,dh,hn(3));
                            psiPS = psiPS + conj(Bm)*P(3);
                            P = circshift(P,-1);
                            P(3) = cosTheta;
                            hn = circshift(hn,-1);
                            hn(3) = hn(2)*((1/mr)-1i);
                            for m = 1:N
                                hf = circshift(hf,-1);
                                hf(3) = computeHX(hf(1:2),m+1,mu);
                                dh = computeDH(hf(2:3),m,mu);
                                Bm = computeBm(m,dh,hn(3));
                                psiPS = psiPS + conj(Bm)*P(3);
                                P = circshift(P,-1);
                                P(3) = computeP(P(1:2),m+1,cosTheta);
                                hn = circshift(hn,-1);
                                hn(3) = computeHX(hn(1:2),m+1,mr);
                            end
                        end
                        tStop = toc(tStart);
                        H_direct{pp,ii}(kk,jj) = -(rho/mu)*exp(1i*mr)*...
                            psiPS;
                        
                        CT_direct{pp,ii}(kk,jj) = tStop/1000;
                    end
                end
            end
        end
    end
    
    % Clear variables generated only by this cell
    varsafter = who;
    newvars = setdiff(varsafter,varsbefore);
    varstokeep = {'H_direct','CT_direct'};
    varstoremove = setdiff(newvars,varstokeep);
    clear(varstoremove{:});
    clearvars varstokeep varstoremove varsafter newvars
end

%% Perform and time rule-of-thumb RS-HRTF computation
%
% This is a direct computation (like the one in the previous cell) but uses 
% the order, N, determined using the rule-of-thumb formula given by 
% Gumerov and Duraiswami. For more information on the latter formula, see
% the paper by Sridhar and Choueiri [1].

if ~loadFlag
    % Get existing variables defined in workspace
    varsbefore = who;
    
    H_gd = cell(numPs,numThetas);
    CT_gd = cell(numPs,numThetas);
    for pp = 1:numPs
        p = pVec(pp);
        for ii = 1:numThetas
            H_gd{pp,ii} = zeros(numMus,numRhos);
            CT_gd{pp,ii} = zeros(numMus,numRhos);
            cosTheta = cosThetaVec(ii);
            for jj = 1:numRhos
                for kk = 1:numMus
                    rho = rhoMat(kk,jj);
                    mu = muMat(kk,jj);
                    mr = muRhoMat(kk,jj);
                    % Although we don't use p here, the order, N, that we 
                    % use depends on p
                    fprintf(['Computing for p = %d, theta = %d, rho =',...
                        ' %d, and mu = %3.1f...\n'],p,thetaVec(ii),rho,...
                        round(mu,1));
                    
                    if rho == inf
                        % Perform 10 untimed repetitions
                        N = round(exp(1)*mu);
                        for qq = 1:10
                            % Initialize hf to store 3 most recent terms
                            hf = zeros(3,1);
                            % Initialize P to store 3 most recent terms
                            P = zeros(3,1);
                            hf(2) = exp(1i*mu)/(1i*mu); % m = 0
                            hf(3) = hf(2)*((1/mu)-1i); % m = 1
                            
                            % Initialize sums
                            psiPS = 0;
                            
                            m = 0;
                            P(3) = 1; % P value for m = 0
                            dh = computeDH(hf(2:3),m,mu); % m = 0
                            Am = computeAm(m,dh); % m = 0
                            psiPS = psiPS + conj(Am)*P(3);
                            P = circshift(P,-1);
                            P(3) = cosTheta;
                            for m = 1:N
                                hf = circshift(hf,-1);
                                hf(3) = computeHX(hf(1:2),m+1,mu);
                                dh = computeDH(hf(2:3),m,mu);
                                Am = computeAm(m,dh);
                                psiPS = psiPS + conj(Am)*P(3);
                                P = circshift(P,-1);
                                P(3) = computeP(P(1:2),m+1,cosTheta);
                            end
                        end
                        
                        % Perform 1000 timed repetitions
                        tStart = tic;
                        N = round(exp(1)*mu);
                        for qq = 1:1000
                            hf = zeros(3,1);
                            P = zeros(3,1);
                            hf(2) = exp(1i*mu)/(1i*mu);
                            hf(3) = hf(2)*((1/mu)-1i);
                            
                            psiPS = 0;
                            
                            m = 0;
                            P(3) = 1;
                            dh = computeDH(hf(2:3),m,mu);
                            Am = computeAm(m,dh);
                            psiPS = psiPS + conj(Am)*P(3);
                            P = circshift(P,-1);
                            P(3) = cosTheta;
                            for m = 1:N
                                hf = circshift(hf,-1);
                                hf(3) = computeHX(hf(1:2),m+1,mu);
                                dh = computeDH(hf(2:3),m,mu);
                                Am = computeAm(m,dh);
                                psiPS = psiPS + conj(Am)*P(3);
                                P = circshift(P,-1);
                                P(3) = computeP(P(1:2),m+1,cosTheta);
                            end
                        end
                        tStop = toc(tStart);
                        H_gd{pp,ii}(kk,jj) = (1/mu^2)*psiPS;
                        
                        CT_gd{pp,ii}(kk,jj) = tStop/1000;
                    else
                        N = round(exp(1)*mu);
                        for qq = 1:10
                            hf = zeros(3,1);
                            P = zeros(3,1);
                            hf(2) = exp(1i*mu)/(1i*mu);
                            hf(3) = hf(2)*((1/mu)-1i);
                            
                            % Initialize sums
                            psiPS = 0;
                            
                            % Init. hn to store 3 most recent terms
                            hn = zeros(3,1);
                            % Note: hn is the variable storing the value of
                            % the spherical-Hankel function (of the first 
                            % kind) with mr as the argument.
                            hn(3) = exp(1i*mr)/(1i*mr); % m = 0
                            
                            m = 0;
                            P(3) = 1; % P for m = 0
                            dh = computeDH(hf(2:3),m,mu); % m = 0
                            Bm = computeBm(m,dh,hn(3)); % m = 0
                            psiPS = psiPS + conj(Bm)*P(3);
                            P = circshift(P,-1);
                            P(3) = cosTheta;
                            hn = circshift(hn,-1);
                            hn(3) = hn(2)*((1/mr)-1i); % m = 1
                            for m = 1:N
                                hf = circshift(hf,-1);
                                hf(3) = computeHX(hf(1:2),m+1,mu);
                                dh = computeDH(hf(2:3),m,mu);
                                Bm = computeBm(m,dh,hn(3));
                                psiPS = psiPS + conj(Bm)*P(3);
                                P = circshift(P,-1);
                                P(3) = computeP(P(1:2),m+1,cosTheta);
                                hn = circshift(hn,-1);
                                hn(3) = computeHX(hn(1:2),m+1,mr);
                            end
                        end
                        
                        tStart = tic;
                        N = round(exp(1)*mu);
                        for qq = 1:1000
                            hf = zeros(3,1);
                            P = zeros(3,1);
                            hf(2) = exp(1i*mu)/(1i*mu);
                            hf(3) = hf(2)*((1/mu)-1i);
                            
                            psiPS = 0;
                            
                            hn = zeros(3,1);
                            hn(3) = exp(1i*mr)/(1i*mr);
                            
                            m = 0;
                            P(3) = 1;
                            dh = computeDH(hf(2:3),m,mu);
                            Bm = computeBm(m,dh,hn(3));
                            psiPS = psiPS + conj(Bm)*P(3);
                            P = circshift(P,-1);
                            P(3) = cosTheta;
                            hn = circshift(hn,-1);
                            hn(3) = hn(2)*((1/mr)-1i);
                            for m = 1:N
                                hf = circshift(hf,-1);
                                hf(3) = computeHX(hf(1:2),m+1,mu);
                                dh = computeDH(hf(2:3),m,mu);
                                Bm = computeBm(m,dh,hn(3));
                                psiPS = psiPS + conj(Bm)*P(3);
                                P = circshift(P,-1);
                                P(3) = computeP(P(1:2),m+1,cosTheta);
                                hn = circshift(hn,-1);
                                hn(3) = computeHX(hn(1:2),m+1,mr);
                            end
                        end
                        tStop = toc(tStart);
                        H_gd{pp,ii}(kk,jj) = -(rho/mu)*exp(1i*mr)*psiPS;
                        
                        CT_gd{pp,ii}(kk,jj) = tStop/1000;
                    end
                end
            end
        end
    end
    
    % Clear variables generated only by this cell
    varsafter = who;
    newvars = setdiff(varsafter,varsbefore);
    varstokeep = {'H_gd','CT_gd'};
    varstoremove = setdiff(newvars,varstokeep);
    clear(varstoremove{:});
    clearvars varstokeep varstoremove varsafter newvars
end

%% Compare computation times

if ~loadFlag
    % Iterative vs direct
    rctMat = cell(numPs,numThetas); % Relative computation time
    rsiMat = cell(numPs,numThetas); % Relative speed increase
    rctMat_gd = cell(numPs,numThetas);
    rsiMat_gd = cell(numPs,numThetas);
    rctMat_gd2 = cell(numPs,numThetas);
    rsiMat_gd2 = cell(numPs,numThetas);
    meanRSIMat = cell(numPs,1);
    minRSIMat = cell(numPs,1);
    meanRSIMat_gd = cell(numPs,1);
    minRSIMat_gd = cell(numPs,1);
    meanRSIMat_gd2 = cell(numPs,1);
    minRSIMat_gd2 = cell(numPs,1);
    globalMeanRSIs = zeros(numPs,1);
    globalMinRSIs = zeros(numPs,1);
    globalMeanRSIs_gd = zeros(numPs,1);
    globalMinRSIs_gd = zeros(numPs,1);
    globalMeanRSIs_gd2 = zeros(numPs,1);
    globalMinRSIs_gd2 = zeros(numPs,1);
    diffVec1 = cell(numPs*numThetas,1);
    diffVec2 = cell(numPs*numThetas,1);
    diffVec3 = cell(numPs*numThetas,1);
    indx = 1;
    for pp = 1:numPs
        meanRSIMat{pp,1} = zeros(numRhos,numThetas);
        minRSIMat{pp,1} = zeros(numRhos,numThetas);
        meanRSIMat_gd{pp,1} = zeros(numRhos,numThetas);
        minRSIMat_gd{pp,1} = zeros(numRhos,numThetas);
        meanRSIMat_gd2{pp,1} = zeros(numRhos,numThetas);
        minRSIMat_gd2{pp,1} = zeros(numRhos,numThetas);
        for ii = 1:numThetas
            rctMat{pp,ii} = CT_direct{pp,ii}./CT_iter{pp,ii};
            rsiMat{pp,ii} = 100*(CT_iter{pp,ii}-CT_direct{pp,ii})./...
                CT_direct{pp,ii};
            meanRSIMat{pp,1}(:,ii) = mean(rsiMat{pp,ii}).';
            minRSIMat{pp,1}(:,ii) = min(rsiMat{pp,ii}).';
            
            rctMat_gd{pp,ii} = CT_gd{pp,ii}./CT_iter{pp,ii};
            rsiMat_gd{pp,ii} = 100*(CT_iter{pp,ii}-CT_gd{pp,ii})./...
                CT_gd{pp,ii};
            meanRSIMat_gd{pp,1}(:,ii) = mean(rsiMat_gd{pp,ii}).';
            minRSIMat_gd{pp,1}(:,ii) = min(rsiMat_gd{pp,ii}).';
            
            rctMat_gd2{pp,ii} = CT_direct{pp,ii}./CT_gd{pp,ii};
            rsiMat_gd2{pp,ii} = 100*(CT_gd{pp,ii}-CT_direct{pp,ii})./...
                CT_direct{pp,ii};
            meanRSIMat_gd2{pp,1}(:,ii) = mean(rsiMat_gd2{pp,ii}).';
            minRSIMat_gd2{pp,1}(:,ii) = min(rsiMat_gd2{pp,ii}).';
            
            diffVec1{indx,1} = reshape(CT_iter{pp,ii}-CT_direct{pp,ii},...
                numMus*numRhos,1);
            diffVec2{indx,1} = reshape(CT_iter{pp,ii}-CT_gd{pp,ii},...
                numMus*numRhos,1);
            diffVec3{indx,1} = reshape(CT_gd{pp,ii}-CT_direct{pp,ii},...
                numMus*numRhos,1);
            indx = indx + 1;
        end
        globalMeanRSIs(pp) = mean(mean(meanRSIMat{pp,1}));
        globalMinRSIs(pp) = min(min(minRSIMat{pp,1}));
        globalMeanRSIs_gd(pp) = mean(mean(meanRSIMat_gd{pp,1}));
        globalMinRSIs_gd(pp) = min(min(minRSIMat_gd{pp,1}));
        globalMeanRSIs_gd2(pp) = mean(mean(meanRSIMat_gd2{pp,1}));
        globalMinRSIs_gd2(pp) = min(min(minRSIMat_gd2{pp,1}));
    end
    diffVec1 = cell2mat(diffVec1);
    diffVec2 = cell2mat(diffVec2);
    diffVec3 = cell2mat(diffVec3);
    
    clearvars pp ii
end

%% Perform paired-sample t-tests
%
% The results generated in this cell are not reported in the paper by
% Sridhar and Choueiri [1].

% [h1,p1,ci1] = ttest(diffVec1,0,'Alpha',0.0001,'Tail','Right');
% [h2,p2,ci2] = ttest(diffVec2,0,'Alpha',0.0001,'Tail','Right');
% [h3,p3,ci3] = ttest(diffVec3,0,'Alpha',0.0001,'Tail','Right');
% [p4,h4,ci4] = signrank(diffVec3,0,'alpha',1e-116,'tail','right');

%% Comments pertaining to results in Sridhar and Choueiri [1].
%
% The RSI values published in Sec. 4 of the paper are calculated as
% mean(globalMeanRSIs) and mean(globalMeanRSIs_gd2).

%% (Optional) Export data

exportFlag = false;

if ~loadFlag && exportFlag
    exportFileName = 'sc2020formula_vs_iter_vs_gd2002'; % CHANGE AS NEEDED
    
    [~,exportFileName,~] = fileparts(exportFileName); % Removes extension
    save(fullfile(dataPath,'Computation',[exportFileName,'.mat']),...
        '-regexp',['^(?!(dataPath|plotsPath|projectPath|',...
        'exportFileName|exportFlag|loadFlag)$).'],'-v7.3');
end

%% Define functions required by script
%
% Do NOT modify the contents of this cell.

function out = computeHX(hx,m,z)
%COMPUTEHX Spherical-Hankel function of the first kind.
%   Y = COMPUTEHX(A,M,Z) computes the value, Y, of the spherical-hankel
%   function of the first kind of order M and with argument Z. A is a
%   2-element vector containing the values of the function for the same
%   argument but for orders M-2 and M-1, stored in this order.

out = (((2*m-1)/z)*hx(2))-hx(1);

end

function out = computeDH(hf,m,z)
%COMPUTEDH Derivative of the spherical-Hankel function of the first kind.
%   Y = COMPUTEDH(A,M,Z) computes the value, Y, of the derivative of the 
%   spherical-hankel function of the first kind of order M and with 
%   argument Z. A is a 2-element vector containing the values of the 
%   spherical-hankel function of the first kind for the same argument, Z,
%   but for orders M and M+1, stored in this order.

out = -hf(2)+((m/z)*hf(1));

end

function out = computeP(P,m,x)
%COMPUTEP Legendre polynomial
%   Y = COMPUTEP(P,M,X) computes the value, Y, of the Legendre polynomial
%   of degree M and with argument X. P is a 2-element vector containing the 
%   values of the Legendre polynomial for the same argument but for degrees
%   M-2 and M-1, stored in this order.

out = (((2*m)-1)/m)*x*P(2)-(((m-1)/m)*P(1));

end

function out = computeAm(m,hp)
%COMPUTEAM Factor in summation term in solution corresponding to infinitely 
%distant source. See the paper by Sridhar and Choueiri [1].

out = ((-1i)^(m-1))*(2*m+1)/hp;

end

function out = computeBm(m,hp,hn)
%COMPUTEBM Factor in summation term in solution corresponding to source at 
%a finite distance. See the paper by Sridhar and Choueiri [1].

out = ((2*m)+1)*hn/hp;

end
