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
    load(fullfile(dataPath,'Computation','cb1980_vs_direct_vs_gd2002.mat'))
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
    TcVec = 10.^(-10:2:-2).';
    numTcs = length(TcVec);
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
% This is based on the algorithm by Cooper and Bauck. For more information, 
% see the publication by Sridhar and Choueiri [1].

if ~loadFlag
    % Get existing variables defined in workspace
    varsbefore = who;
    
    H_iter = cell(numTcs,numThetas);
    CT_iter = cell(numTcs,numThetas);
    N_iter = cell(numTcs,numThetas);
    for pp = 1:numTcs
        Tc = TcVec(pp);
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
                    fprintf(['Computing for Tc = 10e%d, theta = %d,',...
                        ' rho = %d, and mu = %3.1f...\n'],...
                        round(log10(Tc)),thetaVec(ii),rho,round(mu,1));
                    
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
                            psiPSVec = zeros(3,1);
                            
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
                                if termIndx == 3
                                    metricVec = zeros(2,1);
                                    metricVec(1) = abs((psiPSVec(2)-...
                                        psiPSVec(1))/psiPSVec(2));
                                    metricVec(2) = abs((psiPSVec(3)-...
                                        psiPSVec(2))/psiPSVec(3));
                                    lF = any(metricVec > Tc);
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
                            hf = zeros(3,1);
                            P = zeros(3,1);
                            hf(2) = exp(1i*mu)/(1i*mu);
                            hf(3) = hf(2)*((1/mu)-1i);
                            
                            psiPS = 0;
                            psiPSVec = zeros(3,1);
                            
                            m = 0;
                            P(3) = 1;
                            dh = computeDH(hf(2:3),m,mu);
                            Am = computeAm(m,dh);
                            psiPS = psiPS + conj(Am)*P(3);
                            termIndx = 1;
                            psiPSVec(termIndx) = psiPS;
                            lF = true;
                            P = circshift(P,-1);
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
                                if termIndx == 3
                                    metricVec = zeros(2,1);
                                    metricVec(1) = abs((psiPSVec(2)-...
                                        psiPSVec(1))/psiPSVec(2));
                                    metricVec(2) = abs((psiPSVec(3)-...
                                        psiPSVec(2))/psiPSVec(3));
                                    lF = any(metricVec > Tc);
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
                            psiPSVec = zeros(3,1);
                            
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
                                if termIndx == 3
                                    metricVec = zeros(2,1);
                                    metricVec(1) = abs((psiPSVec(2)-...
                                        psiPSVec(1))/psiPSVec(2));
                                    metricVec(2) = abs((psiPSVec(3)-...
                                        psiPSVec(2))/psiPSVec(3));
                                    lF = any(metricVec > Tc);
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
                            hf = zeros(3,1);
                            P = zeros(3,1);
                            hf(2) = exp(1i*mu)/(1i*mu);
                            hf(3) = hf(2)*((1/mu)-1i);
                            
                            psiPS = 0;
                            psiPSVec = zeros(3,1);
                            
                            hn = zeros(3,1);
                            hn(3) = exp(1i*mr)/(1i*mr);
                            
                            m = 0;
                            P(3) = 1;
                            dh = computeDH(hf(2:3),m,mu);
                            Bm = computeBm(m,dh,hn(3));
                            psiPS = psiPS + conj(Bm)*P(3);
                            termIndx = 1;
                            psiPSVec(termIndx) = psiPS;
                            lF = true;
                            P = circshift(P,-1);
                            P(3) = cosTheta;
                            hn = circshift(hn,-1);
                            hn(3) = hn(2)*((1/mr)-1i);
                            while lF
                                m = m + 1;
                                hf = circshift(hf,-1);
                                hf(3) = computeHX(hf(1:2),m+1,mu);
                                dh = computeDH(hf(2:3),m,mu);
                                Bm = computeBm(m,dh,hn(3));
                                psiPS = psiPS + conj(Bm)*P(3);
                                termIndx = termIndx + 1;
                                psiPSVec(termIndx) = psiPS;
                                if termIndx == 3
                                    metricVec = zeros(2,1);
                                    metricVec(1) = abs((psiPSVec(2)-...
                                        psiPSVec(1))/psiPSVec(2));
                                    metricVec(2) = abs((psiPSVec(3)-...
                                        psiPSVec(2))/psiPSVec(3));
                                    lF = any(metricVec > Tc);
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

%% Perform and time direct RS-HRTF computation
%
% This uses the order, N, determined in the iterative calculation above.

if ~loadFlag
    % Get existing variables defined in workspace
    varsbefore = who;
    
    H_direct = cell(numTcs,numThetas);
    CT_direct = cell(numTcs,numThetas);
    for pp = 1:numTcs
        Tc = TcVec(pp);
        for ii = 1:numThetas
            H_direct{pp,ii} = zeros(numMus,numRhos);
            CT_direct{pp,ii} = zeros(numMus,numRhos);
            cosTheta = cosThetaVec(ii);
            for jj = 1:numRhos
                for kk = 1:numMus
                    rho = rhoMat(kk,jj);
                    mu = muMat(kk,jj);
                    mr = muRhoMat(kk,jj);
                    % Although we don't use Tc here, the order, N, that we 
                    % use depends on Tc
                    fprintf(['Computing for Tc = 10e%d, theta = %d,',...
                        ' rho = %d, and mu = %3.1f...\n'],...
                        round(log10(Tc)),thetaVec(ii),rho,round(mu,1));
                    
                    if rho == inf
                        % Perform 10 untimed repetitions
                        N = N_iter{pp,ii}(kk,jj);
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
                        N = N_iter{pp,ii}(kk,jj);
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
                        N = N_iter{pp,ii}(kk,jj);
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
                        N = N_iter{pp,ii}(kk,jj);
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
    
    rho = inf;
    H_gd = zeros(numMus,numThetas);
    CT_gd = zeros(numMus,numThetas);
    for ii = 1:numThetas
        cosTheta = cosThetaVec(ii);
        for kk = 1:numMus
            mu = muMat(kk,1);
            fprintf('Computing for theta = %d and mu = %3.1f...\n',...
                thetaVec(ii),round(mu,1));
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
                    hf(3) = computeHX(hf(1:2),m+1,mu); % m+1, not m
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
            H_gd(kk,ii) = (1/mu^2)*psiPS;
            
            CT_gd(kk,ii) = tStop/1000;
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
    rctMat = cell(numTcs,numThetas); % Relative computation time
    rsiMat = cell(numTcs,numThetas); % Relative speed increase
    meanRSIMat = cell(numTcs,1);
    minRSIMat = cell(numTcs,1);
    globalMeanRSIs = zeros(numTcs,1);
    globalMinRSIs = zeros(numTcs,1);
    diffVec1 = cell(numTcs*numThetas,1);
    indx = 1;
    for pp = 1:numTcs
        meanRSIMat{pp,1} = zeros(numRhos,numThetas);
        minRSIMat{pp,1} = zeros(numRhos,numThetas);
        for ii = 1:numThetas
            rctMat{pp,ii} = CT_direct{pp,ii}./CT_iter{pp,ii};
            rsiMat{pp,ii} = 100*(CT_iter{pp,ii}-CT_direct{pp,ii})./...
                CT_direct{pp,ii};
            meanRSIMat{pp,1}(:,ii) = mean(rsiMat{pp,ii}).';
            minRSIMat{pp,1}(:,ii) = min(rsiMat{pp,ii}).';
            diffVec1{indx,1} = reshape(CT_iter{pp,ii}-CT_direct{pp,ii},...
                numMus*numRhos,1);
            indx = indx + 1;
        end
        globalMeanRSIs(pp) = mean(mean(meanRSIMat{pp,1}));
        globalMinRSIs(pp) = min(min(minRSIMat{pp,1}));
    end
    diffVec1 = cell2mat(diffVec1);
    
    % Iterative vs gd2002
    Tc_indx = find(TcVec == 10^(-6),1);
    rho_indx = find(rhoVec == inf,1);
    if isempty(Tc_indx) || isempty(rho_indx)
        error(['Iterative calculation for T_c = 10^(-6) and rho = inf',...
            ' not found.'])
    end
    rctMat_gd = zeros(numMus,numThetas); % Relative computation time
    rsiMat_gd = zeros(numMus,numThetas); % Relative speed increase
    diffVec2 = cell(numThetas,1);
    for ii = 1:numThetas
        CT_iter_curr = CT_iter{Tc_indx,ii}(:,rho_indx);
        CT_gd_curr = CT_gd(:,ii);
        rctMat_gd(:,ii) = CT_gd_curr./CT_iter_curr;
        rsiMat_gd(:,ii) = 100*(CT_iter_curr-CT_gd_curr)./CT_gd_curr;
        diffVec2{ii,1} = CT_iter_curr-CT_gd_curr;
    end
    meanRSIMat_gd = mean(rsiMat_gd);
    diffVec2 = cell2mat(diffVec2);
    
    clearvars pp ii
end

%% Comments pertaining to results in Sridhar and Choueiri [1].
%
% The RSI values published in Sec. 2.3 of the paper are calculated as
% mean(globalMeanRSIs) and mean(meanRSIMat_gd).

%% (Optional) Export data

exportFlag = false;

if ~loadFlag && exportFlag
    exportFileName = 'cb1980_vs_direct_vs_gd2002'; % CHANGE AS NEEDED
    
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
