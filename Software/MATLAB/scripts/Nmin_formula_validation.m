%% Validation of formula for Nmin
%
% In this script, we validate the formula for Nmin derived using the script
% Nmin_formula_derivation.m. The results generated in this script are not
% directly reported in the paper by Sridhar and Choueiri [1].
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

%% Load model parameter data

load(fullfile(dataPath,'Computation','N_formula_model_params.mat'))

%% Compare Nmin computed using a_Mat, b_Mat, and c_Mat matrices to 
% N_optimal

N_err_spec = cell(rVecLen,numPs);
N_err_max = zeros(rVecLen,numPs);
N_err_mean = zeros(rVecLen,numPs);
N_err_numOver = zeros(rVecLen,numPs);
N_err_numUnder = zeros(rVecLen,numPs);
for ii = 1:rVecLen
    for jj = 1:numPs
        N_test = round(a_Mat(ii,jj) + b_Mat(ii,jj)*muVec(2:nyqIndx).^...
            c_Mat(ii,jj));
        N_err_spec{ii,jj} = N_optimal{ii,jj}(2:nyqIndx)-N_test;
        N_err_max(ii,jj) = max(abs(N_err_spec{ii,jj}));
        N_err_mean(ii,jj) = mean(abs(N_err_spec{ii,jj}));
        N_err_numOver(ii,jj) = sum(double(N_err_spec{ii,jj} < 0));
        N_err_numUnder(ii,jj) = sum(double(N_err_spec{ii,jj} > 0));
    end
end

clearvars ii jj

%% Compare Nmin computed starting from p and q vectors to N_optimal
%
% That is, compare Nmin starting from each formula (for a given p) that is 
% a function of rho only.

N_err_spec2 = cell(rVecLen,numPs);
N_err_max2 = zeros(rVecLen,numPs);
N_err_mean2 = zeros(rVecLen,numPs);
N_err_numOver2 = zeros(rVecLen,numPs);
N_err_numUnder2 = zeros(rVecLen,numPs);
for ii = 1:rVecLen
    for jj = 1:numPs
        if rhoVec(ii) == inf
            aTest = p1_a(jj);
            bTest = p1_b(jj);
            cTest = p1_c(jj);
        else
            aTest = (rhoVec(ii)*p1_a(jj)+p2_a(jj))/(rhoVec(ii)+q1_a(jj));
            bTest = (rhoVec(ii)*p1_b(jj)+p2_b(jj))/(rhoVec(ii)+q1_b(jj));
            cTest = (rhoVec(ii)*p1_c(jj)+p2_c(jj))/(rhoVec(ii)+q1_c(jj));
        end
        
        N_test = round(aTest + bTest*muVec(2:nyqIndx).^cTest);
        N_err_spec2{ii,jj} = N_optimal{ii,jj}(2:nyqIndx)-N_test;
        N_err_max2(ii,jj) = max(abs(N_err_spec2{ii,jj}));
        N_err_mean2(ii,jj) = mean(abs(N_err_spec2{ii,jj}));
        N_err_numOver2(ii,jj) = sum(double(N_err_spec2{ii,jj} < 0));
        N_err_numUnder2(ii,jj) = sum(double(N_err_spec2{ii,jj} > 0));
    end
end

clearvars ii jj

%% Compare Nmin starting from formula that is a function of rho and p.

N_err_spec3 = cell(rVecLen,numPs);
N_err_max3 = zeros(rVecLen,numPs);
N_err_mean3 = zeros(rVecLen,numPs);
N_err_numOver3 = zeros(rVecLen,numPs);
N_err_numUnder3 = zeros(rVecLen,numPs);
for ii = 1:rVecLen
    for jj = 1:numPs
        if rhoVec(ii) == inf
            aTest = 0.3442*pVec(jj) - 0.9584;
            bTest = 1.758*exp(0.1033*pVec(jj));
            cTest = 1.017-0.1041*pVec(jj).^0.5086;
        else
            p1_a = 0.3442*pVec(jj) - 0.9584;
            p2_a = 3.369*pVec(jj) - 3.69;
            q1_a = -(0.7605*pVec(jj) + 1.704)./(pVec(jj) + 0.7461);
            
            p1_b = 1.758*exp(0.1033*pVec(jj));
            p2_b = -(2.043*pVec(jj) - 0.1245);
            p3_b = 2.161*pVec(jj) - 1.964;
            q1_b = 0.001728*pVec(jj).^2.977 - 3.088;
            q2_b = 2.679 - 0.0008429*pVec(jj).^3.343;
            
            p1_c = 1.017-0.1041*pVec(jj).^0.5086;
            p2_c = 0.1267*pVec(jj) - 2.859;
            p3_c = -(0.1227*pVec(jj) - 2.335);
            q1_c = -(0.02058*pVec(jj) + 2.968);
            q2_c = 0.06702*pVec(jj) + 2.169;
            
            aTest = (rhoVec(ii)*p1_a+p2_a)/(rhoVec(ii)+q1_a);
            bTest = (rhoVec(ii)^2*p1_b+rhoVec(ii)*p2_b+p3_b)/...
                (rhoVec(ii)^2+rhoVec(ii)*q1_b+q2_b);
            cTest = (rhoVec(ii)^2*p1_c+rhoVec(ii)*p2_c+p3_c)/...
                (rhoVec(ii)^2+rhoVec(ii)*q1_c+q2_c);
        end
        
        N_test = round(aTest + bTest*muVec(2:nyqIndx).^cTest);
        N_err_spec3{ii,jj} = N_optimal{ii,jj}(2:nyqIndx)-N_test;
        N_err_max3(ii,jj) = max(abs(N_err_spec3{ii,jj}));
        N_err_mean3(ii,jj) = mean(abs(N_err_spec3{ii,jj}));
        N_err_numOver3(ii,jj) = sum(double(N_err_spec3{ii,jj} < 0));
        N_err_numUnder3(ii,jj) = sum(double(N_err_spec3{ii,jj} > 0));
    end
end

clearvars ii jj
