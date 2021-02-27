function [source,head,dsp] = batchComputeSphereHRTFs(source,head,dsp)
%BATCHCOMPUTESPHEREHRTFS Batch compute spherical-head HRTFs.
%   [S,H,D] = BATCHCOMPUTESPHEREHRTFS(S,H,D) takes input structures S, H,
%   and D, computes spherical-head HRTFs and stores the output data in the
%   structures S, H, and D before returning them. The structures contain
%   the following optional input variables (not specifying any will result 
%   in the use of default values):
%
%       In 'S':
%           1. r - source distance in meters
%           2. azDel - azimuth spacing in deg. for computing uniform 
%           azimuths between 0 and 360 deg.
%           3. elDel - elevation spacing in deg. for computing uniform 
%           elevations between -90 and 90 deg.
%           4. azVec - vector of azimuth values in deg.
%           5. elVec - vector of elevation values in deg.
%           6. sDirs - source directions in SOFA spherical coordinates
%           7. sPos - source position in SOFA cartesian coordinates
%       Note: If more than one of 2-7 is specified, then the order of 
%       priority (high to low) is as follows: 7 > 6 > 5/4 > 3/2.
%
%       In 'H':
%           1. a - head radius in meters
%           2. eL - left ear direction in SOFA spherical coordinates
%           specified in degrees as [az,el]
%           3. eR - right ear direction.
%       Note: If only one of eL or eR is specified, HRTFs for only that ear
%       are computed. If neither eL or eR are specified, only HRTFs for the 
%       left ear at [100,-10] are computed.
%
%       In 'D':
%           1. method - method and parameters (see the COMPUTESPHEREHRTF 
%           function in the 3D3A MATLAB Toolbox)
%           2. normloc - HRTF normalization option. Can take values of
%           'center' (default), 'ear', 'center_tdc', and 'center_tdc_corr',
%           where 'center_tdc_corr' is associated with a parameter FMAX 
%           as described in the function COMPUTESPHEREHRTF in the 3D3A 
%           MATLAB Toolbox. In this function, FMAX cannot be specified and 
%           is instead computed automatically based on whether fVec or fS
%           (see point 3 below) is specified.
%           3. fVec or fS - vector (fVec) of frequencies, in Hz, at which 
%           to compute HRTFs or sampling rate (fS) in Hz at which to 
%           compute HRIRs. If both fS and fVec are specified, fS is 
%           preferred. If neither is specified, a default value of fS = 48 
%           kHz is used.
%           
%           If fVec is specified, HRTFs are returned in the struct, 'H'.
%           Otherwise, if fS is specified, HRIRs are returned and the
%           following additional input parameters may be specified to
%           control various aspects of the HRIRs.
%           4. T - HRIR duration in seconds (also see point 6 below)
%           5. causalflag - whether or not to apply a time-shift to ensure 
%           causal IRs are returned. Can take values of true or false. The
%           default value is false.
%           6. causalshift - shift (in ms) to apply if causalflag is true.
%           7. pow2flag - whether or not to make HRIR length a power of 2.
%           Can take values of true or false. If true, then T in point 5
%           above corresponds to the approximate HRIR duration in seconds.
%           The default value is false.
%
%   See also COMPUTESPHEREHRTF.

%   =======================================================================
%   This file is part of the 3D3A Research Lite Toolbox.
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

%% Check input variables and assign default values if needed

narginchk(0,3);

if nargin < 3 || isempty(dsp)
    dsp = struct;
end

if nargin < 2 || isempty(head)
    head = struct;
end

if nargin < 1 || isempty(source)
    source = struct;
end

% Sound source parameters

if ~isfield(source,'r')
    source.r = inf; % Source distance in meters.
end

if ~isfield(source,'sPos') && ~isfield(source,'sDirs')
    if ~isfield(source,'azDel') && ~isfield(source,'azVec')
        source.azDel = 15; % Azimuth delta in degrees.
        azVecMax = ceil((360/source.azDel)-1)*source.azDel;
        azVec = (0:source.azDel:azVecMax).';
    elseif isfield(source,'azDel') && ~isfield(source,'azVec')
        azVecMax = ceil((360/source.azDel)-1)*source.azDel;
        azVec = (0:source.azDel:azVecMax).';
    else
        azVec = source.azVec;
    end
    
    if ~isfield(source,'elDel') && ~isfield(source,'elVec')
        source.elDel = 15; % Elevation delta in degrees.
        elVec = (-90:source.elDel:90).';
    elseif isfield(source,'elDel') && ~isfield(source,'elVec')
        elVec = (-90:source.elDel:90).';
    else
        elVec = source.elVec;
    end
end

% Spherical head parameters

if ~isfield(head,'a')
    head.a = 0.09; % Head radius in meters.
end

if ~isfield(head,'eL')
    % Left ear direction (az,el) in SOFA spherical coordinates (in deg.)
    eLVal = [100,-10];
    lFlag = false;
else
    eLVal = head.eL;
    lFlag = true;
end

if ~isfield(head,'eR')
    eRVal = [260,-10]; % Right ear direction.
    rFlag = false;
else
    eRVal = head.eR;
    rFlag = true;
end

% Computation and DSP parameters

if ~isfield(dsp,'method')
    dsp.method = {'sridharchoueiri2020',inf,4};
end

if ~isfield(dsp,'normloc')
    dsp.normloc = 'center';
end

if ~isfield(dsp,'fS') && ~isfield(dsp,'fVec')
    dsp.fS = 48000; % Sampling rate in Hz.
end

if ~isfield(dsp,'T') && isfield(dsp,'fS')
    dsp.T = 0.01; % Approximate HRIR duration in seconds.
end

if ~isfield(dsp,'causalflag') && isfield(dsp,'fS')
    dsp.causalflag = false;
end

if ~isfield(dsp,'pow2flag') && isfield(dsp,'fS')
    dsp.pow2flag = false;
end

if strcmpi(dsp.normloc,'center_tdc_corr')
    if isfield(dsp,'fS')
        FMAX = dsp.fS/2;
    else
        FMAX = max(dsp.fVec);
    end
else
    FMAX = 0; % Dummy value.
end

%% Pre-compute source position and direction matrices

fprintf('Generating position matrix...');

if source.r == inf
    % Equivalent finite source distance for computing and transforming
    % source positions only.
    efsd = 10;
else
    efsd = source.r;
end

if ~isfield(source,'sPos') && ~isfield(source,'sDirs')
    source.Dirs = makePosMat({azVec,elVec,efsd});
elseif ~isfield(source,'sPos') && isfield(source,'sDirs')
    source.Dirs = source.sDirs;
    source = rmfield(source,'sDirs');
else
    source.Dirs = sofaC2sofaS(source.sPos);
    source = rmfield(source,'sPos');
end

source.numDirs = size(source.Dirs,1);
source.Positions = sofaS2sofaC(source.Dirs);
data_type = class(source.Positions);
switch lower(data_type)
    case 'single'
        round_N = abs(floor(log10(eps('single'))))-2;
    otherwise
        round_N = abs(floor(log10(eps('double'))))-2;
end
[~,uniquePosIndxs,~] = unique(round(source.Positions,round_N),'rows',...
    'stable');
source.Positions = source.Positions(uniquePosIndxs,:);
source.Dirs = sofaC2sofaS(source.Positions);
% Update numDirs after removing duplicate directions
source.numDirs = size(source.Dirs,1);

fprintf('done!\n');

%% Perform other pre-computation tasks

if isfield(dsp,'fS')
    retFlag = 0; % Return HRIRs
    
    if dsp.pow2flag
        dsp.irLen = 2^nextpow2(dsp.T * dsp.fS);
    else
        dsp.irLen = dsp.T * dsp.fS;
    end
    fVec = getFreqVec(dsp.fS,dsp.irLen);
    nyquistIndx = ceil((dsp.irLen+1)/2);
else
    retFlag = 1; % Return HRTFs at specific frequencies
    
    fVec = dsp.fVec;
    nyquistIndx = length(fVec);
end
rho = source.r/(head.a);
muVec = 2*pi*fVec(1:nyquistIndx)*(head.a)/getSoundSpeed();

% Compute angles of incidence from source position data
fprintf('Computing angle(s) of incidence...')
S = source.Dirs(:,1:2);
numPos = size(S,1);
[xS,yS,zS] = sofaS2sofaC(S(:,1),S(:,2),ones(numPos,1));
% Compute left and right HRTFs if both ear positions are specified
if lFlag && rFlag 
    head.eL = eLVal;
    head.eR = eRVal;
    [xEL,yEL,zEL] = sofaS2sofaC(head.eL(1),head.eL(2),1);
    [xER,yER,zER] = sofaS2sofaC(head.eR(1),head.eR(2),1);
    theta{1} = getCentralAngle([xS,yS,zS],repmat([xEL,yEL,zEL],numPos,1));
    theta{2} = getCentralAngle([xS,yS,zS],repmat([xER,yER,zER],numPos,1));
% Compute only right HRTF if only right ear position is specified
elseif ~lFlag && rFlag
    head.eR = eRVal;
    [xER,yER,zER] = sofaS2sofaC(head.eR(1),head.eR(2),1);
    theta{1} = getCentralAngle([xS,yS,zS],repmat([xER,yER,zER],numPos,1));
% Compute only left HRTF if only left ear position/no ear positions are
% specified
else
    head.eL = eLVal;
    [xEL,yEL,zEL] = sofaS2sofaC(head.eL(1),head.eL(2),1);
    theta{1} = getCentralAngle([xS,yS,zS],repmat([xEL,yEL,zEL],numPos,1));
end
fprintf('done!\n')

%% Compute spherical head HRTFs

numEars = numel(theta);
% The formula-based methods are applicable for rho >= 2 only. For rho < 2, 
% the formulas with rho = 2 may be used to obtain satisfactory results.
if rho < 2
    rho = 2;
end

if ~retFlag
    h = cell(numEars,1); 
end
H = cell(numEars,1);
N = cell(numEars,1);
TH = cell(numEars,1);
pL1 = 100/numEars;
pL2 = pL1/numPos;
pL3 = pL2/nyquistIndx;
methodName = dsp.method{1};
switch lower(methodName)
    case 'fixedn'
        if isscalar(dsp.method{2})
            methodParamVec = dsp.method{2}*ones(nyquistIndx,1);
        else
            methodParamVec = dsp.method{2}(1:nyquistIndx);
        end
        convParam = 2; % Dummy value
    case 'maxn'
        methodParamVec = ones(nyquistIndx,1); % Dummy value
        convParam = 2; % Dummy value
    case 'exact'
        methodParamVec = dsp.method{2}*ones(nyquistIndx,1);
        convParam = 2; % Dummy value
    case 'formulan'
        p = dsp.method{2};
        if rho == inf
            alphaVal = 0.3442*p - 0.9584;
            betaVal = 1.758*exp(0.1033*p);
            gammaVal = 1.017-0.1041*p.^0.5086;
        else
            p1_a = 0.3442*p - 0.9584;
            p2_a = 3.369*p - 3.69;
            q1_a = -(0.7605*p + 1.704)./(p + 0.7461);
            
            p1_b = 1.758*exp(0.1033*p);
            p2_b = -(2.043*p - 0.1245);
            p3_b = 2.161*p - 1.964;
            q1_b = 0.001728*p.^2.977 - 3.088;
            q2_b = 2.679 - 0.0008429*p.^3.343;
            
            p1_c = 1.017-0.1041*p.^0.5086;
            p2_c = 0.1267*p - 2.859;
            p3_c = -(0.1227*p - 2.335);
            q1_c = -(0.02058*p + 2.968);
            q2_c = 0.06702*p + 2.169;
            
            alphaVal = (rho*p1_a+p2_a)/(rho+q1_a);
            betaVal = (rho^2*p1_b+rho*p2_b+p3_b)/(rho^2+rho*q1_b+q2_b);
            gammaVal = (rho^2*p1_c+rho*p2_c+p3_c)/(rho^2+rho*q1_c+q2_c);
        end
        methodParamVec = round(alphaVal+(betaVal*muVec.^gammaVal));
        convParam = 2; % Dummy value
        methodName = 'fixedn';
    case 'formulanmin'
        if rho == inf
            alphaVal = 1.4096;
            betaVal = 2.73;
            gammaVal = 0.8027;
        else
            alphaVal = ((1.4096*rho)+3.8985)/(rho-1.3552);
            betaVal = ((2.73*rho)-4.7457)/(rho-1.2067);
            gammaVal = ((0.8027*rho)-1.0104)/(rho-1.4451);
        end
        methodParamVec = round(alphaVal+(betaVal*muVec.^gammaVal));
        convParam = 2; % Dummy value
        methodName = 'fixedn';
    case 'formulanminsimp'
        if rho == inf
            alphaVal = 1.41;
            betaVal = 2.73;
            gammaVal = 0.8;
        else
            alphaVal = (1.41*rho+3.9)/(rho-1.36);
            betaVal = (2.73*rho-4.75)/(rho-1.21);
            gammaVal = (0.8*rho-1.01)/(rho-1.45);
        end
        methodParamVec = round(alphaVal+(betaVal*muVec.^gammaVal));
        convParam = 2; % Dummy value
        methodName = 'fixedn';
    case 'formulagd'
        methodParamVec = round(exp(1)*muVec);
        convParam = 2; % Dummy value
        methodName = 'fixedn';
    otherwise
        methodParamVec = dsp.method{2}*ones(nyquistIndx,1);
        if length(dsp.method) > 2
            convParam = dsp.method{3};
        else
            convParam = 2;
        end
end
cpVal = 0;
fprintf('Computing RS-HRTF..%d%%',cpVal);
for ii = 1:numEars
    H{ii,1} = zeros(nyquistIndx,numPos);
    N{ii,1} = zeros(nyquistIndx,numPos);
    TH{ii,1} = cell(nyquistIndx,numPos);
    for jj = 1:numPos
        for kk = 1:(nyquistIndx-1)
            % Progress
            ppVal = cpVal;
            cpVal = floor((((ii-1)*pL1 + (jj-1)*pL2) + (kk-1)*pL3)/10)*10;
            if cpVal ~= ppVal
                fprintf('..%d%%',cpVal);
            end
            [H{ii,1}(kk,jj),N{ii,1}(kk,jj),TH{ii,1}{kk,jj}] = ...
                computeSphereHRTF(head.a,source.r,theta{ii}(jj),...
                fVec(kk),{methodName,methodParamVec(kk),convParam},...
                dsp.normloc,FMAX);
        end
        kk = nyquistIndx;
        % Progress
        ppVal = cpVal;
        cpVal = floor((((ii-1)*pL1 + (jj-1)*pL2) + (kk-1)*pL3)/10)*10;
        if cpVal ~= ppVal
            fprintf('..%d%%',cpVal);
        end
        [H{ii,1}(kk,jj),N{ii,1}(kk,jj),TH{ii,1}{kk,jj}] = ...
            computeSphereHRTF(head.a,source.r,theta{ii}(jj),fVec(kk),...
            {methodName,methodParamVec(kk),convParam},dsp.normloc,FMAX);
        if ~retFlag
            % If irLen is even, take the absolute value of H at the Nyquist
            % frequency instead of dropping the imaginary part (since the
            % value of H at the Nyquist frequency must be real).
            if mod(dsp.irLen,2) == 0
                H{ii,1}(kk,jj) = abs(H{ii,1}(kk,jj));
            end
        end
    end
    
    if ~retFlag
        if dsp.causalflag
            if source.r == inf
                if ~isfield(dsp,'causalshift')
                    dsp.causalshift = (dsp.irLen/4)*(1000/dsp.fS);
                end
            else
                if ~isfield(dsp,'causalshift')
                    dsp.causalshift = source.r*1000/getSoundSpeed();
                end
            end
            shiftVal = dsp.causalshift*dsp.fS/1000;
            h{ii,1} = shiftSignal(ifft(H{ii,1},dsp.irLen,1,'symmetric'),...
                shiftVal);
        else
            h{ii,1} = ifft(H{ii,1},dsp.irLen,1,'symmetric');
        end
    end
end

switch numEars
    case 1
        if lFlag
            if ~retFlag
                head.hrirL = h{1,1};
            else
                head.hrtfL = H{1,1};
            end
            
            dsp.NMatL = N{1,1};
            dsp.THCellL = TH{1,1};
            source.thetaL = theta{1};
        elseif rFlag
            if ~retFlag
                head.hrirR = h{1,1};
            else
                head.hrtfR = H{1,1};
            end
            
            dsp.NMatR = N{1,1};
            dsp.THCellR = TH{1,1};
            source.thetaR = theta{1};
        else
            if ~retFlag
                head.hrirL = h{1,1};
            else
                head.hrtfL = H{1,1};
            end
            
            dsp.NMatL = N{1,1};
            dsp.THCellL = TH{1,1};
            source.thetaL = theta{1};
        end
    case 2
        if ~retFlag
            head.hrirL = h{1,1};
            head.hrirR = h{2,1};
        else
            head.hrtfL = H{1,1};
            head.hrtfR = H{2,1};
        end
        
        dsp.NMatL = N{1,1};
        dsp.NMatR = N{2,1};
        dsp.THCellL = TH{1,1};
        dsp.THCellR = TH{2,1};
        source.thetaL = theta{1};
        source.thetaR = theta{2};
end
fprintf('\nComputing RS-HRTF..done!\n');

end
