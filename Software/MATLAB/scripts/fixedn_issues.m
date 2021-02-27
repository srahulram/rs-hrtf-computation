%% Problem using fixed truncation order
%
% Illustration of problem using fixed, mu-independent order for computing
% the rigid-sphere HRTF.

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

%% Compute N_max and N_min values for mu = 0.01 and mu = 41, respectively.

[~,N_max,~] = computeSphereHRTF(0.09,inf,0,0.01*343/(2*pi*0.09),...
    {'maxn'});
[~,N_min,~] = computeSphereHRTF(0.09,inf,0,41*343/(2*pi*0.09),...
    {'exact',inf});
[H,~,~] = computeSphereHRTF(0.09,inf,0,0.01*343/(2*pi*0.09),...
    {'fixedn',N_min});

%% Comments
%
% Assuming all calculations are performed with double precision, we see 
% that N_min > N_max indicating that using a fixed N >= N_min cannot be 
% used if HRTFs that span mu from 0.01 to 41 (assuming a head radius of 
% 0.09 m) are to be computed. To see this, we attempt to compute H at mu = 
% 0.01 using the N_min value required to compute H "exactly" at mu = 41.
% We find that the resulting value of H is "not a number" (i.e., NaN).
