function [irw, duration] = warpfir(argin1, argin2, argin3)
% WARPFIR - FIR filter frequency warping.
%
% IRW = warpfir(IR, RHO) warps input impulse responses stored columnwise in IR
% according to the first-order conformal map zeta = (RHO + z)/(1 + RHO*z) to
% produce warped impulse responses IRW.  When RHO has two elements, [Fz Fzeta],
% a conformal map parameter is determined such that the frequency Fz is mapped
% into the frequency Fzeta, both measured in radians/pi.
%
% The call warpfir(IR, RHO, LEVELDB) windows the warped impulse responses to a
% a length given by the greater of the length of the input impulse response
% and that part of any impulse response containing 1 - 10^(LEVELDB/20) of the
% impulse response energy.
%
% See also WARPIIR, WARPFREQ, BCM2FS, FS2BCM, PRONYW. 
% (c) Copyright 2001 Abel Innovations.  All rights reserved.  CONFIDENTIAL.
%
% Jonathan S. Abel
% Created: 26-Apr-2001.
% Revised: 26-Apr-2001, v1.0.
% Version: v1.0
 
 
%% process input
%%
if (nargin <= 1),
     %% not enough input arguments
     disp(['warpfir: ERROR - at least two input arguments required.']);
     return;
 
elseif (nargin == 2),
     ir = argin1;
     if (prod(size(ir)) == length(ir)),
         nsamp = length(ir);
         ir = ir(:);
     else
         nsamp = size(ir,1);
     end; 
    %% warping specification
    temp = argin2; 
 elseif (nargin == 3),
     %% impulse response
     ir = argin1;
     if (prod(size(ir)) == length(ir)),
         nsamp = length(ir);
         ir = ir(:);
     else
         nsamp = size(ir,1);
     end;
 
     %% warping specification
     temp = argin2;
 
     %% impulse response windowing
     level = 10^(argin3/20); 
end;
 
if (prod(size(temp)) == 1),
     %% conformal map parameter specified
     rho = temp; 
else
     %% linear, warped frequencies specified
     rho = sin((temp(1) - temp(2)) * pi/2) / sin((temp(1) + temp(2)) * pi/2); 
end;
 
 
%% initialization
%%
nbinsmax = 65536;        %% maximum fft half bandwidth, bins 

%% warp transfer function
%% 
%% form linear and warped frequency axes
stretch = (1 + abs(rho))/(1 - abs(rho));
nbins = min(nbinsmax, 2^nextpow2(nsamp * stretch));
w = pi*[0:nbins]'/nbins;
z = exp(sqrt(-1) * w);
zeta = (z - rho) ./ (1 - rho*z);
ww = angle(zeta);
 
%% form, interpolate transfer function
temp = fft(ir, 2*nbins);
tf = temp(1:nbins+1,:);
 
temp = interp1(w, tf, ww, 'cubic');
tfw = [temp; conj(flipud(temp(2:nbins,:)))]; 

%% compute warped impulse response
temp = real(ifft(tfw, 2*nbins));
irw = temp(1:nbins,:);
%Orchi-addition
% irw = irw(1:nsamp);

%% window impulse response
if (nargin > 2);
     duration = sum(flipud(cumsum(flipud(irw.^2))) > level^2 * ones(size(irw,1),1)*sum(irw.^2));
     irw = irw(1:max([duration, length(ir)]),:);
end;
 
%% transpose impulse response
if (size(argin1,1) == 1),
     irw = irw';
end;
