function [fmopt, a1mopt, niter, nfunc] = bandwise_td_optimize(ir, fmhat0, a1mhat0, f0, fs, deltaf, room_flag)
%% 
% Time domain mode optimization, but performed frequency bandwise for
% faster convergence
% Inputs
% ir - impulse response
% fmhat0 - estimated mode frequencies
% a1mhat0 - estimated mode dampings
% f0 - fundamental frequency (for centering filters)
% fs - sampling rate
% deltaf - how much frequency deviation is allowed
% room_flag - is input an RIR?
% optim_type - 'td' or 'fd' (time domain or pole optimization)
%
% Outputs
% fmopt, a1mopt - optimized mode frequencies and amplitures
% niter - number of iterations
% nfunc - number of function counts
% Author - Orchisama Das,2021
%%


maxFreq = max(fmhat0);
if room_flag
    nbands = 20; % band count, bands
    ft = (0:nbands)/nbands*fs/2;    % processing band edges, Hz
    wt = 2*pi*ft/fs;
    rho = -round((1.0674 * sqrt(2/pi * atan(0.06583*(fs/1000))) - 0.1916)*1000)/1000;
    wtw = atan2(((1-rho^2)*sin(wt)),((1+rho^2)*cos(wt) - 2*rho));
    ft = wtw*fs/(2*pi);    %warped band edges
    fh = (ft(1:end-1)+ft(2:end))/2; % processing band centers, Hz
    beta = 1.2*diff(ft/2); % subband filter passband bandwidth, Hz

    order = 6;  % subband filter order, poles
    ripplestop = 80;    % stopband ripple, dB
    ripplepass = 1.5;   % passband ripple, dB
    bLv = zeros(nbands, order+1);   %filter zero coeffs
    aLv = zeros(nbands, order+1);   %filter pole coeffs
    
    for b = 1:nbands   
        % design lowpass filter
        [bL, aL] = ellip(order, ripplepass, ripplestop, beta(b)*2/fs);
%        [bL,aL] = butter(order, beta(b)*2/fs);
        bLv(b,:) = bL;
        aLv(b,:) = aL;
    end
else   
    nbands = floor(maxFreq/f0);
    beta = f0/10;   % subband filter passband bandwidth, Hz
    fh = (1:nbands+1)*f0; %processing band centers
    ft = [0, fh + beta]; %processing band edges
    [bL, aL] = butter(4, beta*2/fs); %filter
end

ntaps = length(ir); % impulse response length, taps
t = (0:ntaps-1)'/fs;    % time axis, seconds
fmopt = zeros(length(fmhat0),1);
a1mopt = zeros(length(a1mhat0),1);
niter = 0;
nfunc = 0;

for b = 1:nbands
    
    % heterodyne, filter impulse response
    mu = exp(-1j*fh(b)*2*pi*t);
    irh = mu.*ir;
    
    % lowpass filter
    if room_flag
        irf = filter(bLv(b,:), aLv(b,:), irh);
    else
        irf = filter(bL, aL, irh);
    end
    
    %undo heterodyne
    irf = real(conj(mu).*irf);
    
    %find corresponding mode frequencies in band
    indexb = find((abs(fmhat0) > ft(b)) & (abs(fmhat0) < ft(b+1)));
    fmbhat0 = fmhat0(indexb);
    a1mbhat0 = a1mhat0(indexb);
    
    %optimize now
    if ~isempty(indexb)
        [fmbhat,a1mbhat,iter,func] = optimize_modes_td(irf,fmbhat0,a1mbhat0,[],fs,deltaf);
        fmopt(indexb) =  fmbhat;
        a1mopt(indexb) =  a1mbhat;
        niter = niter + iter;
        nfunc = nfunc + func;
    end
    
    disp(['Frequency band ', num2str(b), ' has been processed']);

end
    
    
    
    
    
end

