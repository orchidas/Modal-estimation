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
% fmopt, a1mopt - optimized mode frequencies and decay rates
% niter - number of iterations
% nfunc - number of function counts
% Author - Orchisama Das,2021
%%


maxFreq = max(fmhat0);
if room_flag
    nbands = 20; % band count, bands
    order = 6;  % subband filter order, poles
    ripplestop = 80;    % stopband ripple, dB
    ripplepass = 1.5;   % passband ripple, dB
    [bLv, aLv, fh, ft, ~] = filter_rir_in_subbands(fs, nbands, order, ripplestop, ripplepass);
    
else   
    nbands = floor(maxFreq/f0);
    order = 4;
    [bL, aL, fh, ft] = filter_harmonic_signal_in_subbands(fs, nbands, order, f0);
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

% return exp(-alpha)
a1mopt = exp(-a1mopt);
% make sure system is stable
a1mopt(a1mopt >= 1.00) = 0.9999;
    
end

