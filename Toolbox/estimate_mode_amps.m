function [gmhat, irhat] = estimate_mode_amps(ir, fmhat, a1mhat, p, t, fs, dur, opt_flag)

%% 
% Estimate mode amplitudes using linear least squares given the frequencies
% and decay rates
% Inputs:
% ir - time domain signal to be modeled
% fmhat - estimated mode frequencies
% a1mhat - estimated mode decay rates (exponentiated : exp(-alpha_m)
% p - index from where to start the least squares estimation
% t - time vector
% opt_flag - have the modes been optimised?
% dur - duration of the signal used for estimation
% Outputs:
% gmhat - the estimated complex mode amplitudes
% irhat - the reconstructed modal signal
%%

if nargin < 7
    dur = length(t);
    opt_flag = 0;
end

%estimate decay rate with least squares

if opt_flag
    basis = exp(t*(1i*2*pi*fmhat - a1mhat*fs).');
    Vs = imag(basis);
    Vc = real(basis);
    V = [Vs(p+1:dur,:), Vc(p+1:dur,:)];
    gmhat = V \ ir(p+1:dur);
    V_full = [Vs(p+1:end,:), Vc(p+1:end,:)];
    irhat = real([zeros(p,1); V_full*gmhat]);

else
    basis = exp(1j*2*pi*t*fmhat' + t*log(a1mhat)'*fs);
    gmhat = basis(p+1:dur,:) \ ir(p+1:dur);
    irhat = real([zeros(p,1); basis(p+1:end,:)*gmhat]);
end

