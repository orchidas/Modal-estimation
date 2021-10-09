function [gmhat, irhat] = estimate_mode_amps(ir, fmhat, a1mhat, p, t, fs, dur, opt_flag)

if nargin < 7
    dur = length(t);
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

