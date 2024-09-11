function [mode_params, irhat] = frequency_warped_modal(ir,fs, f0, opt_flag, room_flag, rhol)
%%
% Frequency warped modal estimation
% Inputs - 
% ir - measured signal
% fs - sampling rate
% f0 - fundamental frequency (if known)
% opt_flag - optimization flag - 0 or 1
% room_flag - is input an RIR? use whenever f0 is not known
% rhol - warping coefficient for low frequencies
% Outputs
% mode params - struct containing mode frequencies, decay rates and amps
% irhat - reconstructed signal in the time domain
% Author - Orchisama Das, 2020
%%

ntaps = length(ir);
t = (0:ntaps-1)'/fs;
dur = min(length(ir),2.0*fs); %data used to calculate mode amplitudes

% modeling parameters
nh = 2048;  % Hankel matrix dimension, columns
p = 1;  % offset Hankel matrix shift, taps
kappa1 = -40;   % minimum mode singular value,
kappa2 = -40;
nbands = 60;

%default warping factor (if not specified is Bark warping)
if nargin < 6
    rhol = -round((1.0674 * sqrt(2/pi * atan(0.06583*(fs/1000))) - 0.1916)*1000)/1000;
end
fc = acos(-rhol)/pi * (fs/2);


%% estimate high frequency mode parameters 


%introduce some extra damping
alpha = 2;
irh = ir .* exp(-alpha*t); %damp the modes
[fmhat0, a1mhat0, ~, ~] = hvmodel_freqs_decay(irh, nh, p ,fs, kappa1);
a1mhat0 = a1mhat0.^exp(alpha*(ntaps-1)/fs);   %undo damping

   
%% estimate low frequency mode parameters 


% warp impulse response
irl_warp = warpfir(ir,rhol);
%a1mhat is the pole radius, not the decay rate = exp(decay rate)
[fmhat,a1mhat, S, nm] = hvmodel_freqs_decay(irl_warp,nh,p,fs,kappa2);
%unwarp warped mode frequencies and decay rates
[fmhat_uw, a1mhat_uw, ~, ~] = unwarp_modes(fmhat, a1mhat, rhol, fs);


%% combine two sets of modes to include high frequency modes

ind1 = find(abs(fmhat0) >= fc);
ind2 = find(abs(fmhat_uw) < fc);
fmhat_com = [fmhat_uw(ind2);fmhat0(ind1)];
a1mhat_com = [a1mhat_uw(ind2); a1mhat0(ind1)];
[fmhat_com,ord] = sort(fmhat_com);
a1mhat_com = a1mhat_com(ord);

if ~room_flag
    lim = find(abs(fmhat_com) < nbands*f0);
    fmhat_com = fmhat_com(lim);
    a1mhat_com = a1mhat_com(lim);
end
% make sure system is stable
a1mhat_com(a1mhat_com >= 1.00) = 0.9999;


%optimize (or not)
if opt_flag
    %a1mhat_opt is the decay rate = log(pole radius)
    deltaf = 2;   % maximum frequency deviation allowed
    [fmhat_opt, a1mhat_opt] = bandwise_td_optimize(ir(1:dur), fmhat_com, ...
        abs(log(a1mhat_com)), f0, fs, deltaf, room_flag);

    [gmhat, irhat] = estimate_mode_amps(ir, fmhat_opt, a1mhat_opt, p, t, fs, dur, opt_flag);
    mode_params.freqs = fmhat_opt;
    mode_params.decay_rate = exp(-a1mhat_opt);
    % make sure system is stable
    mode_params.decay_rate(mode_params.decay_rate >= 1.00) = 0.9999;
    mode_params.amplitude = [gmhat(1:length(fmhat_opt)), gmhat(length(fmhat_opt)+1:end)];
else
    [gmhat, irhat] = estimate_mode_amps(ir, fmhat_com, a1mhat_com, p, t, fs, dur, opt_flag);
    mode_params.freqs = fmhat_com;
    mode_params.decay_rate = a1mhat_com;
    mode_params.amplitude = gmhat;
end



end

