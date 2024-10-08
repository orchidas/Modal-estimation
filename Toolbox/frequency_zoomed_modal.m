function [mode_params, irhat] = frequency_zoomed_modal(ir, fs, f0, r, opt_flag, varargin)

%%
% Frequency-zoomed (subband) modal estimation
% Inputs : 
% ir - measured signal
% fs - sampling rate
% f0 - fundamental frequency (if any)
% r - downsampling factor
% opt_flag - optimization flag, 0 or 1
% Optional:
% room_flag - RIR flag, 0 or 1
% plot - whether to plot intermediate fits
% optim_type - 'td' or 'fd' (time domain or pole optimization)
% Outputs :
% mode params - struct containing mode frequencies, decay rates and amps
% irhat - reconstructed signal using estimated modes
% Author - Orchisama Das, 2020
%%



switch(nargin)
  
    case 6
        room_flag = varargin{1};
        plot_flag = 0;
        opt_type = 'td';
    case 7
        room_flag = varargin{1};
        plot_flag = varargin{2};
        opt_type = 'td';
    case 8
        room_flag = varargin{1};
        plot_flag = varargin{2};
        opt_type = varargin{3};
    otherwise
        room_flag = 0;
        plot_flag = 0;
        opt_type = 'td';

end


ntaps = length(ir); % impulse response length, taps
t = (0:ntaps-1)'/fs;    % time axis, seconds
dur = min(length(ir),2.0*fs); %data used to calculate mode amplitudes

% Hankel-Vandermonde modeling parameters
nhmax = 2048;  % Hankel matrix dimension, columns
p = 1;  % offset Hankel matrix shift, taps
kappa = -40;    % minimum mode singular value, dB
fsr = fs/r; % downsampled impulse response sampling rate, Hz
delf = 2; %frequency deviation allowed in optimization


%% RIR modal estimation
if room_flag
    
    % band processing variables
    nbands = 20; % band count, bands
    nmbmax = 100*ones(nbands,1);   % maximum band mode count, modes
    order = 7;  % subband filter order, poles
    ripplestop = 80;    % stopband ripple, dB
    ripplepass = 1.5;   % passband ripple, dB
    [bLv, aLv, fh, ft, ~] = filter_rir_in_subbands(fs, nbands, order, ripplestop, ripplepass);
    
%% Piano modal estimation
else     
    nbands = 60;
    B = 10^-4;
    order = 4;
    [bL, aL, fh, ft] = filter_harmonic_signal_in_subbands(fs, nbands, order, f0, B);
    nmbmax = 12;   % maximum band mode count, modes
end

%% estimate band mode parameters

% loop through bands
fmhat = [];
a1mhat = [];

for b = (1:nbands)
    % % form downsampled band response

    % heterodyne, filter impulse response
    mu = exp(-1j*fh(b)*2*pi*t);
    irh = mu.*ir;
    
    %low pass filter
    if room_flag
        irf = filter(bLv(b,:), aLv(b,:), irh);
    else
        irf = filter(bL, aL, irh);
    end
    % decimate
    irb = r*irf(1:r:end);
    


    %% estimate frequencies, dampings

    nh = min(nhmax, length(irb)/2);
    [fmbhat, a1mbhat, S, nmb] = hvmodel_freqs_decay(irb, nh, p, fsr, kappa);
    
    if room_flag
        % if mode budget is not completely alloted, redistribute it in remaining bands
        if nmb <= nmbmax(b) 
            ndiff = nmbmax(b) - nmb;
            nmbmax(b+1:end) = nmbmax(b+1:end) + round(ndiff/(nbands-b-1));
        % or compress number of modes 
        else
            nmb = nmbmax(b);
        end
    end
    
    % assign variables
    indexb = find((abs(fmbhat)+fh(b) > ft(b)) & (abs(fmbhat)+fh(b) < ft(b+1)));
    % undo heterodyne - move frequencies up the spectrum
    fmbchat = fmbhat(indexb)+fh(b); 
    a1mbchat = a1mbhat(indexb).^(1/r);
    
    % optimize modes (if flag is true)
    if opt_flag
        % undo heterodyne
        irf = real(conj(mu).*irf);
        if strcmp(opt_type, 'td')
            [fmbchat, a1mbchat] = optimize_modes_td(irf, fmbchat, abs(log(a1mbchat)),[], fs, delf);
        else
            [fmbchat, a1mbchat] = constrained_pole_optimization(irf, fsr, fmbchat, a1mbchat, ft(b), ft(b+1));
        end
    end
        
    fmhat = [fmhat; fmbchat];   %concatenate with previous band estimates
    a1mhat = [a1mhat; a1mbchat];    
    disp(['Frequency band ', num2str(b), ' has been processed']);

    %% plot results
    

    if plot_flag

        rtaps = length(irb);
        tr = (0:rtaps-1)'/fsr;
        index = find(abs(a1mbhat) < 1);
        basis = exp(1j*2*pi*tr*fmbhat(index)' + tr*log(a1mbhat(index)')*fsr);
        gmbhat = basis(p+1:end,:) \ irb(p+1:end);
        irbhat = [zeros(p,1); basis(p+1:end,:)*gmbhat];

        figure(2); clf;
        set(2, 'Position', [716 34 560 917]);

        % plot Hankel, Vandermonde eigenstructure
        subplot(4,1,1);
        plot(1:nh, 20*log10(diag(S)./S(1,1)), '.');hold on; 
        plot(1:nmb, 20*log10(diag(S(1:nmb,1:nmb))./S(1,1)), '.'); hold off; grid;
        title(['Hankel rir matrix singular values, ', int2str(nmb), ' modes, band ', int2str(b)]);
        xlabel('mode index'); ylabel('amplitude, dB');
        ylim([-60 0]);

        % plot mode frequencies, dampings
        rtlims = [2e1 5e4];
        subplot(4,1,2);
        semilogy((fmbhat+fh(b))/1000, log(0.001)./log(a1mbhat)*1000/fsr, '.', ...
            [1; 1]*ft(b:b+1)/1000, rtlims'*[1 1], '-'); grid;
        title(['estimated mode decay times, ', int2str(nmb), ' modes']);
        xlabel('frequency, kHz'); ylabel('60 dB decay time, ms');
        xlim([0 fs/2000]); ylim(rtlims);

        % plot given, estimated responses
        scale = max(abs(irb));
        subplot(4,1,[3 4]);
        plot((p:rtaps-1)/fsr, real([irb(p+1:end) irbhat(p+1:end)]) * [eye(2) [-1 1]']/scale, '-', ...
            (p:rtaps-1)/fsr, imag([irb(p+1:end) irbhat(p+1:end)]) * [eye(2) [-1 1]']/scale + 1); grid;
        title('given, estimated band responses, error');
        xlabel('time, seconds'); ylabel('amplitude');
        xlim([0 1]); ylim([-0.5 1.5]);
        drawnow;
    end


end

%% 

a1mhat(a1mhat >= 1) = 0.9999;
[gmhat,irhat] = estimate_mode_amps(ir(1:dur), fmhat, a1mhat, p, t, fs, dur, opt_flag);

if opt_flag && strcmp(opt_type, 'td')
    mode_params.freqs = fmhat;
    mode_params.decay_rate = a1mhat;
    mode_params.amplitude = [gmhat(1:length(fmhat)), gmhat(length(fmhat)+1:end)];
else
    mode_params.freqs = fmhat;
    mode_params.decay_rate = a1mhat;
    mode_params.amplitude = gmhat;
end



end

