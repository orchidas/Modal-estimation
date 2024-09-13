function [mode_params, irhat] = frequency_zoomed_ARMA(ir, fs, r, opt_flag, room_flag, varargin)

%%
% Estimate modal parameters with frequency-zoomed ARMA. 
% Inputs:
% ir - time-domain signal whose modes are to be estimated
% f - sampling frequency
% r - downsampling factor
% opt_flag - whether to optimise the estimated modes
% room_flag - is the signal a room impulse response
% Optional:
% f0 - fundamental frequency of the signal (provided for harmonic signals)
% Outputs
% mode_params - struct containing the mode frequencies, decay rates and
% amplitudes
%%

p = inputParser;
p.KeepUnmatched = true;
addParameter(p, 'f0', []);
addParameter(p, 'plot_flag', 0)
parse(p,varargin{:});
f0 = p.Results.f0;
plot_flag = p.Results.plot_flag;

ntaps = length(ir); % impulse response length, taps
t = (0:ntaps-1)'/fs;    % time axis, seconds
dur = min(length(ir),2.0*fs); %data used to calculate mode amplitudes
p = 1;
fsr = fs/r; % downsampled impulse response sampling rate, Hz
delf = 2;

if room_flag
    
    % band processing variables
    nbands = 40; % band count, bands
    Nmax = 30;
    Pmax = 150;    %maximum ARMA model order
    order = 7;  % subband filter order, poles
    ripplestop = 80;    % stopband ripple, dB
    ripplepass = 1.5;   % passband ripple, dB

    [bLv, aLv, fh, ft, beta] = filter_rir_in_subbands(fs, nbands, order, ripplestop, ripplepass);

    %model order should depend on filter bandwidth
    betan = beta./max(beta);
    N = max(10*ones(1,length(beta)),round(tanh(3*betan)/tanh(3) * Nmax));
    P = max(20*ones(1,length(beta)),round(tanh(3*betan)/tanh(3) * Pmax));  

else
    
    nbands = 60;
    B = 10^-4;
    order = 4;
    [bL, aL, fh, ft] = filter_harmonic_signal_in_subbands(fs, nbands, order, f0, B);

    N = 12; %model zero order
    P = 12; %mode pole order
end


%% loop through bands
fmhat = [];
a1mhat = [];
for b = 1:nbands
    % % form downsampled band response

    % heterodyne, filter impulse response
    mu = exp(-1j*fh(b)*2*pi*t);
    irh = mu.*ir;
    
    % lowpass filter
    if room_flag
        irf = filter(bLv(b,:), aLv(b,:), irh);
    else        
        irf = filter(bL, aL, irh);
    end
    %decimate
    irb = irf(1:r:end);
    

    % ARMA parameter estimation with Steiglitz McBride
    if room_flag
        [bb,aa] = stmcb(irb, N(b), P(b));
    else
        [bb,aa] = stmcb(irb, N, P);   
    end
    
    if ~(any(isnan(bb)) || any(isnan(aa)))
        [z,poles,K] = tf2zp(bb,aa);
        % poles_uw = (poles.^(1/r)) * exp(1j*2*pi*fh(b)/fs);
        [fmbhat,order] = sort(angle(poles)* fsr/(2*pi));
        a1mbhat = abs(poles(order));

        % assign variables
        indexb = find((fmbhat+fh(b) > ft(b)) & (fmbhat+fh(b) < ft(b+1)));
        fmbchat = fmbhat(indexb)+fh(b); %undo heterodyne - move frequencies up the spectrum
        a1mbchat = a1mbhat(indexb).^(1/r);  %undo damping

        if opt_flag
            % undo heterodyne
            irf = real(conj(mu).*irf);
            [fmbchat, a1mbchat] = optimize_modes(irf, fmbchat, abs(log(a1mbchat)),[], fs, delf);
        end

        fmhat = [fmhat; fmbchat];   %concatenate with previous band estimates
        a1mhat = [a1mhat; a1mbchat];  
        disp(['Frequency band ', num2str(b), ' has been processed']);


        %% plot results
        if (plot_flag)
            rtaps = length(irb);
            tr = (0:rtaps-1)'/fsr;
            index = find(abs(a1mbhat) < 1);
            basis = exp(1j*2*pi*tr*fmbhat(index)' + tr*log(a1mbhat(index)')*fsr);
            gmbhat = basis(p+1:end,:) \ irb(p+1:end);
            irbhat = [zeros(p,1); basis(p+1:end,:)*gmbhat];
    
    
            figure(2); clf;
            set(2, 'Position', [716 34 560 917]);
    
            % plot given, estimated responses
            scale = max(abs(irb));
            plot((p:rtaps-1)/fsr, real([irb(p+1:end) irbhat(p+1:end)]) * [eye(2) [-1 1]']/scale, '-', ...
                (p:rtaps-1)/fsr, imag([irb(p+1:end) irbhat(p+1:end)]) * [eye(2) [-1 1]']/scale + 1); grid;
            title('given, estimated band responses, error');
            xlabel('time, seconds'); ylabel('amplitude');
            xlim([0 1]); ylim([-0.5 1.5]);
        end
    end
    
end

% deal with unstable poles
a1mhat(a1mhat >= 1) = 0.9999;
[gmhat,irhat] = estimate_mode_amps(ir(1:dur), fmhat, a1mhat, p, t, fs, dur, opt_flag);
if opt_flag
    mode_params.freqs = fmhat;
    mode_params.decay_rate = a1mhat;
    mode_params.amplitude = [gmhat(1:length(fmhat)), gmhat(length(fmhat)+1:end)];
else
    mode_params.freqs = fmhat;
    mode_params.decay_rate = a1mhat;
    mode_params.amplitude = gmhat;
end



end