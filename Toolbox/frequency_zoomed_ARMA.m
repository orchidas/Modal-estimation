function [mode_params, irhat] = frequency_zoomed_ARMA(ir, fs, f0, r, opt_flag, room_flag)


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
    ft = (0:nbands)/nbands*fs/2;    % processing band edges, Hz
    wt = 2*pi*ft/fs;
    rho = -round((1.0674 * sqrt(2/pi * atan(0.06583*(fs/1000))) - 0.1916)*1000)/1000;
    wtw = atan2(((1-rho^2)*sin(wt)),((1+rho^2)*cos(wt) - 2*rho));
    ft = wtw*fs/(2*pi);    %warped band edges
    fh = (ft(1:end-1)+ft(2:end))/2; % processing band centers, Hz
    beta = 1.2*diff(ft/2); % subband filter passband bandwidth, Hz
    
    %model order should depend on filter bandwidth
    betan = beta./max(beta);
    N = max(10*ones(1,length(beta)),round(tanh(3*betan)/tanh(3) * Nmax));
    P = max(20*ones(1,length(beta)),round(tanh(3*betan)/tanh(3) * Pmax));  

    order = 7;  % subband filter order, poles
    ripplestop = 80;    % stopband ripple, dB
    ripplepass = 1.5;   % passband ripple, dB
    bLv = zeros(nbands, order+1);   %filter zero coeffs
    aLv = zeros(nbands, order+1);   %filter pole coeffs
    
    for b = 1:nbands   
        % design lowpass filter
        [bL, aL] = ellip(order, ripplepass, ripplestop, beta(b)*2/fs);
        bLv(b,:) = bL;
        aLv(b,:) = aL;
    end
else
    
    % fUP = 16000;    %piano partials won't exceed 16k
    % nbands = floor(fUP/(2*f0));
    nbands = 60;
    B = 10^-4;
    n = 1:nbands+1;
    fh = f0*n.*sqrt(1+B*n.^2); %processing band centers
    beta = f0/10;
    ft = [0, fh + beta];
    N = 12; %model zero order
    P = 12; %mode pole order
    [bL, aL] = butter(4, beta*2/fs);
end



%% estimate band mode parameters

% % %plot subbands
% N = 1000;
% w = linspace(-pi,pi,N+1);
% fHz = (w/pi)*(fs/2);
% H = freqz(bL,aL,w);
% Hcon = [];
% for b = 1:nbands
%     [val,sh] = min(abs(ft(b) - fHz));
%     Hsh = circshift(H,sh);
%     Hcon = [Hcon; Hsh(1:N/2+1)];
% end
% figure(1);
% plot(fHz(N/2+1:end), 20*log10(abs(Hcon)));grid on;
% xlim([0,fUP]);
% ylim([-80,5]);
% xlabel('Frequency(Hz)');
% ylabel('Magnitude(dB)');


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

% deal with unstable poles
a1mhat(a1mhat >= 1) = 0.9999;
[gmhat,irhat] = estimate_mode_amps(ir(1:dur), fmhat, a1mhat, p, t, fs, dur, opt_flag);
if opt_flag
    mode_params = [fmhat, a1mhat, gmhat(1:length(fmhat)), gmhat(length(fmhat)+1:end)];
else
    mode_params = [fmhat, a1mhat, gmhat];
end



end