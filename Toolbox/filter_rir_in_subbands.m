function [bLv, aLv, fh, ft, beta] = filter_rir_in_subbands(fs, nbands, order, ripplestop, ripplepass, varagin)
% Take a broadband RIR in the time domain and filter it into subbands
%   Inputs:
%   fs - sampling rate in Hz
%   nbands - number of subbands
%   order - order of the filter
%   ripplestop - stopband ripple in dB
%   ripplepass - passband ripple in dB
%   plot_flag - whether to plot the filter response
%   Outputs:
%   bLv - matrix of numerator coefficients of size nbands x order+1
%   aLv - matrix of denominator coefficients of size  nbands x order+1
%   fh - band centre frequencies in Hz
%   ft - band edges in Hz
%   beta - subband filter bandwidth in Hz

if (nargin > 5)
    plot_flag = varagin{1};
else
    plot_flag = 0;
end

ft = (0:nbands)/nbands*fs/2;    % processing band edges, Hz
wt = 2*pi*ft/fs;
rho = -round((1.0674 * sqrt(2/pi * atan(0.06583*(fs/1000))) - 0.1916)*1000)/1000;
wtw = atan2(((1-rho^2)*sin(wt)),((1+rho^2)*cos(wt) - 2*rho));
ft = wtw*fs/(2*pi);    %warped band edges
fh = (ft(1:end-1)+ft(2:end))/2; % processing band centers, Hz
beta = 1.2*diff(ft/2); % subband filter passband bandwidth, Hz

%model order should depend on filter bandwidth
bLv = zeros(nbands, order+1);   %filter zero coeffs
aLv = zeros(nbands, order+1);   %filter pole coeffs

for b = 1:nbands   
    % design lowpass filter
    [bL, aL] = ellip(order, ripplepass, ripplestop, beta(b)*2/fs);
    bLv(b,:) = bL;
    aLv(b,:) = aL;
end

% plot subbands
if (plot_flag)
    N = 1000;
    w = linspace(-pi,pi,N+1);
    fHz = (w/pi)*(fs/2);
    H = freqz(bL,aL,w);
    Hcon = [];
    for b = 1:nbands
        [val,sh] = min(abs(ft(b) - fHz));
        Hsh = circshift(H,sh);
        Hcon = [Hcon; Hsh(1:N/2+1)];
    end
    figure(1);
    plot(fHz(N/2+1:end), 20*log10(abs(Hcon)));grid on;
    xlim([0,fUP]);
    ylim([-80,5]);
    xlabel('Frequency(Hz)');
    ylabel('Magnitude(dB)');
end
