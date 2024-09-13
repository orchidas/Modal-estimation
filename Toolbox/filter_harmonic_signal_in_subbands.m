function [bL, aL, fh, ft] = filter_harmonic_signal_in_subbands(fs, nbands, order, f0, varargin)
%Filter a harmonic signal in subbands
%   Inputs:
%   fs - sampling rate
%   nbands - number of subbands
%   f0 - fundamental frequency of signal in Hz
%   B (optional) - inharmonicity factor
%   Outputs:
%   b, a - filter coefficients
%   ft - band edges in Hz
%   fh - band centres in Hz

if nargin > 4
    B = varargin{1};
else
    B = 0;
end

n = 1:nbands+1;
fh = f0*n.*sqrt(1+B*n.^2); %processing band centers
beta = f0/10;
ft = [0, fh + beta];
[bL, aL] = butter(order, beta*2/fs);

