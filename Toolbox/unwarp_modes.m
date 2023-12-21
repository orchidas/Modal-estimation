function [fmhat_uw, a1mhat_uw, modes_w, modes_uw] = unwarp_modes(fmhat, a1mhat, rho, fs)

%% 
% Unwarp warped modes
% Inputs:
% fmhat -  warped frequencies in Hz
% a1mhat - warped decay rates
% rho - warping factor
% fs -sampling frequency
% Outputs:
% fmhat_uw - unwarped mode frequencies
% a1mhat_uw - unwarped mode decay rates
% modes_w - warped modes expressed as a complex number in the unit circle
% modes_uw - unwarped modes expressed as a complex number in the unit circle
%%

modes_w = exp(-1j*2*pi*(fmhat/fs)).* a1mhat;
lambdam = log((-rho + modes_w)./(1 - rho*modes_w)); 
fmhat_uw = imag(lambdam)/(2*pi)*fs;
a1mhat_uw = exp(real(lambdam));

%sort frequencies and amplitudes
[fmhat_uw, order] = sort(fmhat_uw);
a1mhat_uw = a1mhat_uw(order);
modes_uw = exp(1j*2*pi*(fmhat_uw/fs)).* -a1mhat_uw;

end

