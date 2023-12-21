function [H] = parallel_biquad_response(g, fmhat, a1mhat, fs, wvec)

	%%
	% H(z) = \sum_{i=1}^M \frac{g_0,i + g_1,iz^{-1}}{1 - a_1,i z^{-1} - a_2,i
	% z^{-2}}
    % Inputs:
	% g - numerator coefficients
	% fmhat - mode frequencies
	% a1mhat - mode decay rates
	% fs - sampling rate
	% wvec - frequencies over which H is to be calculated
    % Outputs:
    % H - transfer function in the frequency domain calculated over wvec
	%%

	M = length(fmhat);
	p_abs = a1mhat.^2;
	p_ang = 2*pi*fmhat/fs;
	a1 = -2*p_abs.*cos(p_ang);
	a2 = p_abs.^2;
	g0 = g(1:M);
	g1 = g(M+1:end);
	sos = zeros(M,6);
	sos(:,1:3) = [g0, g1, zeros(M,1)];
	sos(:,4:6) = [ones(M,1), a1, a2];
	%some weird matlab issue with freqz
	if size(sos,1) > 1
	    H = freqz(sos, wvec);
	else
	    H = freqz(sos, length(wvec), wvec);
	end
end

