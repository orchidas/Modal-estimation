function [irhat] = resynthesize_signal(fmhat, a1mhat, gmhat, dur, fs)

	%synthesize signal given modes amplitudes, frequencies and decay rates
	% dur is duration in seconds of synthesized response
	% fs is sampling rate
	p = 1;
	t = (0:dur-1)/fs;
	basis = exp(1j*2*pi*t*fmhat' + t*log(a1mhat)'*fs);
	irhat = real([zeros(p,1); basis(p+1:end,:)*gmhat]);

end

