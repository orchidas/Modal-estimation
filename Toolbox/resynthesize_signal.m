function [irhat] = resynthesize_signal(mode_params, dur, fs, varargin)

	%% Synthesize signal in the time domain given modes amplitudes, 
    % frequencies and decay rates.
    % Inputs:
    % mode_params - struct containing mode frequencies, exponentiated decay
    % rates and complex amplitudes
	% dur - duration in samples of synthesized response
	% fs -  sampling rate
    % Optional
    % use_parallel_biquads - if true, the IR will be calculated from the
    %                        parallel combination of biquads, each synthesizing a mode
    % Outputs:
    % irhat - the synthesized time domain signal as a sum of modes.
    
    fmhat = mode_params.freqs;
    a1mhat = mode_params.decay_rate;
    gmhat = mode_params.amplitude;

    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'use_parallel_biquads', 0);
    parse(p,varargin{:});
    use_parallel_biquads = p.Results.use_parallel_biquads;

    if use_parallel_biquads
        num_modes = length(fmhat);
        sos = zeros(num_modes, 6);
        % set the zeros (numerator coefficients)
        if (size(gmhat,2) > 1)
            sos(:, 1) = gmhat(:, 2);
            sos(:, 2) = a1mhat .* (gmhat(:,1) .* sin(2*pi*fmhat/fs) - ...
                gmhat(:, 2) .* cos(2*pi*fmhat/fs)); 
        else
            sos(:, 1) = real(gmhat);
            sos(:, 2) = -a1mhat .* real(gmhat .* exp(-2*1j*pi*fmhat / fs));
        end
        % set the poles (denominator coefficients)
        sos(:, 4) = ones(num_modes, 1);
        sos(:, 5) = -2 * a1mhat .* cos(2*pi*fmhat / fs);
        sos(:, 6) = a1mhat .^ 2;
        % form an impulse and get its response
        t = (0:dur-1)/fs;
        impulse = [1; zeros(length(t)-1,1)];
        irhat = zeros(length(t), 1);
        % parallel biquad filtering
        for k = 1:num_modes
            irhat = irhat + sosfilt(sos(k,:), impulse); 
        end
    else
	    p = 1;
	    t = (0:dur-1)/fs;
        basis = exp((1j*2*pi*fmhat + log(a1mhat)*fs)*t).';

        if (size(gmhat,2) > 1)
            Vs = imag(basis);
            Vc = real(basis);
            V_full = [Vs(p+1:end,:), Vc(p+1:end,:)];
            irhat = real([zeros(p,1); V_full*gmhat(:)]);
        else
	        irhat = real([zeros(p,1); basis(p+1:end,:)*gmhat]);
        end
    end

end

