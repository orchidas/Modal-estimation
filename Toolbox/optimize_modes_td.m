function [fmopt, a1mopt, niter, nfunc] = optimize_modes_td(h, fmhat0, a1mhat0, t60_lims, fs, deltaf)
%% 
% Optimize mode freqs and dampings to fit impulse response
% Inputs
% h - Impulse response
% fmhat0, a1mhat0 - initial set of mode frequencies, decay rates
% t60_lims - lower and upper bound T60 in seconds (if any)
% fs - sampling rate
% deltaf - maximum frequency deviation allowed
% Outputs
% fmopt, a1mopt - optimized mode frequencies, decay rates
% niter, nfunc - Number of optimization iterations/function evaluations
% Notes:
% We are minimizing for decayrate * fs
% make sure a1mhat0 is the decay rate (not pole radius) and is positive
% Author - Orchisama Das, 2020
%%

T = length(h);
M = length(fmhat0);

function[err, grad] = costfn(theta)

    t = (0:T-1)'./fs;
    V = exp(t*(1i*2*pi*theta(M+1:end) - theta(1:M)*fs).'); %OUTER PRODUCT
    Vs = imag(V);
    Vc = real(V);
    V_p = [Vs, Vc];
    gm = V_p\h;
    hhat = V_p*gm;
    err = 0.5*norm(h - hhat)^2;


    figure(1); 
    plot((0:T-1)*1000/fs, [h, hhat]); grid;
    xlabel ('time, milliseconds'); 
    ylabel('amplitude');
    ylim([-max(abs(h)) max(abs(h))]);
    drawnow;


    %calculate gradient for faster computation
    if nargout > 1

        a1m = theta(1:M);
        wm = 2*pi*theta(M+1:end);
        gs = gm(1:M);
        gc = gm(M+1:end);

        gsmat = repmat(gs.',[T,1]);
        gcmat = repmat(gc.', [T,1]);
        tmat = repmat(t,[1,M]);
        D = zeros(T, 2*M);
        D(:,1:M) = -fs*tmat.*exp(-fs*t*a1m').*(sin(t*wm').*gsmat + cos(t*wm').*gcmat);
        D(:,M+1:end) = 2*pi*tmat.*exp(-fs*t*a1m').*(cos(t*wm').*gsmat - sin(t*wm').*gcmat);
        grad = D.'*(hhat - h);
 

    end

end


%optimization stopping criteria
function stop = outfun(theta,optimValues,state) 
    stop = false;
    if optimValues.fval < 1e-3 || optimValues.funccount > 1000 
        stop = true; 
        disp('Stopping, small enough error')
    end
end


%limits on decay rate and frequency
del = deltaf;   %max frequency deviation from initial value
alpha_del = 0.5*a1mhat0;  %max decay rate deviation from initial value

if isempty(t60_lims)
    lb = [a1mhat0 - alpha_del; fmhat0 - del];        
    ub = [a1mhat0 + alpha_del; fmhat0 + del];
else
    alpha_lims = t60_to_decay(t60_lims,fs);
    lb_alpha = [alpha_lims(2)*ones(M,1), a1mhat0 - alpha_del];
    ub_alpha = [alpha_lims(1)*ones(M,1), a1mhat0 + alpha_del];
    lb = [max(lb_alpha,[],2); fmhat0 - del];        
    ub = [min(ub_alpha,[],2); fmhat0 + del];
  
end

%initial point
theta0 = [a1mhat0; fmhat0];
%inequality constraints (to make sure frequencies are ordered)
e = ones(M,1);
A = [zeros(M,M), full(spdiags([e -e], 0:1, M,M))];
b = zeros(M,1);


figure(2);
subplot(211); plot(decay_to_t60(a1mhat0,fs), 'o'); grid;hold on;
subplot(212); plot(abs(fmhat0), 'o'); grid;hold on;
    
options = optimoptions('fmincon','Display','iter','SpecifyObjectiveGradient'...
    ,true,'OutputFcn',@outfun, 'CheckGradient',false, 'MaxIterations', 100);
[theta, fval, exitflag, output] = fmincon(@costfn, theta0, A,b,[],[],lb,ub,[],options);
niter = output.iterations;
nfunc = output.funcCount;
a1mopt = theta(1:M);
fmopt = theta(M+1:end);
    

figure(2);grid on;
subplot(211);plot(decay_to_t60(a1mopt,fs), 'o'); hold off; ylabel('T60 (s)');
subplot(212);plot(abs(fmopt), 'o'); hold off; ylabel('Frequency (Hz)');
drawnow;

%%
function [alpha] = t60_to_decay(t60,fs)
    alpha = log(1000)./t60/fs;
end

function [t60] = decay_to_t60(alpha,fs)
    t60 = log(1000)./alpha/fs;
end


end

