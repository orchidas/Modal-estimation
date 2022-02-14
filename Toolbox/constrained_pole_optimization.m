function [fmhat_opt, a1mhat_opt, niter, nfunc] = constrained_pole_optimization(ir, fs, fmhat, a1mhat, fL, fR)
%%
% Constrained pole optimization according to Maestre et al. "Constrained
% Pole Optimization for Modal Reverberation" - DAFx 2018
% Outputs:
% fmhat_opt, a1mhat_opt     - optimized mode frequencies and pole radii
% Inputs :
% ir        - impulse response
% fs        - sampling rate
% fmhat     - pole frequencies
% a1mhat    - pole radius
%             (no overlap)
% fL        - lower frequency of band (Hz)
% fR        - upper frequency of band (Hz)
% Author :    Orchisama Das, 2020
%% 

M = length(fmhat);
indexb = find((fmhat > fL) & (fmhat < fR)); %index of frequencies that are within the current frequency band
nindexb = setdiff(1:M, indexb);     %index of overlapping frequencies
U = length(indexb);

fftSize = 4096;
H = fft(ir, fftSize);
wL = 2*pi*fL/fs;
wR = 2*pi*fR/fs;
wk = linspace(0, pi, fftSize/2); 
[v1,kL] = min(abs(wk - wL));
[v2, kR] = min(abs(wR - wk));


% optimization cost function
function err = costfn(theta)
    fm = theta(1:U)/(2*pi) * fs;
    a1m = 1 - exp(-theta(U+1:end));
    Q = zeros(K,2*U);
    for u = 1:U
        r0 = parallel_biquad_response([1;0], fm(u), a1m(u), fs, wk(kL:kR));
        Q(:,2*u-1) = r0;
        Q(:,2*u+1) = exp(-1j*wk(kL:kR)).*r0;
    end
    gb = Q\(h-v);
    % plot for sanity check
    hhat = Q*gb + v;
    figure(3); 
    semilogx(wk(kL:kR)*fs/(2*pi), 20*log10(abs([h, hhat]))); grid;
    xlabel('Freq (Hz)'); ylabel('Magnitude (dB)');xlim([20,20000]);
    drawnow;
    
    err = norm(Q*gb + v - h).^2;
end

% measurement vector
h = H(kL:kR);
K = kR - kL + 1;
Q0 = zeros(K,2*M);
for i = 1:M
    r0 = parallel_biquad_response([1;0], fmhat(i), a1mhat(i), fs, wk(kL:kR));
    Q0(:,2*i-1) = r0;
    Q0(:,2*i+1) = exp(-1j*wk(kL:kR)).*r0;
end
g = Q0\h; 
go = [g(nindexb); g(M+nindexb)];
v = parallel_biquad_response(go, fmhat(nindexb), a1mhat(nindexb), fs ,wk(kL:kR)).';

% initialize parameter vector for optimization
theta0 = [2*pi*fmhat(indexb)/fs; -log(1-a1mhat(indexb))];  

% constraints - we are only optimizing non-overlapping modes
e = ones(U,1);
A = [full(spdiags([e, -e],0:1,U,U)), zeros(U,U)]; %f_u-1 < f_u < f_u+1
b = zeros(U,1);
% thresholds for optimization
fth = (fR-fL)/10; %Hz
a1th = 0.5;
sth = -log(1-a1th);
lb = [2*pi*(fmhat(indexb)-fth)/fs; theta0(U+1:end)-sth];
ub = [2*pi*(fmhat(indexb)+fth)/fs; theta0(U+1:end) + sth];

options = optimoptions('fmincon','Algorithm','sqp','MaxIterations',300,'Display','iter');
[theta, val, exitflag, output] = fmincon(@costfn, theta0, A,b,[],[],lb,ub,[],options);
niter = output.iterations;
nfunc = output.funcCount;

% order after optimizing
[fmhat_opt, sind] = sort([theta(1:U)/(2*pi)*fs;fmhat(nindexb)]);
a1mhat_opt = [1-exp(-theta(U+1:end)); a1mhat(nindexb)];
a1mhat_opt = a1mhat_opt(sind);

end

