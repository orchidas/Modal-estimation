function [fmhat, a1mhat, S, nm] = hvmodel_freqs_decay(ir, nh, p, fs, kappa)

%Hankel-Vandermonde modal estimation
%nh - hankel matrix dimensions
%ir - impulse response
%p - offset

if nargin < 5
    kappa = -80;
end


% form offset Hankel matricies
H = hankel(ir(1:nh), ir(nh+(0:nh-1)));
K = hankel(ir(p+(1:nh)), ir(p+nh+(0:nh-1)));

% generate Hankel matrix pseudoinverse
[U, S, Ut] = svd(H);
nm = sum(diag(S) > 10^(kappa/20)*S(1,1));
Hinv = Ut(:,1:nm)*diag(1./diag(S(1:nm,1:nm)))*U(:,1:nm)';

% find offset Hankel matrix generalized eigenstructure
[W, D, Wt] = eig(K*Hinv);

% compute mode parameters, impulse response
nmode = sum(abs(diag(D)) > 10^(kappa/20)*abs(D(1,1)));
[~,knee_point] = knee_pt(20*log10(diag(S)/S(1,1))); 
nmode = min(knee_point, nmode);

lambda = diag(D(1:nmode,1:nmode));
[fmhat, order] = sort(angle(lambda)*fs/(2*pi*p));
a1mhat = abs(lambda(order)).^(1/p);
%a1mhat = abs(log(abs(lambda(order)))/p);

end

