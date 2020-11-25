function [fp, Sig, Sig0] = popMoments(g, sig02, sigb2, Iinf, N, Npow)
%% generates population f' and Sig for give activity parameters
%
% Here, g is the norm of f', sig02 is the baseline variance of Sig0, and
% sigb2 is the maximum spectrum variance of Sig0. The optimal Npow gives
% the power on the 1/N scaling.

if nargin < 6, Npow = 1; end

% f' is random unit vector of length g
fp = randn(1, N);
fp = g * fp / norm(fp);
% Sig0 through random rotation and spectrum
U = RandOrthMat(N);
D = diag(sig02 + sigb2 * (1 ./ (1:N).^Npow));
Sig0 = U * D * U';
% Sig adds information-limiting components
Sig = Sig0 + (1/Iinf) * (fp' * fp);
