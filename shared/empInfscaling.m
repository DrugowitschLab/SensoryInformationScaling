function In = empInfscaling(fp, Sig, norder, T, ds)
%% return empirical info scaling for given f' and Sig using norder
%
% The function orders the data according to norder, and then computes for
% increasing population size the empirical information. It does so
% incrementally by using block matrix identities, which is more efficient
% than computing the information estimate for each subpopulation
% separately.
%
% norder is optional and defaults to 1:N.
%
% If the optional T and ds are given, the function applies bias correction
% to the empirical estimate.

if nargin < 3
    N = size(Sig, 1);
    norder = 1:N;
else
    N = length(norder);
end
In = NaN(1, N);
invSig = NaN(N, N);  % to be filled incrementally

% reorder elements in fp and Sig
fp = fp(norder);
Sig = Sig(norder, norder);

% start with n=1 case
invSig(1,1) = 1 / Sig(1, 1);
In(1) = invSig(1, 1) * fp(1)^2;
if isinf(invSig(1,1))
    invSig(1,1) = 0;
    In(1) = 0;
end
% recursively handle n=2..N cases
for n = 2:N
    % Block matrix pseudo-inv. (Rohde, 1995) to add n'th element to invSig
    a = Sig(1:(n-1), n);
    b = Sig(n, n);
    invSiga = invSig(1:(n-1), 1:(n-1)) * a;
    c = 1 / (b - a' * invSiga);
    if isinf(c)
        % by Rohde (1995), use pseudo-inv. c = 0 for b - a' invSiga = 0
        invSig(1:(n-1), 1:(n-1)) = invSig(1:(n-1), 1:(n-1));
        invSig(1:n, n) = 0;
        invSig(n, 1:(n-1)) = 0;
        In(n) = In(n-1);
    else
        d = - c * invSiga;
        invSig(1:(n-1), 1:(n-1)) = invSig(1:(n-1), 1:(n-1)) + c * (invSiga * invSiga');
        invSig(1:(n-1), n) = d;
        invSig(n, 1:(n-1)) = d';
        invSig(n, n) = c;
        % update information measure by increment
        res = fp(1:(n-1)) * invSiga - fp(n);
        In(n) = In(n-1) + c * res^2;
    end
end

% apply bias correction, if requested
if nargin == 5
    In = (2*T-(1:N)-3)/(2*T-2) .* In - (2/(T*ds^2)) * (1:N);
end
