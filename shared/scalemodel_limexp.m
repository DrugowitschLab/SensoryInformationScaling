function m = scalemodel_limexp(Iincrsamples, Ns)
%% returns structure defining asymptotically limited exp-linear scaling model

if nargin < 2, Ns = 1:size(Iincrsamples, 2); end
assert(length(Ns) == size(Iincrsamples,2));

% Iincr sample statistics
Iincr_mu = mean(Iincrsamples);
N = max(Ns);
Iinfini = max(0,max(cumsum(Iincr_mu)));
cini = max(0,mean(Iincr_mu./diff([0 Ns])));
tauini = 0;
Iinfsd = max(1, 10*Iinfini);
csd = 10*(cini+0.5);

% info scaling functions (N is fixed, as it scales tau)
In = @(x, Ns) Inlimexp(x, Ns, N);
Iincr = @(x, Ns) diff([0 Inlimexp(x, Ns, N)]);
% prior is student's T prior on c, Iinf, empirical mean, SD = 10 x mean,
% student's T prior on tau, zero mean, SD = 1
logp = @(x) -log(1+((x(1)-Iinfini)/Iinfsd)^2) ...
            -log(1+((x(2)-cini)/csd)^2) ...
            -log(1+x(3)^2);

% model properties
m = struct(...
    'name', 'limexp', ...
    'col', [0 0.8 0], ...
    'N', N, ...
    'Ns', Ns, ...
    'Infn', In, ...
    'Iincrfn', Iincr, ...
    'nfracfn', @(x, a) nfrac(x, a, N), ...
    'logp', logp, ...
    'pmin', [0 0 0], ...
    'pmax', [Inf Inf Inf], ...
    'pini', [Iinfini cini tauini], ...
    'pw', [(Iinfsd/5) (csd/20) 10], ...
    'prim', [Iinfini cini 0], ...
    'pris', [Iinfsd csd 1]);
m.pnames = {'Iinf', 'c', 'tau'};


function In = Inlimexp(x, Ns, N)
%% x = [Iinf, c, tau]
if x(3) < 0.01
    I0 = x(2) * Ns;
else
    I0 = x(2) * (Ns + x(3)*N*(exp(-Ns/(x(3)*N)) - 1));
end
% realmin avoids 0 / 0 is Iinf = 0 and I0 = 0
In = I0 .* x(1) ./ max(realmin, x(1) + I0);


function n = nfrac(x, a, N)
%% returns where In reaches a% of Iinf
if x(3) < 0.01
    n = (a./(1-a))*(x(1)/max(realmin,x(2)));
else
    j = length(a);
    n = NaN(size(a));
    Ntau = N * x(3);
    for i = 1:j
        % find a% if Iinf by root finding
        n(i) = fzero(@(n) a(i)*x(1) - ...
            1/(1/(x(2)*(n+Ntau*(exp(-n/Ntau)-1)))+1/x(1)), ...
            (a(i)/(1-a(i)))*(x(1)/max(realmin,x(2))));
    end
end