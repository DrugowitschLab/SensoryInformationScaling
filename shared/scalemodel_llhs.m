function llhs = scalemodel_llhs(Iincrsamples, Ns, llhtype)
%% returns a log-likelihood structure for the given data of the given type
%
% Possible types are 'norm' (default), 'emp', 'gig'

if nargin < 3, llhtype = 'norm'; end
if nargin < 2, Ns = 1:size(Iincrsamples, 2); end
assert(length(Ns) == size(Iincrsamples, 2));
N = max(Ns);


%% settings
empdensitydiscr = 5000;


%% generate log-likelihood functions
switch llhtype
    case 'emp'
        edp = NaN(empdensitydiscr, length(Ns));
        edi1 = NaN(1, length(Ns));
        eddi = NaN(1, length(Ns));
        fprintf('Fitting empirical densities for n =');
        for n = 1:length(Ns)
            if mod(n,50) == 0, fprintf(' %d', Ns(n)); end
            [edp(:,n), is] = ...
                ksdensity(Iincrsamples(:,n), 'NumPoints', empdensitydiscr);
            edi1(n) = is(1);
            eddi(n) = is(2)-is(1);
        end
        fprintf('\n');
        % empirical density: sum of interpolated log-probabilities
        llhfn = @(In) log(lininterp(edi1, eddi, edp, In));
        llhfnZ = @(In) log(lininterp(edi1, eddi, edp, In));
    case 'norm'
        Iincr_mu = mean(Iincrsamples, 1);
        Iincr_var = var(Iincrsamples, [], 1);
        normconst = -0.5*(log(2*pi) + log(Iincr_var));
        llhfn = @(In) -0.5*(In - Iincr_mu).^2 ./ Iincr_var;
        llhfnZ = @(In) normconst - 0.5*(In - Iincr_mu).^2 ./ Iincr_var;
    case 'gig'
        fprintf('Fitting gig for n =');
        gigp = NaN(length(Ns), 4);  %2nd dim: x0 a b p
        for n = 1:length(Ns)
            if mod(n,50) == 0, fprintf(' %d', Ns(n)); end
            gigp(n,:) = fitgig(Iincrsamples(:,n));
        end
        fprintf('\n');
        llhfn = @(In) arrayfun(@(n) gigllh(In(n), gigp(n,:)), 1:length(In));
        llhfnZ = @(In) arrayfun(@(n) gigllhZ(In(n), gigp(n,:)), 1:length(In));
    otherwise
        error('Unknown density type %s', densitytype);
end
llhs = struct(...
    'llhfn', llhfn, ...
    'llhfnZ', llhfnZ, ...
    'llhtype', llhtype, ...
    'N', N, ...
    'Ns', Ns);


function x = ifelse(condition, truex, falsex)
%% helper function for anonymous functions
if condition
    x = truex;
else
    x = falsex;
end


function vs = lininterp(x1s, dxs, vxs, xs)
%% returns linearly interpolated values
%
% vs is an N x K matrix, with each column k containing the values associated
% with x1s(k), x1s(k)+dxs(k), x1s(k)+2 dxs(k), ... The function returns the
% values for each xs(k), using linear interpolation. For extrapolation it
% returns the first/last element in vxs.
%
% It assumes x1s, dxs1, and xs to be row vectors, and returns a row vector.

[N, K] = size(vxs);
% indicies into the columns of vxs
is = (xs - x1s) ./ dxs + 1;
fis = floor(is);
% fix indicies for extrapolation
below1 = is < 1;
aboveN = is >= N;
fis(below1) = 1;    % fis = 1;   is - fis = 0
is(below1) = 1;
fis(aboveN) = N-1;  % fis = N-1; is - fis = 1
is(aboveN) = N;
% perform interpolation
vfis = vxs((0:N:((K-1)*N))+fis);  % pick row fis(k) for each column k
vs = vfis + (vxs((1:N:((K-1)*N+1))+fis) - vfis) .* (is - fis);


function theta = fitgig(x)
%% fits generalized inverse gaussian and returns parameters

opt = optimset('Display','notify');
% determin parameter bounds and initial values
x_min = min(x);
x0 = x_min - 0.01*abs(x_min);
dx = x - x0;
dx_mu = mean(dx);
dx_var = var(dx);
thetaini = [x0 (1/dx_var) (dx_mu^3/dx_var) -0.5];
thetamin = [-Inf 1e-10 1e-10 -Inf];
thetamax = [x_min Inf Inf Inf];
% perform fit
theta = fmincon(...
    @(theta) ifelse(all(theta >= thetamin) && all(theta <= thetamax), ...
                    -gigllhZ(x, theta), Inf), ...
    thetaini, [], [], [], [], thetamin, thetamax, [], opt);


function llh = gigllh(x, theta)
%% returns unnormalized likelihood of generalized inverse gaussian
dx = x - theta(1);
llh = (theta(4)-1)*sum(log(dx)) - (0.5*theta(2))*sum(dx) ...
    - (0.5*theta(3))*sum(1./dx);

function llh = gigllhZ(x, theta)
%% returns likelihood of generalized inverse gaussian
dx = x - theta(1);
llh = length(x) * (0.5*theta(4)*log(theta(2)/theta(3)) ...
    - log(2*besselk(theta(4), sqrt(theta(2)*theta(3))))) ...
    + (theta(4)-1)*sum(log(dx)) - (0.5*theta(2))*sum(dx) ...
    - (0.5*theta(3))*sum(1./dx);