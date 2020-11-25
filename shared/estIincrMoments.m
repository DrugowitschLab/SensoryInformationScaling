function [Iincr_mu, Iincr_var, Iincr] = estIincrMoments(mu, S, T, ds, popsamples, Nmax)
%% estimates the information increase moments for given dataset
%
% mu and S are the empirical f' and Sig estimates, based on T samples with
% stimulus difference ds.
%
% The function uses popsamples subpopulation samples to estimate the info
% estimation variance across different subpopulations, and returns the
% empirical moments of the information increase. They turn out the to dominate
% the overall moments, such that the additional variability arising from
% observing a limited number of trials becomes ignorable.
%
% If Nmax is given (defaults to N), then the function only computes the
% information increase moments for population sizes up to Nmax.
%
% The function returns an estimate of the mean information increase, Iincr_mu,
% and its associated uncertainty, Iincr_var.

% useparfor = usejava('desktop');
useparfor = false; % the line added for debugging

N = length(mu);
if nargin < 6, Nmax = N; end
if Nmax > N
    error('Given Nmax larger than actual population size, %d > %d', Nmax, N);
end

%% sample permutations first
% This ensures that same RNG seed gives same permutation sequence
Norder = NaN(popsamples, Nmax);
for i = 1:popsamples
    Norder(i,:) = randperm(N,Nmax);
end


%% estimate information and estimator variance for different subpopulations
Iest = NaN(popsamples, Nmax);
% avoid using parfor on cluster
if useparfor
    parfor i = 1:popsamples
        if mod(i,50) == 0; fprintf(' %d', i); end
        Iest(i,:) = empInfscaling(mu, S, Norder(i,:), T, ds);
    end
    fprintf('\n');
else
    for i = 1:popsamples
        if mod(i,50) == 0; fprintf(' %d', i); end
        Iest(i,:) = empInfscaling(mu, S, Norder(i,:), T, ds);
    end
    fprintf('\n');
end


%% compute information increase estimate and associated moments
Iincr = diff([zeros(popsamples, 1) Iest], 1, 2);
Iincr_mu = mean(Iincr, 1);
Iincr_var = var(Iincr, [], 1);
