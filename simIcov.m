function simIcov(popsize)
%% simulates population activity to compute info scaling covariances
%
% popsize is the population size of the population to simulate. The
% covariance statistics will be written to
% cov_cache/covsim_N[popsize].mat

%% settings
% (see genPopMoments for moment generation parameters)
N = 300;   % size of subpopulation
T = 1000;
shufsamples = 100;  % sample per shuffle
shufshufs = 10;     % number of shuffles
subsamples = shufsamples * shufshufs;  % number of subsamples
momsamples = 100;   % number of moments
ds = 45*pi/180;
shufpops = 10;
outfile = sprintf('cov_cache%scovsim_N%d.mat', filesep, popsize);
useparfor = true;

% generate varcovpops populations and compute covariance across shuffles
fprintf('Generating %d populations to estimate covariance due to subsampling...\n', shufpops);
shufcov = NaN(N,N,shufpops);
if useparfor
    parfor i = 1:shufpops
        fprintf('%d ', i);
        [fp, Sig, Sig0] = genPopMoments(popsize);
        if i == 1, plotMomentStats(fp, Sig, Sig0, N, T, ds); end
        In = NaN(subsamples, N);
        for j = 1:subsamples
            k = randperm(popsize,N);
            In(j,:) = empInfscaling(fp(k), Sig(k,k), 1:N, T, ds);
        end
        Insubincr = [In(:,1) diff(In,1,2)];
        shufcov(:,:,i) = cov(Insubincr);
    end
else
    for i = 1:shufpops
        fprintf('%d ', i);
        [fp, Sig, Sig0] = genPopMoments(popsize);
        if i == 1, plotMomentStats(fp, Sig, Sig0, N, T, ds); end
        In = NaN(subsamples, N);
        for j = 1:subsamples
            k = randperm(popsize,N);
            In(j,:) = empInfscaling(fp(k), Sig(k,k), 1:N, T, ds);
        end
        Insubincr = [In(:,1) diff(In,1,2)];
        shufcov(:,:,i) = cov(Insubincr);
    end
end
fprintf('\n');


%% compute variance across different suffles and moment samples
fprintf('Compute variance across %d shuffles and %d moment samples each...\n', ...
    shufsamples, momsamples);
[fp, Sig] = genPopMoments(popsize);
Sigchol = chol(Sig);
Inshuf = NaN(momsamples, shufshufs, shufsamples, N);
Insub = NaN(momsamples, shufshufs, shufsamples, N);
% generate shuffles first, to keep them the same across moment samples
kshuf = NaN(shufshufs, shufshufs, N);
ksub = NaN(shufshufs, shufshufs, N);
for i = 1:shufshufs
    k = randperm(popsize, N);
    for j = 1:shufsamples
        kshuf(i,j,:) = k(randperm(N,N));
        ksub(i,j,:) = randperm(popsize,N);
    end
end
% collect info scaling across shuffles and moment samples
if useparfor
    parfor i = 1:momsamples
        if mod(i, 10) == 0, fprintf('%d ', i); end    
        % generate population samples of full population and compute stats
        X1 = bsxfun(@plus, (fp - ds/2), randn(T, popsize) * Sigchol);
        X2 = bsxfun(@plus, (fp + ds/2), randn(T, popsize) * Sigchol);
        mu = (mean(X2) - mean(X1)) / ds;
        S = 0.5*(cov(X1)+cov(X2));
        % compute info scaling for different shuffles
        for j = 1:shufshufs
            for k = 1:shufsamples
                ks = squeeze(kshuf(j,k,:));
                Inshuf(i,j,k,:) = empInfscaling(mu(ks), S(ks,ks), 1:N, T, ds);
                ks = squeeze(ksub(j,k,:));
                Insub(i,j,k,:) = empInfscaling(mu(ks), S(ks,ks), 1:N, T, ds);
            end
        end
    end
else
    for i = 1:momsamples
        if mod(i, 10) == 0, fprintf('%d ', i); end    
        % generate population samples of full population and compute stats
        X1 = bsxfun(@plus, (fp - ds/2), randn(T, popsize) * Sigchol);
        X2 = bsxfun(@plus, (fp + ds/2), randn(T, popsize) * Sigchol);
        mu = (mean(X2) - mean(X1)) / ds;
        S = 0.5*(cov(X1)+cov(X2));
        % compute info scaling for different shuffles
        for j = 1:shufshufs
            for k = 1:shufsamples
                ks = squeeze(kshuf(j,k,:));
                Inshuf(i,j,k,:) = empInfscaling(mu(ks), S(ks,ks), 1:N, T, ds);
                ks = squeeze(ksub(j,k,:));
                Insub(i,j,k,:) = empInfscaling(mu(ks), S(ks,ks), 1:N, T, ds);
            end
        end
    end
end
% variance for different subsamples, across subsamples / moments separately
Insub = reshape(Insub, [momsamples subsamples N]); % collapse 2nd & 3rd dim
Insubincr = cat(3, Insub(:,:,1), diff(Insub,1,3));
subvar = squeeze(var(Insubincr, [], 2));    % subsampling variance
momsubvar = squeeze(var(Insubincr, [], 1)); % moment variance
% variance for different shuffles, across shuffles / moments separately
% those variances are computed as means across shufshuf's
Inshufincr = cat(4, Inshuf(:,:,:,1), diff(Inshuf,1,4));
shufvar = squeeze(mean(var(Inshufincr,[],3),2));    % shuffling variance 
momshufvar = squeeze(mean(var(Inshufincr,[],1),2)); % moment variance
% total variance across subpopulations
minsamples = min(subsamples, momsamples);
totalvar = var(reshape(Insubincr(1:minsamples, 1:minsamples, :), [], N));


%% saving simulations
fprintf('Writing data to %s ...\n', outfile);
save(outfile, 'N', 'T', 'popsize', 'shufcov', ...
    'subvar', 'momsubvar', 'shufvar', 'momshufvar', 'totalvar');


function plotMomentStats(fp, Sig, Sig0, N, T, ds)
%% plots a set of statistics for given moments, and subpopulation size

M = length(fp);
fprintf('Generating moment summary plots\n');
fprintf('Estimate info scaling ');
In_lim = NaN(10,N);
In_nolim = NaN(10,N);
for i = 1:10
    fprintf('%d ', i);
    j = randperm(M,N);
    In_nolim(i,:) = empInfscaling(fp, Sig0, j);
    In_lim(i,:) = empInfscaling(fp, Sig, j);
end
fprintf('\n');
figure('Color', 'white');  hold on;
plot(1:N, mean(In_nolim, 1), 'r-');
plot(1:N, mean(In_lim, 1), 'k-');
xlabel('Population size');
ylabel('Fisher information');

% subsample population and plot various info scaling statistics
i = randperm(M,N);
fpsub = mvnrnd(fp(i), (2/(T*ds^2)) * Sig(i,i));
Ssub = wishrnd(Sig(i,i) / (2*T-2), 2*(T-1));
[Q,D] = eig(Ssub);
D = diag(D);
[D,j] = sort(D,'descend');
Q = Q(:,j);
% variance over principal dimensions
figure('Color', 'white');
subplot(2,2,1);
plot(1:N, D);
ylabel('variance');
xlabel('principal dimension');      xlim([1 N]);
set(gca,'Box','off','YScale','log');
% alignment of fp to principal dimensions
fpsub2 = sum(fpsub.^2);
cos2theta = (fpsub * Q).^2 ./ fpsub2;
subplot(2,2,2);
plot(1:N, cumsum(cos2theta), 'k-');
set(gca,'Box','off');
ylabel('fp alignment');         ylim([0 1]);
xlabel('principal dimension');  xlim([1 N]);
In = NaN(1, N);
for n = 1:N
    Qnfpsub = Q(:,1:n)' * fpsub';
    QnSsub = Q(:,1:n)' * Ssub *  Q(:,1:n);
    % when inverting QnStest, lower-bound eigenvalues to avoid
    % numerical instabilities
    [Qn,Dn] = eig(QnSsub);
    Dninv = 1./diag(Dn);
    Dninv(diag(Dn) < 1e-15) = 0;
    In(n) = Qnfpsub' * bsxfun(@times, Qn, Dninv') * Qn' * Qnfpsub; 
end
subplot(2,2,3);
plot(1:N, In./In(end), 'k-');
xlabel('principal dimension');   xlim([1 N]);
ylabel('fraction information');  ylim([0 1]);


function [fp, Sig, Sig0] = genPopMoments(popsize)
%% generates the activity moments of a population of size popsize

%% settings
Iinf = 100;
avgfp2 = 0.008;

% generate eigenvectors of Sig0 as orthonormal random matrix
[Z,~] = qr(randn(popsize,popsize));
% eigenvalues of Sig0
sigs0 = 0.00005 + 3 ./ (1:popsize).^0.5;
% generate Sig0
Sig0 = bsxfun(@times, Z, sigs0) * Z';

% alignment of fp to Sig0
fpalign = exp(- (1:popsize) / 30) + 0.001;
fpalign = fpalign / sum(fpalign);
% generate actual fp from alignment
fp = sum(bsxfun(@times, fpalign, Z), 2)';
% rescale such that norm(fp)^2 / M = avgfp2
fp = (sqrt(avgfp2 * popsize) / norm(fp)) * fp;  

% generate Sig by adding information-limiting correlations
Sig = Sig0 + (fp'/Iinf)*fp;

