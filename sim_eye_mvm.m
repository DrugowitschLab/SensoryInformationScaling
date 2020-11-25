function sim_eye_mvm()
%% plots figure showing impact of saturating dF/F on simulated data
tic

%% general settings
addpath('..')


%%
pn_list = repmat(0.1: 0.1: 0.9, 1, 3);
% pn_list = [0.9 0.9 0.9];
pt_list = [0.3* ones(1, 9) 0.6* ones(1, 9) ones(1, 9)];
% pt_list = [0.3 0.6 1];

%% model parameters
% network parameters
P = 32;
par = struct(...
    'N', 1000, ...    % number of neurons
    'P', P, ...       % pixels in image (one side)
    'sig', P/5, ...   % gaussian envelope sd
    'lam', P/1.5, ... % pref. wavelength
    'phi', 0, ...     % pref. spatial phase
    'c', 1, ...       % Michelson contrast
    'g', 10, ...      % tuning amplitude
    'sig_a', sqrt(2), ... % tuning ampl. variability
    'sig0', 0.11, ... % input noise (prop. of pixel range)
    'ds', 15 * pi / 180, ...
    'T', 20000);
% simulation parameters
N = par.N;
T = par.T;
theta_sim = 0;
ds = 5 * pi / 180;  % in rad

%% generate population activity
par = OrientedPopulation(par);
fprintf('Simulating %d trials...\n', par.T);
[r1, r2] = simulateTrials(par, theta_sim, ds, 0.5);

% compute info scaling for different types of saturation
pt = 0; 
pn = 0; 

fprintf('Information scaling for spike count data ...\n');
avgInfoScaling(r1, r2, ds, pt, pn);

parfor ii = 1: length(pn_list)
%%
    pn = pn_list(ii);   % fraction of neurons having zero activity
    pt = pt_list(ii);   % fraction of trials having zero activity

    Mrnd1 = rnd_eyemvn(pn, pt, T, N);
    Mrnd2 = rnd_eyemvn(pn, pt, T, N);

    fprintf('Information scaling for spike count data -- eye movement ...\n');
    rem1 = r1.* Mrnd1;
    rem2 = r2.* Mrnd2;

    avgInfoScaling(rem1, rem2, ds, pt, pn);
end

toc

function out = rnd_eyemvn(pn, pt, T, N)

out = ones(T, N);

z_t = rand(T, 1) < pt; 
trial_id = find(z_t == 1);

for ii = 1: length(trial_id)
    z_n = rand(N, 1) < pn; 
    out(trial_id(ii), z_n) = 0; 
end


function par = OrientedPopulation(par)
% generates a population of neurons with parameters par

N = par.N; 
P = par.P;    
alpha = par.sig_a^2 + 1;
par.an = par.g * lognrnd(-log(sqrt(alpha)), sqrt(log(alpha)), N, 1);
ThetaPref = linspace(-pi, pi* (N - 1)/ N, N);
J = NaN(P, P, N);
for ii = 1: N
    J(:, :, ii) = pixelizedGarbor(par, ThetaPref(ii), [0 0]);
end
par.F = standardizevector(reshape(J, P^2, [])); % compute neural filters


function J = pixelizedGarbor(par, theta, mu)
% gernates P x P pixelated Garbor of orientation theta with parameters par,
% centered on mu = [x y] (in pixels);

sig2 = -0.5/ par.sig^2;
P = par.P;
cosTheta = 2*pi *cos(theta)/ par.lam;    
sinTheta = 2* pi *sin(theta)/ par.lam;    
x = linspace(-P/2 + 0.5, P/2 - 0.5, P) + mu(1); 
y = linspace(-P/2 + 0.5, P/2 - 0.5, P) + mu(2);
J = NaN(P, P);
for jj = 1: P  
    J(:, jj) = par.c * ...
        exp(sig2 * (x.^2 + y(jj)^2)) .* ...
        cos(x * cosTheta + y(jj) * sinTheta + par.phi);    
end


function [r1, r2] = simulateTrials(par, theta, ds, dur)
%% simulates a set of trials of two orientations

r1 = NaN(par.T, par.N);
r2 = NaN(par.T, par.N);
theta1 = theta - ds/2;
theta2 = theta + ds/2;
P2 = par.P^2;
sig0 = par.sig0;
ReLU = @(x) (x >= 0) .* x;
for i = 1:par.T
    % response to theta1
    Jtheta = standardizevector(...
        reshape(pixelizedGarbor(par, theta1, [0 0]), P2, []));
    r1(i,:) = poissrnd(ReLU(dur * par.an .* ...
                            (par.F' * (Jtheta + sig0 * randn(P2, 1)))));
    % response to theta2
    Jtheta = standardizevector(...
        reshape(pixelizedGarbor(par, theta2, [0 0]), P2, []));
    r2(i,:) = poissrnd(ReLU(dur * par.an .* ...
                            (par.F' * (Jtheta + sig0 * randn(P2, 1)))));
end


function ftheta = popresp(par, theta, mu)
% computes mean population response of population par to stim theta

P = sqrt(size(par.F, 1));
Jtheta = pixelizedGarbor(par, theta, mu);
JthetaNorm = standardizevector(reshape(Jtheta, P^2, []));
ftheta = par.an.*(par.F'* JthetaNorm);
ftheta(ftheta < 0) = 0;


function SigTheta = poprespcov(par, theta)
% computes population response covariance to stim theta

P = sqrt(size(par.F, 1));    
Jtheta = pixelizedGarbor(par, theta, [0 0]);
JthetaNorm = standardizevector(reshape(Jtheta, P^2, []));
sig02 = par.sig0^2;
ftheta = par.an'.*(JthetaNorm'* par.F);
FF = par.F'*par.F;
ftheta(ftheta < 0) = 0;
FF(FF < 0) = 0;
SigTheta = sig02 * (par.an*par.an') .* FF + diag(ftheta);            

       
function Y = standardizevector(X)
% standardize vector X or columns of matrix X
[n, m] = size(X);
if n==1 || m ==1
    Y = X- mean(X);
    Y = Y/norm(Y);
else
    Y = NaN(size(X));
    for jj = 1: m
        Y(:, jj) = standardizevector(X(:, jj));
    end
end   


function [Imu, Isd] = avgInfoScaling(r1, r2, ds, pt, pn)
%% computes average info scaling across different orders
%
% performs trial-identify shuffling if [optimal] shuf is given and true
[T, N] = size(r1);


% pop statistics
fp = (mean(r1,1) - mean(r2,1)) / ds;
Sig = 0.5*(cov(r1) + cov(r2));

popsamples = 1000; 
[Iincr_mu, Iincr_var, Iincr_samples] = ...
            estIincrMoments(fp, Sig, T, ds, popsamples, N);

Isd = [0 sqrt(cumsum(Iincr_var))];
Imu = [0 cumsum(Iincr_mu)];

% saving the moment file
outfile = ['.' filesep 'moment_cache' filesep 'sim_eye_jitter_pn_' num2str(pn) '_pt_' num2str(pt) '_N_' num2str(N) '.mat'];
fprintf('Writing data to %s\n', outfile);
% Iincr_samples = Iincr;
mu = fp;
S = Sig; 
subdims = N;
% 
save(outfile, 'Iincr_mu', 'Iincr_var', 'Iincr_samples', ...
    'T', 'mu', 'S', 'ds', 'subdims');