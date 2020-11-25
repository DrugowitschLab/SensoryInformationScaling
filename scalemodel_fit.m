function scalemodel_fit(momentfile, llhtype)
%% fits different information scaling models to moments from moment file
%
% momentfile is a filename in 'moment_cache'.
%
% llhtype is one accepted by scalemodel_llhs
%
% The function uses slice sampling to fit the curves.


%% settings
if nargin < 2, llhtype = 'norm'; end
% ML parameters
opt = optimset('Display', 'notify', 'TolX', 1e-10, 'TolFun', 1e-10);
% slice sampling parameters
ssburnin = 100;
sssamples = 100000; 
ssthin = 10;
sschains = 4;
    

%% loading data
infile = ['moment_cache' filesep momentfile '.mat'];
fprintf('Loading %s...\n', infile);
d = load(infile);
if isfield(d, 'Ns'), Ns = d.Ns;
else, Ns = 1:length(d.Iincr_mu); end
N = max(Ns);
Iest_mu = cumsum(d.Iincr_mu);
Iest_var = cumsum(d.Iincr_var);


%% initialize info scaling functions and likelihoods
sms = {scalemodel_lin(d.Iincr_samples, Ns), ...
       scalemodel_limlin(d.Iincr_samples, Ns), ...
       scalemodel_limexp(d.Iincr_samples, Ns)};
llhs = scalemodel_llhs(d.Iincr_samples, Ns, llhtype);


%% ML fits
fprintf('Fitting the following models:\n');
for modeli = 1:length(sms)
    sm = sms{modeli};
    fprintf('- %s\n', sm.name);
    sms{modeli}.ml = llhfit(sm, llhs, opt);
    fprintf('\n');
end


%% posterior samples - initialized with ML fits
fprintf('Drawing %d(+%d burnin; thin %d, %d chains) samples for following models:\n', ...
    sssamples, ssburnin, ssthin, sschains);
for modeli = 1:length(sms)
    sm = sms{modeli};
    fprintf('- %s\n', sm.name);
    sms{modeli}.mc = llhsample(sm, llhs, sssamples, ssburnin, ssthin, sschains);
    fprintf('\n');
end


%% write data to file
outfile = ['fits' filesep momentfile '_' llhtype '.mat'];
fprintf('Writing fit data to %s...\n', outfile);
save(outfile, 'sms');


%% plot fit results (if not on cluster)
if usejava('desktop')
    plotSmFits([momentfile '_' llhtype]);
end



function s = llhfit(sm, llhs, opt)
%% performs max-llh fit with given parameters, and return stats

x = fmincon(@(x) - sum(llhs.llhfn(sm.Iincrfn(x, llhs.Ns))) - sm.logp(x), ...
    sm.pini, [], [], [], [], sm.pmin, sm.pmax, [], opt);
llh = sum(llhs.llhfnZ(sm.Iincrfn(x, llhs.Ns))) + sm.logp(x);
bic = log(length(llhs.Ns))*length(sm.pini) - 2*llh;
aic = 2*length(sm.pini) - 2*llh;

s = struct('llh', llh, 'p', x, 'bic', bic, 'aic', aic);
for i = 1:length(sm.pnames)
    fprintf('  %-9s = %f\n', sm.pnames{i}, x(i));
end
fprintf('  llh       = %f\n', llh);
fprintf(' ~bic       = %f\n', bic);
fprintf('  aic       = %f\n', aic);


function s = llhsample(sm, llhs, sssamples, ssburnin, ssthin, sschains)
%% slice-samples from the posterior, and computes some sample stats

% bounded likelihood function for sampling
boundedllhfn = @(x) ifelse(all(x >= sm.pmin) && all(x <= sm.pmax), ...
    sum(llhs.llhfn(sm.Iincrfn(x, llhs.Ns)))+sm.logp(x), -Inf);
ss = NaN(sssamples, length(sm.pnames), sschains);
pini = sm.pini;
pw = sm.pw;
% avoid using parfor on cluster
if usejava('desktop')
    parfor i = 1:sschains
        ss(:,:,i) = slicesample(pini, sssamples, 'logpdf', boundedllhfn, ...
            'burnin', ssburnin, 'thin', ssthin, 'width', pw);
    end
else
    for i = 1:sschains
        ss(:,:,i) = slicesample(pini, sssamples, 'logpdf', boundedllhfn, ...
            'burnin', ssburnin, 'thin', ssthin, 'width', pw);
    end
end

% Gelman-Rubin potential scale reduction factor
chainEparams = mean(ss, 1);                    % 1 x params x chains
chainvarparams = var(ss, [], 1);               % 1 x params x chains
Echainvarparams = mean(chainvarparams, 3);     % 1 x params
varchainEparams = var(chainEparams, [], 3);    % 1 x params
rhat = sqrt(((1-1/sssamples)*Echainvarparams+varchainEparams)./Echainvarparams);

% per-n likelihood for each parameter sample across chains
ssllhn = NaN(sssamples*sschains, length(llhs.Ns));
for i = 1:sschains
    for j = 1:sssamples
        ssllhn((i-1)*sssamples+j,:) = llhs.llhfnZ(sm.Iincrfn(ss(j,:,i), llhs.Ns));
    end
end
ssllh = sum(ssllhn, 2)';

% information criteria
Eparams = mean(chainEparams,3);
Eparamllh = sum(llhs.llhfnZ(sm.Iincrfn(Eparams, llhs.Ns)));
pDIC = 2*(Eparamllh - mean(ssllh));
pDICalt = 2*var(ssllh);
DIC = -2*Eparamllh + 2*pDIC;
DICalt = -2*Eparamllh + 2*pDICalt;

lppdn = log(mean(exp(ssllhn),1));
lppd = sum(lppdn);
pWAIC1 = 2*sum(lppdn - mean(ssllhn,1));
pWAIC2 = sum(var(ssllhn,[],1));
WAIC1 = -2*lppd + 2*pWAIC1;
WAIC2 = -2*lppd + 2*pWAIC2;

% store and output summary
s = struct('Eparams', Eparams, 'Eparamllh', Eparamllh, 'lppd', lppd, ...
    'ss', ss, 'rhat', rhat, ...
    'pDIC', pDIC, 'pDICalt', pDICalt, 'DIC', DIC, 'DICalt', DICalt, ...
    'pWAIC1', pWAIC1, 'pWAIC2', pWAIC2, 'WAIC1', WAIC1, 'WAIC2', WAIC2);
for i = 1:length(sm.pnames)
    fprintf('  E%-8s = %f (rhat=%6.3f)\n', sm.pnames{i}, Eparams(i), rhat(i)); 
end
fprintf('  Eparamllh = %f\n', Eparamllh);
fprintf('  lppd      = %f\n', lppd);
fprintf('  pDIC      = %f\n', pDIC);
fprintf('  pDICalt   = %f\n', pDICalt);
fprintf('  DIC       = %f\n', DIC);
fprintf('  DICalt    = %f\n', DICalt);
fprintf('  pWAIC1    = %f\n', pWAIC1);
fprintf('  pWAIC2    = %f\n', pWAIC2);
fprintf('  WAIC1     = %f\n', WAIC1);
fprintf('  WAIC2     = %f\n', WAIC1);


function x = ifelse(condition, truex, falsex)
%% helper function for anonymous functions
if condition
    x = truex;
else
    x = falsex;
end
