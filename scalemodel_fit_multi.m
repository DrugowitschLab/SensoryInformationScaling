function scalemodel_fit_multi(varargin)
%% fits different information scaling models to moments from multiple files
%
% The function is called as
% scalemodel_fit_multi(outfile, momentfile1, momentfile2, ...)
%
% Where the outfile is written to 'fit' and the momentfiles are read from
% 'moment_cache'.
%
% The llhtype of the fit is fixed to 'norm'.
%
% The function uses slice sampling to jointly fit the curves for all moment
% files.


%% settings
if nargin < 2
    error('Require at least two arguments');
end
outfilebase = varargin{1};
momentfiles = varargin(2:end);
K = nargin - 1;
smfns = {@scalemodel_lin, @scalemodel_limlin, @scalemodel_limexp};
smn = length(smfns);
% ML parameters
opt = optimset('Display', 'notify');
llhtype = 'norm';
% slice sampling parameters
ssburnin = 100;
sssamples = 100000;
ssthin = 10;
sschains = 4;
    

%% loading data
ds = cell(1, K);
Ns = cell(1, K);
Iest_mu = [];
Iest_var = [];
fprintf('Loading moments data...\n');
for k = 1:K
    infile = ['moment_cache' filesep momentfiles{k}];
    fprintf('%s\n', infile);
    ds{k} = load(infile);
    if isfield(ds{k}, 'Ns'), Ns{k} = ds{k}.Ns;
    else, Ns{k} = 1:length(ds{k}.Iincr_mu); end
    Iest_mu = cat(2, Iest_mu, cumsum(ds{k}.Iincr_mu));
    Iest_var = cat(2, Iest_var, cumsum(ds{k}.Iincr_var));
end
fprintf('\n');


%% initialize info scaling functions and likelihoods, priors
% first initialize for invidiual datasets to get prior settings
prims = cell(1, smn);
priss = cell(1, smn);
pinis = cell(1, smn);
pws = cell(1, smn);
sms = cell(smn, K);
llhs = cell(1, K);
for k = 1:K
    for smi = 1:smn
        sms{smi,k} = smfns{smi}(ds{k}.Iincr_samples, Ns{k});
        if k == 1
            prims{smi} = sms{smi,k}.prim;
            priss{smi} = sms{smi,k}.pris;
            pinis{smi} = sms{smi,k}.pini;
            pws{smi} = sms{smi,k}.pw;
        else
            prims{smi} = cat(1, prims{smi}, sms{smi,k}.prim);
            priss{smi} = cat(1, priss{smi}, sms{smi,k}.pris);
            pinis{smi} = cat(1, pinis{smi}, sms{smi,k}.pini);
            pws{smi} = cat(1, pws{smi}, sms{smi,k}.pw);
        end
    end
    llhs{k} = scalemodel_llhs(ds{k}.Iincr_samples, Ns{k}, llhtype);
end
% prior is uses average prior settings across datasets
fits = cell(1, smn);
for smi = 1:smn
    prim = mean(prims{smi},1);
    pris = mean(priss{smi},1);
    % not using struct(.) to avoid struct array creation
    fits{smi}.name = sms{smi,1}.name;
    fits{smi}.pini = mean(pinis{smi},1);
    fits{smi}.pw = mean(pws{smi},1);
    fits{smi}.pmin = sms{smi,1}.pmin;
    fits{smi}.pmax = sms{smi,1}.pmax;
    fits{smi}.prim = prim;
    fits{smi}.pris = pris;
    fits{smi}.logp = @(x) -sum(log(1 + ((x - prim)./pris).^2));
    fits{smi}.Ns = Ns;
    fits{smi}.pnames = sms{smi,1}.pnames;
end


%% ML fits
fprintf('Fitting the following models:\n');
for smi = 1:smn
    fprintf('- %s\n', fits{smi}.name);
    fits{smi}.ml = llhfit(fits{smi}, sms(smi,:), llhs, opt);
    fprintf('\n');
end


%% posterior samples - initialized with ML fits
fprintf('Drawing %d(+%d burnin; thin %d, %d chains) samples for following models:\n', ...
    sssamples, ssburnin, ssthin, sschains);
for smi = 1:smn
    fprintf('- %s\n', fits{smi}.name);
    fits{smi}.mc = llhsample(fits{smi}, sms(smi,:), llhs, ...
        sssamples, ssburnin, ssthin, sschains);
    fprintf('\n');
end


%% write data to file
outfile = ['fits' filesep outfilebase '.mat'];
fprintf('Writing fit data to %s...\n', outfile);
save(outfile, 'sms', 'fits');


%% plot fit results (if not on cluster)
%if usejava('desktop')
%    plotSmFits([momentfile '_' llhtype]);
%end



function s = llhfit(fitd, sm, llhs, opt)
%% performs max-llh fit with given parameters, and return stats

x = fmincon(@(x) - sum(combinedllhs(x, sm, llhs, fitd.Ns, false)) - fitd.logp(x), ...
    fitd.pini, [], [], [], [], fitd.pmin, fitd.pmax, [], opt);
llh = sum(combinedllhs(x, sm, llhs, fitd.Ns, true)) + fitd.logp(x);
bic = log(sum(cellfun(@length, fitd.Ns)))*length(fitd.pini) - 2*llh;
aic = 2*length(fitd.pini) - 2*llh;

s = struct('llh', llh, 'p', x, 'bic', bic, 'aic', aic);
for i = 1:length(fitd.pnames)
    fprintf('  %-9s = %f\n', fitd.pnames{i}, x(i));
end
fprintf('  llh       = %f\n', llh);
fprintf(' ~bic       = %f\n', bic);
fprintf('  aic       = %f\n', aic);


function s = llhsample(fitd, sm, llhs, sssamples, ssburnin, ssthin, sschains)
%% slice-samples from the posterior, and computes some sample stats

llhfn = @(x) sum(combinedllhs(x, sm, llhs, fitd.Ns, false))+fitd.logp(x);
% bounded likelihood function for sampling
boundedllhfn = @(x) ifelse(all(x >= fitd.pmin) && all(x <= fitd.pmax), ...
    llhfn(x), -Inf);
ss = NaN(sssamples, length(fitd.pnames), sschains);
pini = fitd.pini;
pw = fitd.pw;
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
Nlens = cellfun(@length, fitd.Ns);
ssllhn = NaN(sssamples*sschains, sum(Nlens));
for i = 1:sschains
    for j = 1:sssamples
        ssllhn((i-1)*sssamples+j,:) = ...
            combinedllhs(ss(j,:,i), sm, llhs, fitd.Ns, true);
    end
end
ssllh = sum(ssllhn, 2)';

% information criteria
Eparams = mean(chainEparams,3);
Eparamllh = sum(combinedllhs(Eparams, sm, llhs, fitd.Ns, false));
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
for i = 1:length(fitd.pnames)
    fprintf('  E%-8s = %f (rhat=%6.3f)\n', fitd.pnames{i}, Eparams(i), rhat(i)); 
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


function l = combinedllhs(x, sm, llhs, Ns, normalize)
%% returns the combined vector of likelihoods for parameters x
Nlens = cellfun(@length, Ns);
l = NaN(1, sum(Nlens));
Ni = [0 cumsum(Nlens)];
if normalize
    for k = 1:length(Ns)
        l((Ni(k)+1):Ni(k+1)) = llhs{k}.llhfnZ(sm{k}.Iincrfn(x, Ns{k}));
    end
else
    for k = 1:length(Ns)
        l((Ni(k)+1):Ni(k+1)) = llhs{k}.llhfn(sm{k}.Iincrfn(x, Ns{k}));
    end
end


function x = ifelse(condition, truex, falsex)
%% helper function for anonymous functions
if condition
    x = truex;
else
    x = falsex;
end
