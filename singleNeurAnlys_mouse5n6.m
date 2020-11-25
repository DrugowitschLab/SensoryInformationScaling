function singleNeurAnlys_mouse5n6(dataset,subsample)

tic;

%% (general) settings
if nargin < 2, subsample = 'none'; end

cachefolder = ['.' filesep 'tuning_fits' filesep dataset filesep];
if ~exist(cachefolder, 'dir')
   mkdir(cachefolder)
end
TC_coef_path = [cachefolder 'TuningCoef_Combo_Subsample_' subsample '.mat'];


tcfitpath = ['.' filesep 'tuning_fits'];

% exampleids = [1 3 6];  % example neurons [non-tuned ori dir]
eventRateScale  = 30;

%% general plot settings
set(groot, 'DefaultFigureColor', 'white');
set(groot, 'defaultAxesLabelFontSizeMultiplier', 1, 'defaultAxesTitleFontSizeMultiplier', 1);
set(groot, 'defaultAxesFontSize', 6, 'defaultAxesFontName', 'Helvetica');
set(groot, 'defaultAxesFontSizeMode', 'manual');
set(groot, 'defaultTextFontSize', 6, 'defaultTextFontName', 'Helvetica');
set(groot, 'defaultAxesTickDir', 'out', 'defaultAxesTickDirMode', 'manual');
set(groot, 'defaultAxesXColor', [0 0 0], 'defaultAxesYColor', [0 0 0]);
set(groot, 'defaultAxesBox', 'off');  % overridden by plot(.)
set(groot, 'defaultAxesLayer', 'top');


%% Loading data

fprintf('Loading %s ...\n', dataset');
[sr, oris, cons, ~] = loaddata(dataset,subsample);
sr = sr* eventRateScale; % to be consistent with previous

%% loading parameters of individual tuning fits
tcfitfile = [tcfitpath filesep dataset filesep 'TuningCoef_Subsample_' subsample '.mat'];
fprintf('Loading %s ...\n', tcfitfile);
tcd = load(tcfitfile);


%%
conus           = unique(cons);
ncons           = length(conus);
conid           = 1; % doesn't matter for von mises function

vmfun.OriTun_vonMises = tcd.coef.OriTun_vonMises{conid};
vmfun.DirTun_vonMises = tcd.coef.DirTun_vonMises{conid};
vmfun.NullTun         = tcd.coef.NullFunc{conid};


for con = 1: ncons
    x_OriTun(:, :, con)    = tcd.coef.x_OriTun{con,1};
    x_DirTun(:, :, con)    = tcd.coef.x_DirTun{con,1};
    x_null(:, :, con)      = tcd.coef.x_null{con,1};
end

%% select neuron tuning based on high-contrast stim
conid = 2; 

nonsigneurons = find(tcd.coef.alpha_OriVsNull{conid} >= 0.05);
orineurons = find(tcd.coef.alpha_OriVsNull{conid} < 0.05 & ...
    tcd.coef.alpha_OriVsDir{conid} > 0.05);
dirneurons = find(tcd.coef.alpha_OriVsNull{conid} < 0.05 & ...
    tcd.coef.alpha_OriVsDir{conid} < 0.05);

P_hi.OriTun    = tcd.coef.x_OriTun{conid,1};
P_hi.DirTun    = tcd.coef.x_DirTun{conid,1};
P_hi.Null      = tcd.coef.x_null{conid,1};


%% von-mises functions are the same for both contrast

conid = 1;

vmfun.OriTun_vonMises   = tcd.coef.OriTun_vonMises{conid};
vmfun.DirTun_vonMises   = tcd.coef.DirTun_vonMises{conid};
vmfun.NullTun           = tcd.coef.NullFunc{conid};


P_lo.OriTun    = tcd.coef.x_OriTun{conid,1};
P_lo.DirTun    = tcd.coef.x_DirTun{conid,1};
P_lo.Null      = tcd.coef.x_null{conid,1};

%% compare individual fit with combo fit f_hi(\theta) = a+ b* f_lo(\theta)

sr_lo = sr(cons == conus(1), :);
oris_lo = oris(cons == conus(1), :); 

sr_hi = sr(cons == conus(2), :);
oris_hi = oris(cons == conus(2), :); 


%% Some Manual Tuning done in other works for tuning fit as well (citation to be added!)
dummy_OriTun  = log(0.5)/(cos(0.8* pi/4)- 1); 
lb_OriTun = [0     0       0       -inf];     % Hard lower bounds Model1
ub_OriTun = [inf  inf   dummy_OriTun   inf];    % Hard upper bounds

dummy_DirTun  = log(0.5)/(cos(0.45* pi/4)- 1); 
lb_DirTun = [0     0       0       -inf   0];     % Hard lower bounds Model1
ub_DirTun = [inf  inf   dummy_DirTun   inf  inf];
lb         = [-inf  -inf ];                % Hard lower bounds Model1
ub         = [ inf   inf];                % Hard upper bounds


MultStartNum    = 200;          % number of run for fmincon solver

x0 = 0.01 + 0.2* rand(1, 2);        % Starting point

nNeurons = size(sr, 2);
x_combo = cell(nNeurons, 1);
f_combo  = nan(nNeurons, 1);
for ii = 1: nNeurons             

    %% Optimizing using fmincon with multistart
    opts = optimoptions(@fmincon,'Algorithm','sqp');
    if sum(dirneurons == ii)
        lsf_combo  = @(x) fitMultip(x, vmfun, oris_lo, oris_hi, sr_lo(:, ii), sr_hi(:, ii), 'DirTun');
        prblm = createOptimProblem('fmincon','objective',lsf_combo,...
            'x0', [P_hi.DirTun(ii, :) x0], 'lb', [lb_DirTun lb], ...
            'ub', [ub_DirTun ub], 'options', opts);
    
        ls_run_flag = true;
    elseif sum(orineurons == ii)
        lsf_combo  = @(x) fitMultip(x, vmfun, oris_lo, oris_hi, sr_lo(:, ii), sr_hi(:, ii), 'OriTun');
        prblm = createOptimProblem('fmincon','objective',lsf_combo,...
            'x0', [P_hi.OriTun(ii, :) x0], 'lb', [lb_OriTun lb], ...
            'ub', [ub_OriTun ub], 'options', opts);
        ls_run_flag = true; 
    else % no need to run the code for untuned neurons
%         lsf_combo  = @(x) fitMultip(x, vmfun, oris_lo, oris_hi, sr_lo(:, ii), sr_hi(:, ii), 'null');
        ls_run_flag = false;
    end
    
    if ls_run_flag
        
        ms = MultiStart('Display', 'off');
        [x_combo{ii},f_combo(ii)] = run(ms, prblm, MultStartNum);
    
        fprintf('fitting for neuron %d with %d random start by fmincon \n', ii, MultStartNum)
    end
end

SSE_combo   = nan(nNeurons, 1);
SSE_ind     = nan(nNeurons, 1);
aic_combo   = nan(nNeurons, 1);
aic_ind     = nan(nNeurons, 1);

bic_combo   = nan(nNeurons, 1);
bic_ind     = nan(nNeurons, 1);

aicFunc  = @(SSE, k, n) (2*k+ n* log(SSE) );
bicFunc  = @(SSE, k, n) (log(n)*k+ n* log(SSE) );

nTrials = size(sr_hi, 1);


aic_win = zeros(3, 1);
bic_win = zeros(3, 1);

nDirOriNull = [length(dirneurons); length(orineurons); length(nonsigneurons)];

for ii = 1: nNeurons
    % residuals for max likelihood meathod        
    if sum(dirneurons == ii)
        SSE_combo(ii)   = f_combo(ii);
        SSE_ind(ii)     = fitVonMises(P_lo.DirTun(ii, :), oris_lo, sr_lo(:, ii), 'DirTun')+ ...
            fitVonMises(P_hi.DirTun(ii, :), oris_hi, sr_hi(:, ii), 'DirTun');
        aic_combo(ii)   = aicFunc(SSE_combo(ii), 7, nTrials);
        aic_ind(ii)     = aicFunc(SSE_ind(ii), 10, nTrials);
        
        bic_combo(ii)   = bicFunc(SSE_combo(ii), 7, nTrials);
        bic_ind(ii)     = bicFunc(SSE_ind(ii), 10, nTrials);
        
        if aic_combo(ii) < aic_ind(ii)
            aic_win(1)  = aic_win(1)+ 1;
        end
        if bic_combo(ii) < bic_ind(ii)
            bic_win(1)  = bic_win(1)+ 1;
        end
        
    elseif sum(orineurons == ii)
        SSE_combo(ii)   = f_combo(ii); 
        SSE_ind(ii)     = fitVonMises(P_lo.OriTun(ii, :), oris_lo, sr_lo(:, ii), 'OriTun')+ ...
            fitVonMises(P_hi.OriTun(ii, :), oris_hi, sr_hi(:, ii), 'OriTun');
        aic_combo(ii)   = aicFunc(SSE_combo(ii), 6, nTrials);
        aic_ind(ii)     = aicFunc(SSE_ind(ii), 8, nTrials);
        
        bic_combo(ii)   = bicFunc(SSE_combo(ii), 6, nTrials);
        bic_ind(ii)     = bicFunc(SSE_ind(ii), 8, nTrials);
        
        if aic_combo(ii) < aic_ind(ii)
            aic_win(2)  = aic_win(2)+ 1;
        end
        if bic_combo(ii) < bic_ind(ii)
            bic_win(2)  = bic_win(2)+ 1;
        end
    end
end


%% save combo fit stats

coef.SSE_combo = SSE_combo;
coef.SSE_ind = SSE_ind;

coef.aic_combo = aic_combo;
coef.aic_ind = aic_ind;

coef.bic_combo = bic_combo;
coef.bic_ind = bic_ind;

coef.x_combo = x_combo;
coef.f_combo = f_combo;

coef.aic_win = aic_win;
coef.bic_win = bic_win;

coef.dirneurons = dirneurons;
coef.orineurons = orineurons;

coef.nDirOriNull = nDirOriNull;

%% write data to file
fprintf('Writing TC fitted coeficients to %s\n', TC_coef_path);
save(TC_coef_path, 'coef');
T_total = toc;
fprintf('total running time = %0.2f hours \n', T_total/3600);

    
%% direct fitting of the means ignoring the variances
function Residual = fitMultip(P, vmfun, theta_lo, theta_hi, r_lo, r_hi, Tuning)
               % model with 4 or 5 parameters
    if strcmp(Tuning, 'DirTun')
        mu_hi = vmfun.DirTun_vonMises(P(1: 5), theta_hi);
        mu_lo = P(6)+ P(7)* vmfun.DirTun_vonMises(P(1: 5), theta_lo);
    elseif strcmp(Tuning, 'OriTun')
        mu_hi = vmfun.OriTun_vonMises(P(1: 4), theta_hi);
        mu_lo = P(5)+ P(6)* vmfun.OriTun_vonMises(P(1: 4), theta_lo);
    elseif strcmp(Tuning, 'null')
        mu_hi = vmfun.NullTun(P(1), theta_hi);
        mu_lo = P(2)+ P(3)* vmfun.NullTun(P(12), theta_lo);          
    end
    Residual  = sum((r_lo- mu_lo).^2)+ sum((r_hi- mu_hi).^2); 


function Residual = fitVonMises(P, theta, r, Tuning)
               % model with 4 or 5 parameters
    if strcmp(Tuning, 'DirTun')
        mu = P(1)+ P(2)* exp(P(3)*cos((theta-P(4))*pi/180)) + P(5)* exp(-P(3)*cos((theta-P(4))*pi/180));
    elseif strcmp(Tuning, 'OriTun')
        mu = P(1)+ P(2)* exp(P(3)*cos(2*(theta-P(4))*pi/180));
    elseif strcmp(Tuning, 'null')
        mu = repmat(P(1), length(theta), 1);            
    end
    Residual  = sum((r- mu).^2); 
