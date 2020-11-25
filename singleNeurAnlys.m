function singleNeurAnlys(dataset,subsample)
%% plots running speed statistics for given dataset

addpath('shared')

%% load and pre-process data
if nargin < 2, subsample = 'none'; end
cachefolder     = ['.' filesep 'tuning_fits' filesep dataset filesep];

if ~exist(cachefolder, 'dir')
   mkdir(cachefolder)
end

[sr, oris, cons, spds] = loaddata(dataset,subsample);
uoris = unique(oris);
ucons = unique(cons);
orin = length(uoris);
conn = length(ucons);
spds = spds * 100; % turn into cm/s


%% initializations for tuning curve fitting

% addpath(genpath('shared'));

TC_coef_path = [cachefolder 'TuningCoef_Subsample_' subsample '.mat'];

eventRateScale  = 30;          % the number that we got to get event rate
MultStartNum    = 200;          % number of run for fmincon solver

% Three models: Orientation Tuning; Direction Tuning and Null model with no
% tuning with their # of parameters in the model
OriTun_vonMises = @(P, theta) P(1)+ P(2)* exp(P(3)*cos(2*(theta-P(4))*pi/180));
nParam_OriTun  = 4; % number of parameters in orientation tuning von mises function

DirTun_vonMises = @(P, theta) P(1)+ P(2)* exp(P(3)*cos((theta-P(4))*pi/180)) + P(5)* exp(-P(3)*cos((theta-P(4))*pi/180));    
nParam_DirTun  = 5;

NullFunc = @(P, theta) repmat(P(1), length(theta), 1);
nParam_null  = 1;                           
                         
Theta_plt = linspace(uoris(1), uoris(end), 91);   % orientation in degree for ploting tuning curves

x0_OriTun = 0.01 + 0.2* rand(1, nParam_OriTun);    % Starting point
x0_DirTun = 0.01 + 0.2* rand(1, nParam_DirTun);    % Starting point
x0_null = 0.01 + 0.2* rand(1, nParam_null);        % Starting point

nNeurons = size(sr, 2);
rnd_plt = randperm(nNeurons, 16); % plotting the tuning curves for 16 random neurons

pval05  = 0.05;
pval01  = 0.01;
pval001 = 0.001;

%% 
for coni = 1: conn  % different levels of contrasts

    indx_cont = cons==ucons(coni); 

    averageSpkRateOverTrial = sr(indx_cont, :)* eventRateScale;       % event rate
    orientationForTrials    = oris(indx_cont);           % stim orientations
    orientations            = unique(orientationForTrials);

    [nTrial, nNeurons]      = size(averageSpkRateOverTrial);
    x_OriTun  = zeros(nNeurons, nParam_OriTun);
    f_OriTun  = zeros(nNeurons, 1);
    
    x_DirTun  = zeros(nNeurons, nParam_DirTun);
    f_DirTun  = zeros(nNeurons, 1);
    
    x_null  = zeros(nNeurons, nParam_null);
    f_null  = zeros(nNeurons, 1);

    SSE_OriTun      = zeros(nNeurons, 1);
    SSE_DirTun      = zeros(nNeurons, 1);
    SSE_null        = zeros(nNeurons, 1);
    
    R2_OriTun       = zeros(nNeurons, 1);
    R2_DirTun       = zeros(nNeurons, 1);    
    R2_null         = zeros(nNeurons, 1);
    
    %% Some Manual Tuning done in other works for tuning fit as well (citation to be added!)
    dummy_OriTun  = log(0.5)/(cos(0.8* pi/4)- 1); 
    lb1_OriTun = [0     0       0       -inf];     % Hard lower bounds Model1
    ub1_OriTun = [inf  inf   dummy_OriTun   inf];    % Hard upper bounds

    dummy_DirTun  = log(0.5)/(cos(0.45* pi/4)- 1); 
    lb1_DirTun = [0     0       0       -inf   0];     % Hard lower bounds Model1
    ub1_DirTun = [inf  inf   dummy_DirTun   inf  inf];
    
    for ii = 1: nNeurons             

        cur_neuron    = averageSpkRateOverTrial(:, ii);       

        %% Optimizing using fmincon with multistart
        opts = optimoptions(@fmincon,'Algorithm','sqp');
        
        % OriTun model with all the parameters -- least square fit
        lsf_OriTun  = @(x) fitVonMises(x, orientationForTrials, cur_neuron, 'OriTun');
        problem_OriTun = createOptimProblem('fmincon','objective',... 
        lsf_OriTun,'x0',x0_OriTun, 'lb', lb1_OriTun, 'ub', ub1_OriTun, 'options',opts);
        ms_OriTun = MultiStart('Display', 'off');
        [x_OriTun(ii, :),f_OriTun(ii)] = run(ms_OriTun,problem_OriTun, MultStartNum);
        
        % DirTun model with all the parameters -- least square fit
        lsf_DirTun  = @(x) fitVonMises(x, orientationForTrials, cur_neuron, 'DirTun');
        problem_DirTun = createOptimProblem('fmincon','objective',... 
        lsf_DirTun,'x0',x0_DirTun, 'lb', lb1_DirTun, 'ub', ub1_DirTun, 'options',opts);
        ms_DirTun = MultiStart('Display', 'off');
        [x_DirTun(ii, :),f_DirTun(ii)] = run(ms_DirTun,problem_DirTun, MultStartNum);

        % NULL model with constant mean and variance -- least square fit
        lsf_null  = @(x) fitVonMises(x, orientationForTrials, cur_neuron, 'null');
        problem_null = createOptimProblem('fmincon','objective',... 
        lsf_null,'x0',x0_null, 'lb', 0, 'ub', inf, 'options',opts);
        ms_null = MultiStart('Display', 'off');
        [x_null(ii, :),f_null(ii)] = run(ms_null,problem_null, MultStartNum);


        fprintf('fitting for neuron %d and contrast = %0.2f with %d random start by fmincon \n', ii, ucons(coni), MultStartNum) 

    end
    x_OriTun(:, 4) = mod(x_OriTun(:, 4), 360);   % preferred orientation in the range of [0, 2*pi]
    x_DirTun(:, 4) = mod(x_DirTun(:, 4), 360);   % preferred orientation in the range of [0, 2*pi]

    for nn = 1: nNeurons
        % residuals for max likelihood meathod        
        
        SSE_OriTun(nn) = sum((OriTun_vonMises(x_OriTun(nn, :), orientationForTrials)- averageSpkRateOverTrial(:, nn)).^2);
        SSE_DirTun(nn) = sum((DirTun_vonMises(x_DirTun(nn, :), orientationForTrials)- averageSpkRateOverTrial(:, nn)).^2);
        SSE_null(nn) = sum((NullFunc(x_null(nn, :), orientationForTrials)- averageSpkRateOverTrial(:, nn)).^2);
        
        R2_OriTun(nn)   = corr(OriTun_vonMises(x_OriTun(nn, :), orientationForTrials),  averageSpkRateOverTrial(:, nn))^2;
        R2_DirTun(nn)   = corr(DirTun_vonMises(x_DirTun(nn, :), orientationForTrials), averageSpkRateOverTrial(:, nn))^2;
        R2_null(nn)     = corr(NullFunc(x_null(nn, :), orientationForTrials), averageSpkRateOverTrial(:, nn))^2;
    end
    
    
    %%  test orientation tuning vs direction tuning
    v1 = nParam_DirTun- nParam_OriTun;
    v2 = nTrial- nParam_DirTun;
    f_ratio_OriVsDir = v2*(SSE_OriTun-SSE_DirTun)./ (v1* SSE_DirTun);
    alpha_OriVsDir  = 1- fcdf(f_ratio_OriVsDir, v1, v2);

    
    %%  test orientation tuning vs non-tuning (null) model
    v1 = nParam_OriTun- nParam_null;
    v2 = nTrial- nParam_OriTun;
    f_ratio_OriVsNull = v2*(SSE_null-SSE_OriTun)./ (v1* SSE_OriTun);
    alpha_OriVsNull  = 1- fcdf(f_ratio_OriVsNull, v1, v2);

    nDirNeur05 = sum(alpha_OriVsDir <pval05/2 & alpha_OriVsNull<pval05/2);
    nOriNeur05 = sum(alpha_OriVsDir >pval05/2 & alpha_OriVsNull<pval05/2);    
    
    nDirNeur01 = sum(alpha_OriVsDir <pval01/2 & alpha_OriVsNull<pval01/2);
    nOriNeur01 = sum(alpha_OriVsDir >pval01/2 & alpha_OriVsNull<pval01/2);    
    
    nDirNeur001 = sum(alpha_OriVsDir <pval001/2 & alpha_OriVsNull<pval001/2);
    nOriNeur001 = sum(alpha_OriVsDir >pval001/2 & alpha_OriVsNull<pval001/2);
    
    
    figure('Position', [50, 50, 700, 900]); hold on;
    subplot(2, 1, 1); plot(alpha_OriVsDir, '*r', 'MarkerSize', 4);
    ylabel('\alpha'); xlabel('neuron number'); xlim([1, nNeurons]);
    title(sprintf('F test, OriTun vs DirTun Bonferoni-corrected (0.05:%d, 0.01:%d, 0.001:%d)', ...
        nDirNeur05,nDirNeur01,nDirNeur001));
    
    subplot(2, 1, 2); plot(alpha_OriVsNull, '*r', 'MarkerSize', 4);
    ylabel('\alpha'); xlabel('neuron number'); xlim([1, nNeurons]);
    % applying bonferoni correction due to multiple comparsions
    title(sprintf('F test, least square, # of tunned neurons Bonferoni-corrected (0.05:%d, 0.01:%d, 0.001:%d)',...
        nOriNeur05,nOriNeur01,nOriNeur001));
        
    %% plot of fitted tuning curves for some of neurons chosen randomly   
    figure('Position', [50, 50, 1300, 800]); 
    for plt = 1: length(rnd_plt)
        subplot(4,4, plt); hold on;  
        plot(Theta_plt, DirTun_vonMises(x_DirTun(rnd_plt(plt), :), Theta_plt), 'LineWidth', 2, 'Color', [0.75, 0.75, 0]);
        plot(Theta_plt, OriTun_vonMises(x_OriTun(rnd_plt(plt), :), Theta_plt), 'LineWidth', 2, 'Color', 'r'); 
        plot(Theta_plt, NullFunc(x_null(rnd_plt(plt), :), Theta_plt), ...
            'LineWidth', 2, 'Color', 'k'); hold on;   
        meanOverOri = zeros(size(orientations));
        for ori = 1: length(orientations)
            meanOverOri(ori, 1) = mean(averageSpkRateOverTrial(orientationForTrials== orientations(ori), rnd_plt(plt)));
        end     
        [uoris_plt, sr_ori_ave, sr_ori_neg, sr_ori_pos] = errbar_util(...
            orientationForTrials, averageSpkRateOverTrial(:, rnd_plt(plt)));         
        errorbar(uoris_plt, sr_ori_ave, sr_ori_ave-sr_ori_neg, sr_ori_pos-sr_ori_ave, 'LineStyle', 'none', 'Color', 'b', 'LineWidth', 2); 
        
        scatter(uoris_plt, sr_ori_ave, 30,[0.3010, 0.7450, 0.9330],'filled') 
        ylim([0 max(averageSpkRateOverTrial(:, rnd_plt(plt)))])
        xticks(orientations);
        
        if plt == 1
            legend('DirTun von-Mises fit','OriTun von-Mises fit', 'null model', 'data (90% interval)', 'data (means)')
        end
        ylim([0, inf]);
        title(sprintf('neur # %d (%0.3f, %0.3f)', rnd_plt(plt), ...
            alpha_OriVsDir(rnd_plt(plt)), alpha_OriVsNull(rnd_plt(plt))));
    end    
    suplabel('orientation', 'x'); suplabel('rate', 'y');
    suplabel('direction tuned neurons (\alpha_{Ori Vs Dir}, \alpha_{Ori Vs Null})', 't')


    coef.x_OriTun{coni,1} = x_OriTun;
    coef.x_DirTun{coni,1} = x_DirTun;
    coef.x_null{coni,1} = x_null;
    
    coef.f_OriTun{coni,1} = f_OriTun;
    coef.f_DirTun{coni,1} = f_DirTun;
    coef.f_null{coni,1} = f_null;

    coef.DirTun_vonMises{coni,1} = DirTun_vonMises;
    coef.OriTun_vonMises{coni,1} = OriTun_vonMises;
    coef.NullFunc{coni,1} = NullFunc;
    
    coef.f_ratio_OriVsDir{coni,1} = f_ratio_OriVsDir;
    coef.f_ratio_OriVsNull{coni,1} = f_ratio_OriVsNull;
    
    coef.alpha_OriVsDir{coni,1} = alpha_OriVsDir;
    coef.alpha_OriVsNull{coni,1} = alpha_OriVsNull;
    
    coef.R2_OriTun  = R2_OriTun;
    coef.R2_DirTun  = R2_DirTun;
    coef.R2_null    = R2_null;

end 
%% write data to file
fprintf('Writing TC fitted coeficients to %s\n', TC_coef_path);
save(TC_coef_path, 'coef');
T_total = toc;
fprintf('total running time = %0.2f hours \n', T_total/3600);

close all;

%% direct fitting of the means ignoring the variances
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


function sr_avg = SS_total_func(ori, sr)
    uoris           = sort(unique(ori), 'ascend');
    orin            = length(uoris);
    sr_avg          = nan(length(sr), 1);
    
    for ii = 1: orin
        sr_ori_tmp  = sr(ori == uoris(ii));
        sr_avg(ori == uoris(ii))  = mean(sr_ori_tmp);
    end

    
function [uoris, sr_ori_ave, sr_ori_neg, sr_ori_pos] = errbar_util(ori, sr)
    uoris           = sort(unique(ori), 'ascend');
    orin            = length(uoris);
    sr_ori_ave      = nan(orin, 1);
    sr_ori_neg      = nan(orin, 1);
    sr_ori_pos      = nan(orin, 1);    
    for ii = 1: orin
        sr_ori_tmp  = sr(ori == uoris(ii));
        sr_ori_ave(ii)  = mean(sr_ori_tmp);
        sr_ori_neg(ii)  = prctile(sr_ori_tmp, 5);
        sr_ori_pos(ii)  = prctile(sr_ori_tmp, 95);
    end
