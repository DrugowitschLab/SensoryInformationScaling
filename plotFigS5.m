function plotFigS5
%% plots detailed analysis of which factors contribute to speed-dep info boost

%% (general) settings
singlefig = true;   % set to true to have all planels in single figure
momentpath = ['.' filesep 'moment_cache'];
fitpath = ['.' filesep 'fits'];
locol = [31 120 180] / 255;
hicol = [227 26 28] / 255;
plotprctiles = [5 25 50 75 95];
afrac = 0.95;

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
if singlefig
    f = figure('Name', 'Figure', 'Units', 'centimeters', 'Position', [0 0 14.5 4.8]);
    f.Renderer = 'Painter'; % ensure vector output
end

%% info comparison for low vs high speed across sessions/mice

datasets = {'m25a', 'm25b', 'm26a', 'm26b', ...
    'aj42a', 'aj42b', 'aj42c', 'aj42d', 'aj42e', ...
    'aj43a', 'aj43b', 'aj43c', 'aj43d', 'aj43e', 'aj43f', 'aj43g'};
d45discrs = {'o1-2','o2-3','o3-4','o4-5','o5-6','o6-7','o7-8','o1-8'};
dsmouse = [1 1 2 2 3 3 3 3 3 4 4 4 4 4 4 4];
dssession = [1 2 1 2 1 2 3 4 5 1 2 3 4 5 6 7];
dslabels = {'1A','1B', '2A', '2B', '3A','3B','3C','3D', '3E', '4A','4B','4C','4D','4E', '4F', '4G'};

ndataset = length(datasets);

%%
if exist('plotFigS5.mat')
    load('plotFigS5.mat')
else
    for imouse = 1: length(datasets)
        for iori = 1: length(d45discrs)
            momentfile = [datasets{imouse} '_' d45discrs{iori} '_c1'];

            fprintf('Loading %s ...\n', [momentpath filesep momentfile '_hispd.mat']);
            dhi = load([momentpath filesep momentfile '_hispd.mat']);

            fprintf('Loading %s ...\n', [momentpath filesep momentfile '_losdp.mat']);
            dlo = load([momentpath filesep momentfile '_lospd.mat']);   
        
            
            % varhi(imouse, iori) = sum(diag(dhi.S));
            % varlo(imouse, iori) = sum(diag(dlo.S));
            

            [Ihi(imouse, iori), Ihi_sd(imouse, iori)] = empTotalInf(dhi);
            [Ilo(imouse, iori), Ilo_sd(imouse, iori)] = empTotalInf(dlo);      
            
            pval(imouse, iori)    = normcdf(0, Ihi(imouse, iori)- Ilo(imouse, iori), ...
                sqrt(Ihi_sd(imouse, iori)^2+ Ilo_sd(imouse, iori)^2));


            [Imu_tmp, Isd_tmp]              = populationInfo(dhi.mu, dlo.S, dhi.T, dlo.T, dhi.ds); 
            Ihi_fhi_siglo(imouse, iori)     = Imu_tmp(end);
            Ihi_fhi_siglo_sd(imouse, iori)  = Isd_tmp(end);
            pval_fhi_siglo(imouse, iori)    = normcdf(0, Ihi_fhi_siglo(imouse, iori)- Ilo(imouse, iori), ...
                sqrt(Ihi_fhi_siglo_sd(imouse, iori)^2+ Ilo_sd(imouse, iori)^2));


            [Imu_tmp, Isd_tmp]              = populationInfo(dlo.mu, dhi.S, dlo.T, dhi.T, dhi.ds); 
            Ihi_flo_sighi(imouse, iori)     = Imu_tmp(end);
            Ihi_flo_sighi_sd(imouse, iori)  = Isd_tmp(end);
            pval_flo_sighi(imouse, iori)    = normcdf(0, Ihi_flo_sighi(imouse, iori)- Ilo(imouse, iori), ...
                sqrt(Ihi_flo_sighi_sd(imouse, iori)^2+ Ilo_sd(imouse, iori)^2));

            corrlo                          = corrcov(dlo.S); 
            varlo                           = diag(dlo.S);  
            corrhi                          = corrcov(dhi.S); 
            varhi                           = diag(dhi.S); 
           
            S_corrlo                        = corr2cov(sqrt(varhi), corrlo);
            S_corrlo(isnan(S_corrlo))       = 0; % remove nan instances
            
            S_varlo                         = corr2cov(sqrt(varlo), corrhi);
            S_varlo(isnan(S_varlo))         = 0; % remove nan instances

            [Imu_tmp, Isd_tmp]              = populationInfo(dhi.mu, S_corrlo, dhi.T, dhi.T, dhi.ds); 
            Ihi_fhi_corrlo(imouse, iori)    = Imu_tmp(end);
            Ihi_fhi_corrlo_sd(imouse, iori) = Isd_tmp(end);
            pval_fhi_corrlo(imouse, iori)    = normcdf(0, Ihi_fhi_corrlo(imouse, iori)- Ilo(imouse, iori), ...
                sqrt(Ihi_fhi_corrlo_sd(imouse, iori)^2+ Ilo_sd(imouse, iori)^2));

            [Imu_tmp, Isd_tmp]              = populationInfo(dhi.mu, S_varlo, dhi.T, dhi.T, dhi.ds); 
            Ihi_fhi_varlo(imouse, iori)     = Imu_tmp(end);
            Ihi_fhi_varlo_sd(imouse, iori)  = Isd_tmp(end);
            pval_fhi_varlo(imouse, iori)    = normcdf(0, Ihi_fhi_varlo(imouse, iori)- Ilo(imouse, iori), ...
                sqrt(Ihi_fhi_varlo_sd(imouse, iori)^2+ Ilo_sd(imouse, iori)^2));
            
            
            [Imu_tmp, Isd_tmp]              = populationInfo(dlo.mu, S_varlo, dhi.T, dlo.T, dlo.ds); 
            Ihi_flo_varlo(imouse, iori)    = Imu_tmp(end);
            Ihi_flo_varlo_sd(imouse, iori)  = Isd_tmp(end);
            pval_flo_varlo(imouse, iori)    = normcdf(0, Ihi_flo_varlo(imouse, iori)- Ilo(imouse, iori), ...
                sqrt(Ihi_flo_varlo_sd(imouse, iori)^2+ Ilo_sd(imouse, iori)^2));
            
            [Imu_tmp, Isd_tmp]              = populationInfo(dlo.mu, S_corrlo, dhi.T, dlo.T, dhi.ds); 
            Ihi_flo_corrlo(imouse, iori)    = Imu_tmp(end);
            Ihi_flo_corrlo_sd(imouse, iori) = Isd_tmp(end);
            pval_flo_corrlo(imouse, iori)    = normcdf(0, Ihi_flo_corrlo(imouse, iori)- Ilo(imouse, iori), ...
                sqrt(Ihi_flo_corrlo_sd(imouse, iori)^2+ Ilo_sd(imouse, iori)^2));
        end
    end

%%
    save('plotFigS_speed_split_factors_data.mat', 'Ilo', 'Ihi', 'Ilo_sd', 'pval', 'Ihi_fhi_siglo', 'Ihi_fhi_siglo_sd', ...
            'pval_fhi_siglo', 'Ihi_flo_sighi', 'Ihi_flo_sighi_sd', 'pval_flo_sighi', ...
            'Ihi_fhi_corrlo', 'Ihi_fhi_corrlo_sd', 'pval_fhi_corrlo', ...
            'Ihi_fhi_varlo', 'Ihi_fhi_varlo_sd', 'pval_fhi_varlo', ...
            'Ihi_flo_varlo', 'Ihi_flo_varlo_sd', 'pval_flo_varlo',...
            'Ihi_flo_corrlo', 'Ihi_flo_corrlo_sd', 'pval_flo_corrlo', ...
            'datasets', 'd45discrs');
end
    
sgnfcns_level = 0.05; 

% indices for sessions in which Ihi is significantly higher than Ilo
nonsgnfcns_idx = pval > 0.05;
Ilo(nonsgnfcns_idx) = nan;
Ihi(nonsgnfcns_idx) = nan;
Ilo_sd(nonsgnfcns_idx) = nan;
Ihi_flo_sighi(nonsgnfcns_idx) = nan;
Ihi_fhi_siglo(nonsgnfcns_idx) = nan;
Ihi_fhi_corrlo(nonsgnfcns_idx) = nan;
Ihi_fhi_varlo(nonsgnfcns_idx) = nan;
Ihi_flo_corrlo(nonsgnfcns_idx) = nan;
Ihi_flo_varlo(nonsgnfcns_idx) = nan;


%% change in f'
if singlefig
    subplotcm([1 0.8 1.5 3]);  hold on;
    text(-1,3,'a','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
    plot([0 90], [0 90], '-', 'Color', [1 1 1]*0.5);
else
    figure;  hold on;
    plot([0 120], [0 120], '-', 'Color', [1 1 1]*0.5);
end
for i = 1:length(datasets)* length(d45discrs)
    plot(Ilo(i)+[-1 1]*Ilo_sd(i), [1 1]*Ihi_fhi_siglo(i), '-', ...
        'Color', [1 1 1]*0.6, 'LineWidth', 0.25);
    plot([1 1]*Ilo(i), Ihi_fhi_siglo(i)+[-1 1]*Ihi_fhi_siglo_sd(i), '-', ...
        'Color', [1 1 1]*0.6, 'LineWidth', 0.25);
end
% then centers
i = pval_fhi_siglo(:) >= sgnfcns_level;
plot(Ilo(i), Ihi_fhi_siglo(i), 'o', 'MarkerSize', 3, ...
    'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [1 1 1]*0.2);
plot(Ilo(~i), Ihi_fhi_siglo(~i), 'o', 'MarkerSize', 3, ...
    'MarkerFaceColor', [1 1 1]*0.2, 'MarkerEdgeColor', 'none');
if singlefig, xlim([0 90]); ylim([0 180]);
else, xlim([0 120]);  ylim([0 180]); end
set(gca,'Box','off','XTick',0:50:100,'YTick',0:50:150);
xlabel('I_{lo}');
ylabel('Fisher information, f_{hi}');
title('change of f''')

[~,p,~,stats] = ttest(Ilo(:), Ihi_fhi_siglo(:));

fprintf(['stats for testing Ilo vs Ihi (flo) %d dataset with total %d '...
    'independent tasks;\n p-val = %0.2e \n tstats = %0.2f \n df = %d\n '...
    'sd = %0.2f\n'],...
    ndataset, length(Ilo(:)), p, stats.tstat, stats.df, stats.sd);



%% change in sigma
if singlefig
    subplotcm([3 0.8 1.5 3]);  hold on;
    plot([0 90], [0 90], '-', 'Color', [1 1 1]*0.5);
else
    figure;  hold on;
    plot([0 120], [0 120], '-', 'Color', [1 1 1]*0.5);
end
for i = 1:length(datasets)* length(d45discrs)
    plot(Ilo(i)+[-1 1]*Ilo_sd(i), [1 1]*Ihi_flo_sighi(i), '-', ...
        'Color', [1 1 1]*0.6, 'LineWidth', 0.25);
    plot([1 1]*Ilo(i), Ihi_flo_sighi(i)+[-1 1]*Ihi_flo_sighi_sd(i), '-', ...
        'Color', [1 1 1]*0.6, 'LineWidth', 0.25);
end
% then centers
i = pval_flo_sighi(:) >= sgnfcns_level;
plot(Ilo(i), Ihi_flo_sighi(i), 'o', 'MarkerSize', 3, ...
    'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [1 1 1]*0.2);
plot(Ilo(~i), Ihi_flo_sighi(~i), 'o', 'MarkerSize', 3, ...
    'MarkerFaceColor', [1 1 1]*0.2, 'MarkerEdgeColor', 'none');
if singlefig
    xlim([0 90]); ylim([0 180]);
    set(gca,'YColor','none');
else
    xlim([0 120]);  ylim([0 180]);
    xlabel('I_{lo}');
    ylabel('Fisher information, \Sigma_{hi}');
end
set(gca,'Box','off','XTick',0:50:100,'YTick',0:50:150);
title('change of covariance');

[~,p,~,stats] = ttest(Ilo(:), Ihi_flo_sighi(:));

fprintf(['stats for testing Ilo vs Ihi (siglo) %d dataset with total %d '...
    'independent tasks;\n p-val = %0.2e \n tstats = %0.2f \n df = %d\n '...
    'sd = %0.2f\n'],...
    ndataset, length(Ilo(:)), p, stats.tstat, stats.df, stats.sd);


%% change in variance with flo
if singlefig
    subplotcm([5 0.8 1.5 3]);  hold on;
    plot([0 90], [0 90], '-', 'Color', [1 1 1]*0.5);
else
    figure;  hold on;
    plot([0 120], [0 120], '-', 'Color', [1 1 1]*0.5);
end
for i = 1:length(datasets)* length(d45discrs)
    plot(Ilo(i)+[-1 1]*Ilo_sd(i), [1 1]*Ihi_flo_corrlo(i), '-', ...
        'Color', [1 1 1]*0.6, 'LineWidth', 0.25);
    plot([1 1]*Ilo(i), Ihi_flo_corrlo(i)+[-1 1]*Ihi_flo_corrlo_sd(i), '-', ...
        'Color', [1 1 1]*0.6, 'LineWidth', 0.25);
end
% then centers
i = pval_fhi_corrlo(:) >= sgnfcns_level;
plot(Ilo(i), Ihi_flo_corrlo(i), 'o', 'MarkerSize', 3, ...
    'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [1 1 1]*0.2);
plot(Ilo(~i), Ihi_flo_corrlo(~i), 'o', 'MarkerSize', 3, ...
    'MarkerFaceColor', [1 1 1]*0.2, 'MarkerEdgeColor', 'none');
if singlefig
    xlim([0 90]); ylim([0 180]);
    set(gca,'YColor','none');    
else
    xlim([0 120]);  ylim([0 180]);
    xlabel('I_{lo}');
    ylabel('Fisher information, \sigma^2_{hi}');    
end
set(gca,'Box','off','XTick',0:50:100,'YTick',0:50:150);
title('change of variance, flo')

[~,p,~,stats] = ttest(Ilo(:), Ihi_flo_corrlo(:));

fprintf(['stats for testing Ilo vs Ihi (varlo, flo) %d dataset with total %d '...
    'independent tasks;\n p-val = %0.2e \n tstats = %0.2f \n df = %d\n '...
    'sd = %0.2f\n'],...
    ndataset, length(Ilo(:)), p, stats.tstat, stats.df, stats.sd);


%% change in correlation with flo
if singlefig
    subplotcm([7 0.8 1.5 3]);  hold on;
    plot([0 90], [0 90], '-', 'Color', [1 1 1]*0.5);
else
    figure;  hold on;
    plot([0 120], [0 120], '-', 'Color', [1 1 1]*0.5);
end
plot([0 120], [0 120], '-', 'Color', [1 1 1]*0.5);
for i = 1:length(datasets)* length(d45discrs)
    plot(Ilo(i)+[-1 1]*Ilo_sd(i), [1 1]*Ihi_flo_varlo(i), '-', ...
        'Color', [1 1 1]*0.6, 'LineWidth', 0.25);
    plot([1 1]*Ilo(i), Ihi_flo_varlo(i)+[-1 1]*Ihi_flo_varlo_sd(i), '-', ...
        'Color', [1 1 1]*0.6, 'LineWidth', 0.25);
end
% then centers
i = pval_fhi_varlo(:) >= sgnfcns_level;
plot(Ilo(i), Ihi_flo_varlo(i), 'o', 'MarkerSize', 3, ...
    'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [1 1 1]*0.2);
plot(Ilo(~i), Ihi_flo_varlo(~i), 'o', 'MarkerSize', 3, ...
    'MarkerFaceColor', [1 1 1]*0.2, 'MarkerEdgeColor', 'none');
if singlefig
    xlim([0 90]); ylim([0 180]);
    set(gca,'YColor','none');    
else
    xlim([0 120]);  ylim([0 180]);
    xlabel('I_{lo}');
    ylabel('Fisher information, \sigma^2_{hi}');    
end
set(gca,'Box','off','XTick',0:50:100,'YTick',0:50:150);
title('change of correlation with flo')

[~,p,~,stats] = ttest(Ilo(:), Ihi_flo_varlo(:));

fprintf(['stats for testing Ilo vs Ihi (corrlo, flo) %d dataset with total %d '...
    'independent tasks;\n p-val = %0.2e \n tstats = %0.2f \n df = %d\n '...
    'sd = %0.2f\n'],...
    ndataset, length(Ilo(:)), p, stats.tstat, stats.df, stats.sd);



%% summary impact analyses
Ihi_flo_sighi_chng = 100* (Ihi_flo_sighi- Ilo)./(Ihi- Ilo);
Ihi_fhi_corrl_chng = 100* (Ihi_fhi_corrlo- Ilo)./(Ihi- Ilo);
Ihi_fhi_siglo_chng = 100* (Ihi_fhi_siglo- Ilo)./(Ihi- Ilo);
Ihi_flo_corrl_chng = 100* (Ihi_flo_corrlo- Ilo)./(Ihi- Ilo);
Ihi_flo_varlo_chng = 100* (Ihi_flo_varlo- Ilo)./(Ihi- Ilo);
Ihi_fhi_varlo_chng = 100* (Ihi_fhi_varlo- Ilo)./(Ihi- Ilo);

data = [Ihi_fhi_siglo_chng(:) Ihi_flo_sighi_chng(:) Ihi_flo_corrl_chng(:) ...
    Ihi_fhi_corrl_chng(:) Ihi_flo_varlo_chng(:) Ihi_fhi_varlo_chng(:)  ];

xtk = 1:6;

      
if singlefig
    subplotcm([9.5 0.8 4 3]);  hold on;
    text(-1,3,'b','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
    colors = [27 158 119; 217 95 2; 117 112 179; ...
              231 41 138; 102 166 30; 230 171 2] / 255;
else
    figure;  hold on;
    colors = [0.6520    0.2754    0.8801; ...
              0.0662    0.6818    0.8443; ...
              0.1       0.1       0.8801; ...
              0.8       0.8       0.1; ...
              1         0         1; ...
              0         1         1];
end

for i = xtk
    boxchart(xtk(i)*ones(size(data(:,i))), data(:,i), 'BoxFaceColor', colors(i,:), 'MarkerStyle','none')
end
ylabel('percentage change')
title('impact analysis');
ylim([0, 600])
xlim([0.5 6.5]);
set(gca,'Box','off','XTick',xtk,'XTickLabels',...
    {'f_{hi},\Sigma_{lo}', 'f_{lo},\Sigma_{hi}', 'f_{lo},\sigma_{hi},R_{lo}', ...
     'f_{hi},\sigma_{hi},R_{lo}', 'f_{lo},\sigma_{lo},R_{hi}', 'f_{hi},\sigma_{lo},R_{lo}'});

%% save figure
if singlefig
    fprintf('\nWriting figure to figS5.pdf\n');
    print(['figs' filesep 'figS5'], '-dpdf');
end



function [Imu, Isd] = populationInfo(mu, S, Tf, Ts, ds)
N = length(mu);
ordn = 100;
In = NaN(ordn, N);
for ordi = 1:ordn
    if mod(ordi, 10) == 0, fprintf('%d ', ordi); end
    In(ordi,:) = empInfscaling_spd(mu, S, randperm(N), Tf, Ts, ds);
end
fprintf('\n');

% statistics
Iincr = diff([zeros(ordn,1) In], 1, 2);
Imu = [0 mean(cumsum(Iincr, 2))];
Isd = [0 sqrt(cumsum(var(Iincr, [], 1)))];

function [In, sdIn] = empTotalInf(d)
%% return total information in recorded population in given dataset
In = sum(mean(d.Iincr_samples,1));
varIn = sum(var(d.Iincr_samples,[],1));
sdIn = sqrt(varIn);

