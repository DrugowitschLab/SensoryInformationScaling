function plotFig3
%% plots panels of figure 3


%% (general) settings
singlefig = true;   % set to true to have all planels in single figure
momentpath = ['.' filesep 'moment_cache'];
shufcol = [0.89 0.1 0.11];

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
    figure('Name', 'Figure', 'Units', 'centimeters', 'Position', [0 0 10.9 9.6]);
end
unreswarning = warning('query', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle');


%% info scaling for example discrimiation, measured vs. shuffled + p(correct)
momentfile = 'm25b_o3-4_c1';
% load data and compute stats
fprintf('Loading %s ...\n', [momentpath filesep momentfile '.mat']);
d = load([momentpath filesep momentfile '.mat']);
fprintf('Loading %s ...\n', [momentpath filesep momentfile '_shuf.mat']);
dshuf = load([momentpath filesep momentfile '_shuf.mat']);
Iincr_mu = mean(d.Iincr_samples, 1);
Iincr_var = var(d.Iincr_samples, [], 1);
I_mu = [0 cumsum(Iincr_mu)];
I_var = [0 cumsum(Iincr_var)];
pc_mu = normcdf(0.5*sqrt(I_mu*d.ds^2));
pc_up = normcdf(0.5*sqrt((I_mu+sqrt(I_var))*d.ds^2));
pc_lo = normcdf(0.5*sqrt(max(0,I_mu-sqrt(I_var))*d.ds^2));
Iincr_mu_shuf = mean(dshuf.Iincr_samples, 1);
Iincr_var_shuf = var(dshuf.Iincr_samples, [], 1);
I_mu_shuf = [0 cumsum(Iincr_mu_shuf)];
I_var_shuf = [0 cumsum(Iincr_var_shuf)];
pc_mu_shuf = normcdf(0.5*sqrt(I_mu_shuf*d.ds^2));
pc_up_shuf = normcdf(0.5*sqrt((I_mu_shuf+sqrt(I_var_shuf))*d.ds^2));
pc_lo_shuf = normcdf(0.5*sqrt(max(0,I_mu_shuf-sqrt(I_var_shuf))*d.ds^2));
% plot info scaling
if singlefig
    subplotcm([1 5.3 4.5 3.5]);  hold on;
    text(-1,4,'a','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
alpha(patch([0:(length(I_mu)-1) fliplr(0:(length(I_mu)-1))], ...
        [(I_mu + sqrt(I_var)) fliplr(I_mu - sqrt(I_var))], 1, ...
        'FaceColor',[1 1 1]*0.2,'EdgeColor','none'), 0.2);
alpha(patch([0:(length(I_mu_shuf)-1) fliplr(0:(length(I_mu_shuf)-1))], ...
        [(I_mu_shuf + sqrt(I_var_shuf)) fliplr(I_mu_shuf - sqrt(I_var_shuf))], 1, ...
        'FaceColor',shufcol,'EdgeColor','none'), 0.2);
plot(0:(length(I_mu)-1), I_mu, '-', 'Color', [1 1 1]*0.2, 'LineWidth', 1);
plot(0:(length(I_mu_shuf)-1), I_mu_shuf, '-', 'LineWidth', 1, 'Color', shufcol);
xlim([0 (length(I_mu)-1)]);  ylim([0 45]);
set(gca,'Box','off','XTick',[0 100 200 300],'YTick',[0 10 20 30 40]);
xlabel('number of neurons N');
ylabel('Fisher information [rad^{-2}]');
% plot p(correct) scaling
if singlefig
    subplotcm([1.7 7.4 1.8 1.4], 'centimeters', true);  hold on;
else
    figure;  hold on;
end
alpha(patch([0:(length(pc_mu)-1) fliplr(0:(length(pc_mu)-1))], ...
        [pc_up fliplr(pc_lo)], 1, ...
        'FaceColor',[1 1 1] * 0.2,'EdgeColor','none'), 0.2);
alpha(patch([0:(length(pc_mu_shuf)-1) fliplr(0:(length(pc_mu_shuf)-1))], ...
        [pc_up_shuf fliplr(pc_lo_shuf)], 1, ...
        'FaceColor',shufcol,'EdgeColor','none'), 0.2);
plot(0:(length(pc_mu)-1), pc_mu, '-', 'Color', [1 1 1]*0.2, 'LineWidth', 1);
plot(0:(length(pc_mu_shuf)-1), pc_mu_shuf, '-', 'LineWidth', 1, 'Color', shufcol);
xlim([0 (length(pc_mu)-1)]);  ylim([0.5 1]);
set(gca,'Box','off','XTick',[0 100 200 300],'YTick',[0.5 1]);
xlabel('number of neurons N');
ylabel('p(correct)');


%% mapping between information and orientation discrimination threshold (80% corr)
If2thresh = @(If) 180/pi * sqrt(2) * norminv(0.8) / sqrt(If);
Imesh = logspace(1,4,100);
thresh = nan(size(Imesh));
for ii = 1: length(Imesh)
    thresh(ii) = If2thresh(Imesh(ii));
end
if singlefig
    subplotcm([1 0.8 4.5 3.5]);  hold on;
    text(-1,4,'c','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
plot(Imesh, thresh, '-', 'Color', [1 1 1]*0.2, 'LineWidth', 1);
xlim([10^1 10^4]);  ylim([0 25]);
set(gca,'Box','off','XScale','log','YTick',[0 5 10 15 20 25]);
xlabel('Fisher information [rad^{-2}]'); 
ylabel('direction discrimination threshold [deg]');


%% orientation discrimination thresholds for =45deg across different drift 
% direction pairs, measures vs. shuffled.
dorifiles = {'m25b_o1-8_c1.mat' ...
    'm25b_o1-2_c1.mat' ...
    'm25b_o2-3_c1.mat' ...
    'm25b_o3-4_c1.mat' ...
    'm25b_o4-5_c1.mat' ...
    'm25b_o5-6_c1.mat' ...
    'm25b_o6-7_c1.mat' ...
    'm25b_o7-8_c1.mat'};
dorinames = {'0\circ vs. 45\circ', '45\circ vs. 90\circ', ...
    '90\circ vs. 135\circ', '135\circ vs. 180\circ', ...
    '180\circ vs. 225\circ', '225\circ vs. 270\circ', ...
    '270\circ vs. 315\circ', '315\circ vs. 0\circ'};
if singlefig
    subplotcm([6.9 0.8 3.5 3.5]);  hold on;
    text(-1,4,'d','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
for i = 1:length(dorifiles)
    fprintf('Loading %s ...\n', [momentpath filesep dorifiles{i}]);
    d = load([momentpath filesep dorifiles{i}]);
    [In, varIn] = empTotalInf(d);
    plot([1 1]*i, [If2thresh(In+sqrt(varIn)) If2thresh(In-sqrt(varIn))], ...
        '-', 'LineWidth', 0.25, 'Color', [1 1 1]*0.6);
    plot(i, If2thresh(In), 'o', 'MarkerSize', 3, ...
        'MarkerFaceColor', [1 1 1]*0.2, 'MarkerEdgeColor', 'none');
end
xlim([0.75 length(dorifiles)+0.25]);  ylim([0 20]);
ylabel('direction discrimination threshold [deg]');
set(gca,'Box','off','YTick',[0 5 10 15 20],...
    'XTick',1:length(dorifiles),'XTickLabels',dorinames,'XTickLabelRotation',45);


%% Bootstrap p-value for comparing I vs I_shuf 
% collect data and comp
datasets = {'m25a', 'm25b', 'm26a', 'm26b', ...
        'aj42a', 'aj42b', 'aj42c', 'aj42d', 'aj42e', ...
        'aj43a', 'aj43b', 'aj43c', 'aj43d', 'aj43e', 'aj43f', 'aj43g', ...
%         'aj60a', 'aj60b', 'aj60c', 'aj60d', ...
%         'aj61a', 'aj61b', 'aj61c'
        };
statdiscrs = {'o1-2','o3-4','o5-6','o7-8'};  % to collect for stats

dori_val = 45;
sgnfcns_level = 0.05; 

ndataset = length(datasets);
dataid = [];
In_all = [];
var_all = [];
In_shuf_all = [];
var_shuf_all = [];
pval_all = [];
In_stats = [];
var_stats = [];
In_shuf_stats = [];
fprintf('Orientation difference %5.1f\n', dori_val);
for dataseti = 1: ndataset
    dataset = datasets{dataseti};
    di = dataInfo(datasets{dataseti});
    oricomb = di.oricomb(1:2, di.oricomb(3,:) == dori_val);
    oricombn = size(oricomb, 2);
    for coni = 1: length(di.cons)
        for oricompi = 1: oricombn
            prefix = sprintf('_o%d-%d_c%d', oricomb(1, oricompi), ...
                oricomb(2, oricompi), coni);
            momentfile = [dataset prefix];
            d = load([momentpath filesep momentfile '.mat']);
            dshuf = load([momentpath filesep momentfile '_shuf.mat']);
            [In, varIn] = empTotalInf(d);
            [In_shuf, varIn_shuf] = empTotalInf(dshuf);
            pval = normcdf(0, In_shuf- In, sqrt(varIn_shuf+ varIn));
            fprintf('p-value comparing I vs I_shuf for %s is %.5e \n', ...
                momentfile, pval)
            In_all = cat(2, In_all, In);
            var_all = cat(2, var_all, varIn);
            In_shuf_all = cat(2, In_shuf_all, In_shuf);
            var_shuf_all = cat(2, var_shuf_all, varIn_shuf);
            pval_all = cat(2, pval_all, pval);
            if strcmp(momentfile, 'm25b_o3-4_c1')
               In_o3_4 = In;
               In_shuf_o3_4 = In_shuf;
               varIn_o3_4 = varIn;
               varIn_shuf_o3_4 = varIn_shuf;
            end
            if any(strcmp(prefix(2:5), statdiscrs))
                In_stats = cat(2, In_stats, In);
                In_shuf_stats = cat(2, In_shuf_stats, In_shuf);
                var_stats = cat(2, var_stats, varIn);
                dataid = cat(2, dataid, dataseti);
            end
        end
    end
end
% plot results
if singlefig
    subplotcm([6.9 5.3 3.5 3.5]);  hold on;
    text(-1,4,'b','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
plot([0 120], [0 120], '-', 'Color', [1 1 1]*0.5);
% error bars first
for i = 1:length(In_all)
    plot(In_all(i)+[-1 1]*sqrt(var_all(i)), [1 1]*In_shuf_all(i), '-', ...
        'Color', [1 1 1]*0.6, 'LineWidth', 0.25);
    plot([1 1]*In_all(i), In_shuf_all(i)+[-1 1]*sqrt(var_shuf_all(i)), '-', ...
        'Color', [1 1 1]*0.6, 'LineWidth', 0.25);
end
% then centers
i = pval_all >= sgnfcns_level;
plot(In_all(i), In_shuf_all(i), 'o', 'MarkerSize', 3, ...
    'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [1 1 1]*0.2);
plot(In_all(~i), In_shuf_all(~i), 'o', 'MarkerSize', 3, ...
    'MarkerFaceColor', [1 1 1]*0.2, 'MarkerEdgeColor', 'none');
xlim([0 120]);  ylim([0 120]);
set(gca,'Box','off','XTick',0:20:140,'YTick',0:20:140);
xlabel('Fisher information [rad^{-2}]');
ylabel('Fisher information, shuffled [rad^{-2}]');


%% stats for comparing I vs I_shuf for 135 vs 180 discrimination in m25b dataset
pval_o3_4 = normcdf(0, In_shuf_o3_4-In_o3_4, sqrt(varIn_o3_4+varIn_shuf_o3_4));
fprintf(['stats for comparing I vs I_shuf for 135 vs 180 discrimination '...
    'task \n I-->(ave=%0.2f, var=%0.2f) \n I_shuf--> '...
    '(ave=%0.2f, var=%0.2f)\n pval=%0.2e\n'], In_o3_4, varIn_o3_4, ...
    In_shuf_o3_4, varIn_shuf_o3_4, pval_o3_4);

fprintf(['discrimination threshold for 135 vs 180 discrimination '...
    'task \n is equal to %0.2f\n'], If2thresh(In_o3_4));

%%

[~,p,~,stats] = ttest(In_stats', In_shuf_stats');

fprintf(['stats for testing I vs I_shuf across %d dataset with total %d '...
    'independent tasks;\n p-val = %0.2e \n tstats = %0.2f \n df = %d\n '...
    'sd = %0.2f\n'],...
    ndataset, length(In_stats), p, stats.tstat, stats.df, stats.sd);

fprintf('Testing if In differs across %d different discriminations:\n', ...
    length(statdiscrs));
for i = 1:length(datasets)
    j = dataid == i;
    stats = mdbstrp_stat(In_stats(j), var_stats(j)); 
    fprintf('%s: p=%f\n', datasets{i}, stats.pvalg);
end


%% save figure
warning(unreswarning.state, 'MATLAB:dispatcher:UnresolvedFunctionHandle');
if singlefig
    fprintf('\nWriting figure to fig3.pdf\n');
    print(['figs' filesep 'fig3'], '-dpdf');
end


%% utility functions

function [In, varIn] = empTotalInf(d)
%% return total information in recorded population in given dataset
In = sum(mean(d.Iincr_samples,1));
varIn = sum(var(d.Iincr_samples,[],1));

function out = rmv_shared_dataset(input)
% removing the discriminations that share dataset to perform independent
% statistical test
cur_datasets = reshape(input(:, 1), 1, []);
cur_col = 1; 

for ii = 2: size(input, 2)
    tmp = reshape(input(:, ii), 1, []);
    if isempty(intersect(tmp, cur_datasets))
        cur_datasets = [cur_datasets, tmp];
        cur_col = [cur_col, ii];
    end
end
out = input(:, cur_col);


function stats = mdbstrp_stat(In, varIn)
MC = 1e5; 
I = mvnrnd(In, diag(varIn), MC);
Iave_all = mean(In); 
I_null = mvnrnd(Iave_all* ones(size(In)), diag(varIn), MC);
Tstat = sum((I- mean(I, 1)).^2, 2);
Tstat_null = sum((I_null- mean(I_null, 1)).^2, 2);
stats.pvalg = sum(Tstat > Tstat_null)/ MC;
