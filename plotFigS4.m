function plotFigS4
%% plots panels of supplementary Fig. showing running speed median split

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
    figure('Name', 'Figure', 'Units', 'centimeters', 'Position', [0 0 13.333 11.5]);
end

%% info scaling for example discrimiation, slow vs. fast
momentfile = 'm25b_o3-4_c1';
% load data and compute stats
fprintf('Loading %s ...\n', [momentpath filesep momentfile '.mat']);
d = load([momentpath filesep momentfile '.mat']);
fprintf('Loading %s ...\n', [momentpath filesep momentfile '_shuf.mat']);
dshuf = load([momentpath filesep momentfile '_shuf.mat']);
fprintf('Loading %s ...\n', [momentpath filesep momentfile '_hispd.mat']);
dhi = load([momentpath filesep momentfile '_hispd.mat']);
fprintf('Loading %s ...\n', [momentpath filesep momentfile '_hishuf.mat']);
dhishuf = load([momentpath filesep momentfile '_hishuf.mat']);
fprintf('Loading %s ...\n', [momentpath filesep momentfile '_losdp.mat']);
dlo = load([momentpath filesep momentfile '_lospd.mat']);
fprintf('Loading %s ...\n', [momentpath filesep momentfile '_loshuf.mat']);
dloshuf = load([momentpath filesep momentfile '_loshuf.mat']);
I_mu = [0 cumsum(mean(d.Iincr_samples, 1))];
I_mushuf = [0 cumsum(mean(dshuf.Iincr_samples, 1))];
I_mu_hi = [0 cumsum(mean(dhi.Iincr_samples, 1))];
I_var_hi = [0 cumsum(var(dhi.Iincr_samples, [], 1))];
I_mu_hishuf = [0 cumsum(mean(dhishuf.Iincr_samples, 1))];
I_mu_lo = [0 cumsum(mean(dlo.Iincr_samples, 1))];
I_var_lo = [0 cumsum(var(dlo.Iincr_samples, [], 1))];
I_mu_loshuf = [0 cumsum(mean(dloshuf.Iincr_samples, 1))];
if singlefig
    subplotcm([1 6.4 4 3.5]);  hold on;
    text(-1,4,'a','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
alpha(patch([0:(length(I_mu_hi)-1) fliplr(0:(length(I_mu_hi)-1))], ...
        [(I_mu_hi + sqrt(I_var_hi)) fliplr(I_mu_hi - sqrt(I_var_hi))], 1, ...
        'FaceColor',hicol,'EdgeColor','none'), 0.2);
alpha(patch([0:(length(I_mu_lo)-1) fliplr(0:(length(I_mu_lo)-1))], ...
        [(I_mu_lo + sqrt(I_var_lo)) fliplr(I_mu_lo - sqrt(I_var_lo))], 1, ...
        'FaceColor',locol,'EdgeColor','none'), 0.2);
plot(0:(length(I_mushuf)-1), I_mushuf, '--', 'Color', [1 1 1]*0.2, 'LineWidth', 1);
plot(0:(length(I_mu_loshuf)-1), I_mu_loshuf, '--', 'Color', locol, 'LineWidth', 1);
plot(0:(length(I_mu_hishuf)-1), I_mu_hishuf, '--', 'Color', hicol, 'LineWidth', 1);
plot(0:(length(I_mu)-1), I_mu, '-', 'Color', [1 1 1]*0.2, 'LineWidth', 1);
plot(0:(length(I_mu_hi)-1), I_mu_hi, '-', 'Color', hicol, 'LineWidth', 1);
plot(0:(length(I_mu_lo)-1), I_mu_lo, '-', 'LineWidth', 1, 'Color', locol);
xlim([0 (length(I_mu_hi)-1)]);  ylim([0 60]);
set(gca,'Box','off','XTick',[0 100 200 300],'YTick',[0 10 20 30 40 50 60]);
xlabel('number of neurons N');
ylabel('Fisher information [rad^{-2}]');


%% orientation discrimination thresholds for 45deg across different drift 
If2thresh = @(If) 180/pi * sqrt(2) * norminv(0.8) / sqrt(If);
% direction pairs, measures vs. shuffled.
dorifiles = {'m25b_o1-8_c1', 'm25b_o1-2_c1', 'm25b_o2-3_c1', 'm25b_o3-4_c1' ...
             'm25b_o4-5_c1', 'm25b_o5-6_c1', 'm25b_o6-7_c1', 'm25b_o7-8_c1'};
dorinames = {'0\circ vs. 45\circ', '45\circ vs. 90\circ', ...
    '90\circ vs. 135\circ', '135\circ vs. 180\circ', ...
    '180\circ vs. 225\circ', '225\circ vs. 270\circ', ...
    '270\circ vs. 315\circ', '315\circ vs. 0\circ'};
if singlefig
    subplotcm([6 6.4 3 3.5]);  hold on;
    text(-1,4,'b','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
for i = 1:length(dorifiles)
    fprintf('Loading %s ...\n', [momentpath filesep dorifiles{i} '.mat']);
    d = load([momentpath filesep dorifiles{i} '.mat']);
    fprintf('Loading %s ...\n', [momentpath filesep dorifiles{i} '_lospd.mat']);
    dlo = load([momentpath filesep dorifiles{i} '_lospd.mat']);
    fprintf('Loading %s ...\n', [momentpath filesep dorifiles{i} '_hispd.mat']);
    dhi = load([momentpath filesep dorifiles{i} '_hispd.mat']);
    [In,~] = empTotalInf(d);
    [In_lo, varIn_lo] = empTotalInf(dlo);
    [In_hi, varIn_hi] = empTotalInf(dhi);
    plot([1 1]*i, [If2thresh(In_lo+sqrt(varIn_lo)) If2thresh(In_lo-sqrt(varIn_lo))], ...
        '-', 'LineWidth', 0.25, 'Color', 0.5+0.5*locol);
    plot([1 1]*i, [If2thresh(In_hi+sqrt(varIn_hi)) If2thresh(In_hi-sqrt(varIn_hi))], ...
        '-', 'LineWidth', 0.25, 'Color', 0.5+0.5*hicol);
    plot(i, If2thresh(In), 'o', 'MarkerSize', 3, ...
        'MarkerFaceColor', [1 1 1]*0.2, 'MarkerEdgeColor', 'none');
    plot(i, If2thresh(In_lo), 'o', 'MarkerSize', 3, ...
        'MarkerFaceColor', locol, 'MarkerEdgeColor', 'none');
    plot(i, If2thresh(In_hi), 'o', 'MarkerSize', 3, ...
        'MarkerFaceColor', hicol, 'MarkerEdgeColor', 'none');
end
xlim([0.75 length(dorifiles)+0.25]);  ylim([0 22]);
ylabel('direction discrimination threshold [deg]');
set(gca,'Box','off','YTick',[0 5 10 15 20],...
    'XTick',1:length(dorifiles),'XTickLabels',dorinames,'XTickLabelRotation',45);


%% Bootstrap p-value for comparing I_lo vs I_hi 
% collect data and comp
%datasets = {'m25a', 'm25b'};
datasets = {'m25a', 'm25b', 'm26a', 'm26b', ...
        'aj42a', 'aj42b', 'aj42c', 'aj42d', 'aj42e', ...
        'aj43a', 'aj43b', 'aj43c', 'aj43d', 'aj43e', 'aj43f', 'aj43g', ...
        };
statdiscrs = {'o1-2','o3-4','o5-6','o7-8'};  % to collect for stats

dori_val = 45;
sgnfcns_level = 0.05; 

ndataset = length(datasets);
dataid = [];
In_lo_all = [];
var_lo_all = [];
In_hi_all = [];
var_hi_all = [];
pval_all = [];
In_lo_stats = [];
var_lo_stats = [];
In_hi_stats = [];
var_hi_stats = [];
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
            dlo = load([momentpath filesep momentfile '_lospd.mat']);
            dhi = load([momentpath filesep momentfile '_hispd.mat']);
            [In_lo, varIn_lo] = empTotalInf(dlo);
            [In_hi, varIn_hi] = empTotalInf(dhi);
            pval = normcdf(0, In_hi- In_lo, sqrt(varIn_hi+ varIn_lo));
            fprintf('p-value comparing I_lo vs I_hi for %s is %.5e \n', ...
                momentfile, pval)
            In_lo_all = cat(2, In_lo_all, In_lo);
            var_lo_all = cat(2, var_lo_all, varIn_lo);
            In_hi_all = cat(2, In_hi_all, In_hi);
            var_hi_all = cat(2, var_hi_all, varIn_hi);
            pval_all = cat(2, pval_all, pval);
            if any(strcmp(prefix(2:5), statdiscrs))
                In_lo_stats = cat(2, In_lo_stats, In_lo);
                In_hi_stats = cat(2, In_hi_stats, In_hi);
                var_lo_stats = cat(2, var_lo_stats, varIn_lo);
                var_hi_stats = cat(2, var_hi_stats, varIn_hi);
                dataid = cat(2, dataid, dataseti);
            end
        end
    end
end
% plot results
if singlefig
    subplotcm([10 6.4 2.333 3.5]);  hold on;
    text(-1,4,'c','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
plot([0 120], [0 120], '-', 'Color', [1 1 1]*0.5);
% error bars first
for i = 1:length(In_lo_all)
    plot(In_lo_all(i)+[-1 1]*sqrt(var_lo_all(i)), [1 1]*In_hi_all(i), '-', ...
        'Color', [1 1 1]*0.6, 'LineWidth', 0.25);
    plot([1 1]*In_lo_all(i), In_hi_all(i)+[-1 1]*sqrt(var_hi_all(i)), '-', ...
        'Color', [1 1 1]*0.6, 'LineWidth', 0.25);
end
% then centers
i = pval_all >= sgnfcns_level;
plot(In_lo_all(i), In_hi_all(i), 'o', 'MarkerSize', 3, ...
    'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [1 1 1]*0.2);
plot(In_lo_all(~i), In_hi_all(~i), 'o', 'MarkerSize', 3, ...
    'MarkerFaceColor', [1 1 1]*0.2, 'MarkerEdgeColor', 'none');
xlim([0 120]);  ylim([0 180]);
set(gca,'Box','off','XTick',0:50:100,'YTick',0:50:150);
xlabel('Fisher information [rad^{-2}]');
ylabel('Fisher information, shuffled [rad^{-2}]');

[~,p,~,stats] = ttest(In_lo_stats', In_hi_stats');

fprintf(['stats for testing I_lo vs I_hi across %d dataset with total %d '...
    'independent tasks;\n p-val = %0.2e \n tstats = %0.2f \n df = %d\n '...
    'sd = %0.2f\n'],...
    ndataset, length(In_hi_stats), p, stats.tstat, stats.df, stats.sd);

fprintf('Testing if In_lo differs across %d different discriminations:\n', ...
    length(statdiscrs));
for i = 1:length(datasets)
    j = dataid == i;
    stats = mdbstrp_stat(In_lo_stats(j), var_lo_stats(j)); 
    fprintf('%s: p=%f\n', datasets{i}, stats.pvalg);
end
fprintf('Testing if In_hi differs across %d different discriminations:\n', ...
    length(statdiscrs));
for i = 1:length(datasets)
    j = dataid == i;
    stats = mdbstrp_stat(In_hi_stats(j), var_hi_stats(j)); 
    fprintf('%s: p=%f\n', datasets{i}, stats.pvalg);
end


%% per-animal across-session plots, Iinf
sessfiles = {...
    {{'m25a_dori1a_c1',  'm25a_dori2a_c1',  'm25a_dori3a_c1'}, ...  % m25
     {'m25b_dori1a_c1',  'm25b_dori2a_c1',  'm25b_dori3a_c1'}}, ...
    {{'m26a_dori1a_c1',  'm26a_dori2a_c1',  'm26a_dori3a_c1'}, ...  % m26
     {'m26b_dori1a_c1',  'm26b_dori2a_c1',  'm26b_dori3a_c1'}}, ...
    {{'aj42a_dori1a_c1', 'aj42a_dori2a_c1', 'aj42a_dori3a_c1'}, ... % aj42
     {'aj42b_dori1a_c1', 'aj42b_dori2a_c1', 'aj42b_dori3a_c1'}, ...
     {'aj42c_dori1a_c1', 'aj42c_dori2a_c1', 'aj42c_dori3a_c1'}, ...
     {'aj42d_dori1a_c1', 'aj42d_dori2a_c1', 'aj42d_dori3a_c1'}, ...
     {'aj42e_dori1a_c1', 'aj42e_dori2a_c1', 'aj42e_dori3a_c1'}}, ...
    {{'aj43a_dori1a_c1', 'aj43a_dori2a_c1', 'aj43a_dori3a_c1'}, ... % aj43
     {'aj43b_dori1a_c1', 'aj43b_dori2a_c1', 'aj43b_dori3a_c1'}, ...
     {'aj43c_dori1a_c1', 'aj43c_dori2a_c1', 'aj43c_dori3a_c1'}, ...
     {'aj43d_dori1a_c1', 'aj43d_dori2a_c1', 'aj43d_dori3a_c1'}, ...
     {'aj43e_dori1a_c1', 'aj43e_dori2a_c1', 'aj43e_dori3a_c1'}, ...
     {'aj43f_dori1a_c1', 'aj43f_dori2a_c1', 'aj43f_dori3a_c1'}, ...
     {'aj43g_dori1a_c1', 'aj43g_dori2a_c1', 'aj43g_dori3a_c1'}}};
animals = length(sessfiles);
sessdorinames = {'45\circ','90\circ','135\circ'};

unreswarning = warning('query', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
pw = (11.333 - 0.5*(animals-1))/animals;
ph = 2.2;
totalsessions = sum(cellfun(@(s) length(s), sessfiles));
Iinfmeds = NaN(totalsessions, length(sessdorinames), 2);
N95meds = NaN(totalsessions, length(sessdorinames), 2);
s = 1;
for i = 1:animals
    sessions = length(sessfiles{i});
    sxlo = @(j) 0.04*(j-1)-0.02*(sessions-1)-0.01;
    sxhi = @(j) 0.04*(j-1)-0.02*(sessions-1)+0.01;
    locols = bsxfun(@times, locol, linspace(1,0.4,sessions)') + ...
        bsxfun(@times, [1 1 1], linspace(0,0.6,sessions)');
    hicols = bsxfun(@times, hicol, linspace(1,0.4,sessions)') + ...
        bsxfun(@times, [1 1 1], linspace(0,0.6,sessions)');
    
    % Iinf estimates
    if singlefig
        subplotcm([(1+(i-1)*(pw+0.5)) (1.0+ph) pw ph]);  hold on;
        if i == 1
            text(-1,2.7,'d','Units','centimeters','FontWeight','bold',...
                'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
        end
    else
        figure;  hold on;
    end
    for j = 1:sessions
        meds = NaN(1, length(sessdorinames));
        for k = 1:length(sessdorinames)
            fprintf('Loading %s ...\n', [fitpath filesep sessfiles{i}{j}{k} '_lospd.mat']);
            m = load([fitpath filesep sessfiles{i}{j}{k} '_lospd']);
            sms = m.fits{2};
            assert(strcmp(sms.name, 'limlin'));  % make sure to get the right model
            assert(strcmp(sms.pnames{1}, 'Iinf'));
            ss = reshape(sms.mc.ss(:,1,:),1,[]);
            meds(k) = prctile(ss, plotprctiles(3));
            plotPrctiles(k+sxlo(j), ss, locols(j,:));
        end
        Iinfmeds(s,:,1) = meds;  s = s + 1;
        plot((1:length(sessdorinames))+sxlo(j), meds, '-', 'Color', locols(j,:));
    end
    s = s - sessions;
    for j = 1:sessions
        meds = NaN(1, length(sessdorinames));
        for k = 1:length(sessdorinames)
            fprintf('Loading %s ...\n', [fitpath filesep sessfiles{i}{j}{k} '_hispd.mat']);
            m = load([fitpath filesep sessfiles{i}{j}{k} '_hispd']);
            sms = m.fits{2};
            assert(strcmp(sms.name, 'limlin'));  % make sure to get the right model
            assert(strcmp(sms.pnames{1}, 'Iinf'));
            ss = reshape(sms.mc.ss(:,1,:),1,[]);
            meds(k) = prctile(ss, plotprctiles(3));
            plotPrctiles(k+sxhi(j), ss, hicols(j,:));
        end
        Iinfmeds(s,:,2) = meds;  s = s + 1;
        plot((1:length(sessdorinames))+sxhi(j), meds, '-', 'Color', hicols(j,:));
    end
    s = s - sessions;
    % plot formatting
    set(gca,'Box','off','YScale','log',...
        'XTick',1:length(sessdorinames),'XTickLabel',[],...
        'YTick',[1:9 10:10:90 100:100:1000],'YTickLabel',...
        {'10^0','','','','','','','','',...
         '10^1','','','','','','','','',...
         '10^2','','','','','','','','','10^3'});
    if i == 1
        ylabel('Fisher information [rad^{-2}]');
    end
    xlim([0.75 length(sessdorinames)+0.25]);  ylim([1 5000]);
    
    % N95 estimates
    if singlefig
        subplotcm([(1+(i-1)*(pw+0.5)) 0.8 pw ph]);  hold on;
        if i == 1
            text(-1,2.7,'e','Units','centimeters','FontWeight','bold',...
                'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
        end
    else
        figure;  hold on;
    end
    for j = 1:sessions
        meds = NaN(1, length(sessdorinames));
        for k = 1:length(sessdorinames)
            fprintf('Loading %s ...\n', [fitpath filesep sessfiles{i}{j}{k} '_lospd.mat']);
            m = load([fitpath filesep sessfiles{i}{j}{k} '_lospd']);
            sms = m.fits{2};
            assert(strcmp(sms.name, 'limlin'));  % make sure to get the right model
            assert(strcmp(sms.pnames{1}, 'Iinf'));
            assert(strcmp(sms.pnames{2}, 'c'));
            Iinfss = reshape(sms.mc.ss(:,1,:),1,[]);
            css = reshape(sms.mc.ss(:,2,:),1,[]);
            N95ss = (afrac/(1-afrac)) * Iinfss ./ css;
            meds(k) = prctile(N95ss, plotprctiles(3));
            plotPrctiles(k+sxlo(j), N95ss, locols(j,:));
        end
        N95meds(s,:,1) = meds;  s = s + 1;
        plot((1:length(sessdorinames))+sxlo(j), meds, '-', 'Color', locols(j,:));
    end
    s = s - sessions;
    for j = 1:sessions
        meds = NaN(1, length(sessdorinames));
        for k = 1:length(sessdorinames)
            fprintf('Loading %s ...\n', [fitpath filesep sessfiles{i}{j}{k} '_hispd.mat']);
            m = load([fitpath filesep sessfiles{i}{j}{k} '_hispd']);
            sms = m.fits{2};
            assert(strcmp(sms.name, 'limlin'));  % make sure to get the right model
            assert(strcmp(sms.pnames{1}, 'Iinf'));
            assert(strcmp(sms.pnames{2}, 'c'));
            Iinfss = reshape(sms.mc.ss(:,1,:),1,[]);
            css = reshape(sms.mc.ss(:,2,:),1,[]);
            N95ss = (afrac/(1-afrac)) * Iinfss ./ css;
            meds(k) = prctile(N95ss, plotprctiles(3));
            plotPrctiles(k+sxhi(j), N95ss, hicols(j,:));
        end
        N95meds(s,:,2) = meds;  s = s + 1;
        plot((1:length(sessdorinames))+sxhi(j), meds, '-', 'Color', hicols(j,:));
    end
    % plot formatting
    set(gca,'Box','off','YScale','log',...
        'XTick',1:length(sessdorinames),'XTickLabel',sessdorinames,...
        'YTick',[1000:1000:9000 10000:10000:90000 100000:100000:500000],...
        'YTickLabel',{'10^3','','','','','','','','',...
                      '10^4','','','','','','','','',...
                      '10^5','','','',''});
    if i == 1
        xlabel('drift direction difference');
        ylabel('N_{95}');
    else
        set(gca, 'YTickLabel',[]);
    end
    xlim([0.75 length(sessdorinames)+0.2]);  ylim([1000 500000]);
    text(length(sessdorinames), 1000, sprintf('mouse %d', i), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');    
end
warning(unreswarning.state, 'MATLAB:dispatcher:UnresolvedFunctionHandle');

fprintf('Comparing Iinf(hi) vs. Iinf(lo) for different dtheta\n');
for k = 1:length(sessdorinames)
    Iinfdiff = squeeze(Iinfmeds(:,k,2) - Iinfmeds(:,k,1));
    [p,~,stats] = signrank(Iinfdiff);
    fprintf('%8s: <Ihi - Ilo> = %f  p-val = %0.2e zval = %d\n', ...
        sessdorinames{k}, mean(Iinfdiff), p, stats.zval);
end

fprintf('Comparing N95(hi) vs. N95(lo) for different dtheta\n');
for k = 1:length(sessdorinames)
    N95diff = squeeze(N95meds(:,k,2) - N95meds(:,k,1));
    [p,~,stats] = signrank(N95diff);
    fprintf('%8s: <N95hi - N95lo> = %f  p-val = %0.2e zval = %d\n', ...
        sessdorinames{k}, mean(N95diff), p, stats.zval);
end


%% save figure
if singlefig
    fprintf('\nWriting figure to figS4.pdf\n');
    print(['figs' filesep 'figS4'], '-dpdf');
end


function [In, varIn] = empTotalInf(d)
%% return total information in recorded population in given dataset
In = sum(mean(d.Iincr_samples,1));
varIn = sum(var(d.Iincr_samples,[],1));

function stats = mdbstrp_stat(In, varIn)
MC = 1e5; 
I = mvnrnd(In, diag(varIn), MC);
Iave_all = mean(In); 
I_null = mvnrnd(Iave_all* ones(size(In)), diag(varIn), MC);
Tstat = sum((I- mean(I, 1)).^2, 2);
Tstat_null = sum((I_null- mean(I_null, 1)).^2, 2);
stats.pvalg = sum(Tstat > Tstat_null)/ MC;
