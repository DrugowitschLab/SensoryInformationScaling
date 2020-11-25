function plotFigS1_S2
%% plots supplementary Fig. showing tuning curves and pairwise corr.

%% (general) settings
singlefig = true;   % set to true to have all planels in single figure
tcfitpath = ['.' filesep 'tuning_fits'];
dataset = 'm25b';
datasets = {{'m25a', 'm25b'}, {'m26a', 'm26b'}, ...
    {'aj42a', 'aj42b', 'aj42c', 'aj42d', 'aj42e'}, ...
    {'aj43a', 'aj43b', 'aj43c', 'aj43d', 'aj43e', 'aj43f', 'aj43g'}};

subsample = 'none';
conid = 1;  % contrast to use
eventRateScale  = 30;
ylimmax = 1;   % limit for the y-axis for tuning plots
ymulimmax = 1; % limit for the y-axis for sequential effects plots
uncol = [27 158 119] / 255;
dircol = [217 95 2] / 255;
oricol = [117 112 179] / 255;
oricols = [247 244 249; 129 144 190;  12  44 132;  78  22  98; ...
           145   0  63; 161   0  50; 177   0  38; 212 122 150] / 255;
tuningrows = 5;
tuningcols = 2;
cdfsubsamples = 10000;


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
set(groot, 'defaultFigureRenderer', 'painters');
if singlefig
    figure('Name', 'Figure', 'Units', 'centimeters', 'Position', [0 0 21 11.5]);
end


%% loading data
% activity data
fprintf('Loading %s ...\n', dataset');
[sr, oris, cons] = loaddata(dataset,subsample);
uoris = unique(oris);
ucons = unique(cons);
% tuning curve fits
tcfitfile = [tcfitpath filesep dataset filesep 'TuningCoef_Subsample_' subsample '.mat'];
fprintf('Loading %s ...\n', tcfitfile);
tcd = load(tcfitfile);
OriTun_vonMises = tcd.coef.OriTun_vonMises{conid};
DirTun_vonMises = tcd.coef.DirTun_vonMises{conid};
NullTun         = tcd.coef.NullFunc{conid};
x_OriTun    = tcd.coef.x_OriTun{conid,1};
x_DirTun    = tcd.coef.x_DirTun{conid,1};
orientationForTrials = oris;
averageSpkRateOverTrial = sr; 
nNeurons = size(averageSpkRateOverTrial, 2);
% compute R2 of fits
R2_OriTun = NaN(1, nNeurons);
R2_DirTun = NaN(1, nNeurons);
f = NaN(1, length(uoris));
for nn = 1: nNeurons
    for orii = 1:length(uoris)
        f(orii) = mean(sr(oris == uoris(orii) & cons == ucons(conid), nn));
    end
    f_ori = OriTun_vonMises(x_OriTun(nn, :), uoris)' / eventRateScale;
    f_dir = DirTun_vonMises(x_DirTun(nn, :), uoris)' / eventRateScale;
    R2_OriTun(nn) = 1 - sum((f_ori - f).^2) / sum((f - mean(f)).^2);
    R2_DirTun(nn) = 1 - sum((f_dir - f).^2) / sum((f - mean(f)).^2);
    R2_OriTun(nn)   = corr(OriTun_vonMises(x_OriTun(nn, :), orientationForTrials),  averageSpkRateOverTrial(:, nn))^2;
    R2_DirTun(nn)   = corr(DirTun_vonMises(x_DirTun(nn, :), orientationForTrials), averageSpkRateOverTrial(:, nn))^2;
end
% separate neurons into categories by fit significance
nonsigneurons = find(tcd.coef.alpha_OriVsNull{conid} >= 0.05);
orineurons = find(tcd.coef.alpha_OriVsNull{conid} < 0.05 & ...
    tcd.coef.alpha_OriVsDir{conid} > 0.05);
dirneurons = find(tcd.coef.alpha_OriVsNull{conid} < 0.05 & ...
    tcd.coef.alpha_OriVsDir{conid} < 0.05);


%% plot activity / tuning curves of untuned and tuned neurons for example data
k = tuningrows * tuningcols;
thetas          = linspace(0, 360, 361);
exampleids_untund = 1:k; % randsample([1: length(nonsigneurons)],k);
[~, exampleids_ori] = sort(R2_OriTun(orineurons), 'descend');
exampleids_ori = exampleids_ori(1:k); % randsample(exampleids_ori(1: k), k);
[~, exampleids_dir] = sort(R2_DirTun(dirneurons), 'descend');
exampleids_dir = exampleids_dir(1:k); % randsample(exampleids_dir(1: k), k);
if ~singlefig; figure; hold on; end
for jj = 1: k
    jr = mod(jj-1,tuningrows)+1;
    jc = ceil(jj/tuningrows);
    if singlefig
        subplotcm([(1+(jc-1)*2) (0.8+(tuningrows-jr+1)*1.5) 1.5 1]);  hold on;
    else, subplot(tuningrows,ceil(k/tuningrows), jj);
    end
    i = nonsigneurons(exampleids_untund(jj));
    plotntuning(sr, oris, cons, i, conid, thetas, ...
        NullTun(tcd.coef.x_null{conid}(i,:), thetas) / eventRateScale, uncol);
    ylim([0 ylimmax]);
    set(gca,'Box','off','XTick',[0 max(uoris)],'XTickLabels',[],'YTick',[0 0.1]);
    if jc > 1, set(gca, 'YColor','none');  end
end
if ~singlefig; figure; hold on; end
for jj = 1: k
    jr = mod(jj-1,tuningrows)+1;
    jc = ceil(jj/tuningrows);
    if singlefig
        subplotcm([(1+(tuningcols+jc-1)*2) (0.8+(tuningrows-jr+1)*1.5) 1.5 1]);  hold on;
    else, subplot(tuningrows,ceil(k/tuningrows), jj);
    end
    i = orineurons(exampleids_ori(jj));
    plotntuning(sr, oris, cons, i, conid, thetas, ...
        OriTun_vonMises(tcd.coef.x_OriTun{conid}(i,:), thetas) / eventRateScale, oricol);
    ylim([0 ylimmax]);
    set(gca,'Box','off','XTick',[0 max(uoris)],'XTickLabels',[],'YColor','none');
end
if ~singlefig; figure; hold on; end
for jj = 1: k
    jr = mod(jj-1,tuningrows)+1;
    jc = ceil(jj/tuningrows);
    if singlefig
        subplotcm([(1+(2*tuningcols+jc-1)*2) (0.8+(tuningrows-jr+1)*1.5) 1.5 1]);  hold on;
    else, subplot(tuningrows,ceil(k/tuningrows), jj);
    end
    i = dirneurons(exampleids_dir(jj));
    plotntuning(sr, oris, cons, i, conid, thetas, ...
        DirTun_vonMises(tcd.coef.x_DirTun{conid}(i,:), thetas) / eventRateScale, dircol);
    ylim([0 ylimmax]);
    set(gca,'Box','off','XTick',[0 max(uoris)],'XTickLabels',[],'YColor','none');
end

%% plot adaptation
poris = [0; uoris];
ushift = [-4.5 -3.5 -2.5 -1.5 -0.5 0.5 1.5 2.5 3.5 4.5] * 2;
if ~singlefig; figure; hold on; end
for jj = 1:tuningcols
    if singlefig
        subplotcm([(1+(jj-1)*2) 0.8 1.5 1]);  hold on;        
    else, subplot(2,1,jj);  hold on;
    end
    [msr, ssr] = oriresp(sr, oris, nonsigneurons(exampleids_untund((jj-1)*tuningrows+1)));
    for ii = 1:length(uoris)
        for kk=1:(length(uoris)+1)
            if kk==1, m = msr(ii,end) + [-1 1]*ssr(ii,end);
            else, m = msr(ii,kk-1) + [-1 1]*ssr(ii,kk-1); end
            plot((poris(kk)+ushift(ii))*[1 1], m, '-', 'Color', oricols(ii,:), 'LineWidth', 0.25);
        end
    end
    for ii = 1:(length(uoris)+1)
        if ii > length(uoris), c = [0 0 0]; lw = 1;
        else, c = oricols(ii,:); lw = 0.5; end
        plot(poris, [msr(ii,end) msr(ii,:)], '-', 'Color', c, 'LineWidth', lw);
    end
    xlim([-3 363]);  ylim([0 ymulimmax]);
    set(gca,'Box','off','XTick',[0 max(uoris)],'XTickLabels',{'0','360'},'YTick',[0 0.1]);
    if jj == 1, xlabel('drift direction');  else, set(gca, 'YColor','none'); end
end
if ~singlefig; figure; hold on; end
for jj = 1:2
    if singlefig
        subplotcm([(1+(tuningcols+jj-1)*2) 0.8 1.5 1]);  hold on;        
    else, subplot(2,1,jj);  hold on;
    end
    [msr, ssr] = oriresp(sr, oris, orineurons(exampleids_ori((jj-1)*tuningrows+1)));
    for ii = 1:length(uoris)
        for kk=1:(length(uoris)+1)
            if kk==1, m = msr(ii,end) + [-1 1]*ssr(ii,end);
            else, m = msr(ii,kk-1) + [-1 1]*ssr(ii,kk-1); end
            plot((poris(kk)+ushift(ii))*[1 1], m, '-', 'Color', oricols(ii,:), 'LineWidth', 0.25);
        end
    end
    for ii = 1:(length(uoris)+1)
        if ii > length(uoris), c = [0 0 0]; lw = 1;
        else, c = oricols(ii,:); lw = 0.5; end
        plot(poris, [msr(ii,end) msr(ii,:)], '-', 'Colo', c, 'LineWidth', lw);
    end
    xlim([-3 363]);  ylim([0 ymulimmax]);
    set(gca,'Box','off','XTick',[0 max(uoris)],'XTickLabels',{'0','360'},'YColor','none');
end
if ~singlefig; figure; hold on; end
for jj = 1:2
    if singlefig
        subplotcm([(1+(tuningcols*2+jj-1)*2) 0.8 1.5 1]);  hold on;        
    else, subplot(2,1,jj);  hold on;
    end
    [msr, ssr] = oriresp(sr, oris, dirneurons(exampleids_dir((jj-1)*tuningrows+1)));
    for ii = 1:length(uoris)
        for kk=1:(length(uoris)+1)
            if kk==1, m = msr(ii,end) + [-1 1]*ssr(ii,end);
            else, m = msr(ii,kk-1) + [-1 1]*ssr(ii,kk-1); end
            plot((poris(kk)+ushift(ii))*[1 1], m, '-', 'Color', oricols(ii,:), 'LineWidth', 0.25);
        end
    end
    for ii = 1:(length(uoris)+1)
        if ii > length(uoris), c = [0 0 0]; lw = 1;
        else, c = oricols(ii,:); lw = 0.5; end
        plot(poris, [msr(ii,end) msr(ii,:)], '-', 'Colo', c, 'LineWidth', lw);
    end
    xlim([-3 363]);  ylim([0 ymulimmax]);
    set(gca,'Box','off','XTick',[0 max(uoris)],'XTickLabels',{'0','360'},'YColor','none');
end


%% plot pairwise correlations and R2 across animals
animals = length(datasets);
fiwidth = 2;
ledge = 1.5+6*tuningcols;           % left edge of animals panel
fheight = 1.5*(tuningrows+1)-0.5;  % overall height across animals
fiheight = (fheight-(animals-1)*0.5)/animals;  % height of per-animal panel
for animal = 1:animals
    sessions = length(datasets{animal});
    bedge = 0.8+fheight-animal*fiheight-(animal-1)*0.5;
    % load and process data
    c = cell(1,sessions);
    cshuf = cell(1,sessions);
    cun = cell(1,sessions);
    cori = cell(1,sessions);
    cdir = cell(1,sessions);
    R2ori = cell(1,sessions);
    R2dir = cell(1,sessions);
    for sessi = 1:sessions
        dsi = datasets{animal}{sessi};
        [sr, oris, cons] = loaddata(dsi,subsample);
        tcfitfile = [tcfitpath filesep dsi filesep 'TuningCoef_Subsample_' subsample '.mat'];
        uoris = unique(oris);
        ucons = unique(cons);
        fprintf('Loading %s ...\n', tcfitfile);
        tcd = load(tcfitfile);
        OriTun_vonMises = tcd.coef.OriTun_vonMises{conid};
        DirTun_vonMises = tcd.coef.DirTun_vonMises{conid};
        x_OriTun    = tcd.coef.x_OriTun{conid,1};
        x_DirTun    = tcd.coef.x_DirTun{conid,1};
        orientationForTrials = oris;
        averageSpkRateOverTrial = sr; 
        nNeurons = size(averageSpkRateOverTrial, 2);
        % compute R2 of fits
        R2_OriTun = NaN(1, nNeurons);
        R2_DirTun = NaN(1, nNeurons);
        f = NaN(1, length(uoris));
        for nn = 1: nNeurons
            for orii = 1:length(uoris)
                f(orii) = mean(sr(oris == uoris(orii) & cons == ucons(conid), nn));
            end
            f_ori = OriTun_vonMises(x_OriTun(nn, :), uoris)' / eventRateScale;
            f_dir = DirTun_vonMises(x_DirTun(nn, :), uoris)' / eventRateScale;
            R2_OriTun(nn) = 1 - sum((f_ori - f).^2) / sum((f - mean(f)).^2);
            R2_DirTun(nn) = 1 - sum((f_dir - f).^2) / sum((f - mean(f)).^2);
            %R2_OriTun(nn)   = corr(OriTun_vonMises(x_OriTun(nn, :), orientationForTrials),  averageSpkRateOverTrial(:, nn))^2;
            %R2_DirTun(nn)   = corr(DirTun_vonMises(x_DirTun(nn, :), orientationForTrials), averageSpkRateOverTrial(:, nn))^2;
        end
        % separate neurons into categories by fit significance
        nonsigneurons = find(tcd.coef.alpha_OriVsNull{conid} >= 0.05);
        orineurons = find(tcd.coef.alpha_OriVsNull{conid} < 0.05 & ...
            tcd.coef.alpha_OriVsDir{conid} > 0.05);
        dirneurons = find(tcd.coef.alpha_OriVsNull{conid} < 0.05 & ...
            tcd.coef.alpha_OriVsDir{conid} < 0.05);
        % compute statistics
        c{sessi} = pairwisecorr(sr);
        cshuf{sessi} = pairwisecorr(shufsr(sr, oris, cons));
        cun{sessi} = pairwisecorr(sr(:, nonsigneurons));
        cori{sessi} = pairwisecorr(sr(:, orineurons));
        cdir{sessi} = pairwisecorr(sr(:, dirneurons));
        R2ori{sessi} = R2_OriTun(orineurons);
        R2dir{sessi} = R2_DirTun(dirneurons);
    end

    % pairwise correlations, all vs. shuffled
    if singlefig
        subplotcm([ledge bedge fiwidth fiheight]);  hold on;
    else
        figure;  hold on;
    end
    for sessi = 1:sessions
        if strcmp(dataset, datasets{animal}{sessi}), lw = 1; else, lw = 0.5; end
        [f, x] = subsampledecdf(c{sessi}, cdfsubsamples);
        stairs([0 x' 1],[0 f' 1],'Color',[0 0 0],'LineWidth',lw);
        [f, x] = subsampledecdf(cshuf{sessi}, cdfsubsamples);
        stairs([0 x' 1],[0 f' 1],'Color',[0.8 0 0],'LineWidth',lw);
        xlim([-0.2 0.4]);  ylim([0 1]);
        plot([1 1]*mean(c{sessi}), ylim, '-', 'LineWidth', 0.5, 'Color', [1 1 1]*0.5);
        plot([1 1]*mean(cshuf{sessi}), ylim, '-', 'LineWidth', 0.5, 'Color', [0.9 0.5 0.5]);
    end
    xlim([-0.2 0.4]);  ylim([0 1]);
    set(gca,'Box','off','XTick',[-0.2 0 0.2 0.4],'YTick',[0 1]);
    text(0.4, 0.05, sprintf('mouse %d', animal), ...
        'HorizontalAlignment','right','VerticalAlignment','bottom');
    if animal == animals
        xlabel('correlation coefficient \rho');
        ylabel('cumulative probability');
    else, set(gca,'XTickLabels',[]);
    end
    
    
    % pairwise correlations, untuned & tuned
    if singlefig
        subplotcm([ledge+0.5+fiwidth bedge fiwidth fiheight]);
        hold on;
    else
        figure;  hold on;
    end
    for sessi = 1:sessions
        if strcmp(dataset, datasets{animal}{sessi}), lw = 1; else, lw = 0.5; end
        [f, x] = subsampledecdf(cun{sessi}, cdfsubsamples);
        stairs([0 x' 1],[0 f' 1],'Color',uncol,'LineWidth',lw);
        [f, x] = subsampledecdf(cori{sessi}, cdfsubsamples);
        stairs([0 x' 1],[0 f' 1],'Color',oricol,'LineWidth',lw);
        [f, x] = subsampledecdf(cdir{sessi}, cdfsubsamples);
        stairs([0 x' 1],[0 f' 1],'Color',dircol,'LineWidth',lw);
        xlim([-0.2 0.4]);  ylim([0 1]);
        plot([1 1]*mean(cun{sessi}), ylim, '-', 'LineWidth', 0.5, 'Color', 0.5+0.5*uncol);
        plot([1 1]*mean(cori{sessi}), ylim, '-', 'LineWidth', 0.5, 'Color', 0.5+0.5*oricol);
        plot([1 1]*mean(cdir{sessi}), ylim, '-', 'LineWidth', 0.5, 'Color', 0.5+0.5*dircol);
    end
    set(gca,'Box','off','XTick',[-0.2 0 0.2 0.4],'YTick',[0 1]);
    if animal == animals, xlabel('correlation coefficient \rho');
    else set(gca,'XTickLabels',[]); end
    
    % R2 for fitted neurons
    if singlefig
        subplotcm([ledge+1+2*fiwidth bedge fiwidth fiheight]);
        hold on;
    else
        figure;  hold on;
    end
    for sessi = 1:sessions
        if strcmp(dataset, datasets{animal}{sessi}), lw = 1; else, lw = 0.5; end
        [f,x] = ecdf(R2ori{sessi});
        stairs([0 x' 1],[0 f' 1],'Color',oricol,'LineWidth',lw);
        [f,x] = ecdf(R2dir{sessi});
        stairs([0 x' 1],[0 f' 1],'Color',dircol,'LineWidth',lw);
        plot([1 1]*mean(R2ori{sessi}), ylim, '-', 'LineWidth', 0.5, 'Color', 0.5+0.5*oricol);
        plot([1 1]*mean(R2dir{sessi}), ylim, '-', 'LineWidth', 0.5, 'Color', 0.5+0.5*dircol);
    end
    xlim([0 1]);  ylim([0 1]);        
    set(gca,'Box','off','XTick',[0 0.25 0.5 0.75 1],'YTick',[0 1]);
    if animal == animals, xlabel('R^2');
    else set(gca,'XTickLabels',[]); end
end


%% write plot to file
if singlefig
    fprintf('\nWriting figure to figS1_S2.pdf\n');
    print(['figs' filesep 'figS1_S2'], '-dpdf');
end


function plotntuning(sr, oris, cons, neurid, conid, thetas, fthetas, col)
%% plots neural responses as well as fitted tuning curve

jittersd = 1;
shadedcol = col*0.8+[1 1 1]*0.2;

[sr_ori_ave, ~, sr_ori_neg, sr_ori_pos] = errbar_util(sr(:, neurid), oris);

consu = unique(cons);
orisu = [0 unique(oris)'];
hold on;
for i = 1:length(orisu)
    j = cons == consu(conid) & mod(oris, 360) == mod(orisu(i), 360);
    plot(orisu(i) + jittersd * randn(1, sum(j)), sr(j,neurid)','o',...
        'MarkerFaceColor',shadedcol, ...
        'MarkerEdgeColor','None','MarkerSize',0.5); alpha(0.1);
    if i == 1, j = length(orisu) - 1; else j = i-1; end
    plot(orisu(i)*[1 1], [sr_ori_neg(j) sr_ori_pos(j)], '-', ...
        'Color', col, 'LineWidth', 1.5);
end

plot(orisu, [sr_ori_ave(end); sr_ori_ave], 'o', 'MarkerSize', 3, ...
    'MarkerFaceColor', col, 'MarkerEdgeColor', 'None');
plot(thetas, fthetas, '-', 'LineWidth', 0.5, 'Color', shadedcol);
xlim([(min(thetas)-3*jittersd) (max(thetas)+3*jittersd)]);

    
function [sr_ori_ave, sr_ori_std, sr_ori_neg, sr_ori_pos] = errbar_util(sr, ori)
uoris           = sort(unique(ori), 'ascend');
orin            = length(uoris);
sr_ori_ave      = nan(orin, 1);
sr_ori_std      = nan(orin, 1);
sr_ori_neg      = nan(orin, 1);
sr_ori_pos      = nan(orin, 1);    
for ii = 1: orin
    sr_ori_tmp  = sr(ori == uoris(ii));
    sr_ori_ave(ii)  = mean(sr_ori_tmp);
    sr_ori_std(ii)  = std(sr_ori_tmp);
    sr_ori_neg(ii)  = prctile(sr_ori_tmp, 25);
    sr_ori_pos(ii)  = prctile(sr_ori_tmp, 75);
end


function sr = shufsr(sr, oris, cons)
%% trial-shuffle data within each condition
N = size(sr, 2);
ucons = unique(cons);  conn = length(ucons);
uoris = unique(oris);  orin = length(uoris);
for coni = 1:conn
    for orii = 1:orin
        condtrials = find((oris == uoris(orii)) & (cons == ucons(coni)));
        T = length(condtrials);
        for n = 2:N
            sr(condtrials,n) = sr(condtrials(randperm(T)),n);
        end
    end
end


function [f,x] = subsampledecdf(x, subsamples)
if length(x) <= subsamples
    [f, x] = ecdf(x);
else
    [f, x] = ecdf(x(randperm(length(x), subsamples)));
end


function c = pairwisecorr(sr)
c = corr(sr); 
c = c(logical(triu(ones(size(sr, 2)), 1)));


function [msr, ssr] = oriresp(sr, oris, n)
%% computes mean responses for neuron n conditional on previous orientation

uoris = unique(oris);  orin = length(uoris);
msr = NaN(orin + 1, orin);
ssr = NaN(orin + 1, orin);
prevoris = [-1; oris(1:(end-1))];
for i = 1:(orin+1)  % prev. orientation
    if i > orin, previ = true(size(sr,1),1);  % last i across all trials
    else, previ = prevoris == uoris(i);  end
    for j = 1:orin  % curr. orientation
        msr(i,j) = mean(sr(oris == uoris(j) & previ, n));
        ssr(i,j) = sqrt(var(sr(oris == uoris(j) & previ, n)) / sum(oris == uoris(j) & previ));
    end
end
