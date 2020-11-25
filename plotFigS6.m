function plotFigS6
%% plots SI figure showing model comparison for info scaling fits

%% (general) settings
singlefig = true;
fitpath = ['.' filesep 'fits'];
datasets = {'m25a', 'm25b', 'm26a', 'm26b', ...
    'aj42a', 'aj42b', 'aj42c', 'aj42d', 'aj42e', ...
    'aj43a', 'aj43b', 'aj43c', 'aj43d', 'aj43e', 'aj43f', 'aj43g'};
discriminations = {'o1-2_c1', 'o1-8_c1', 'o2-3_c1', 'o3-4_c1', ...
                   'o4-5_c1', 'o5-6_c1', 'o6-7_c1', 'o7-8_c1'};               
waic_field = 'WAIC1';


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


%% collect statistics
unreswarning = warning('query', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
waics = NaN(length(datasets), length(discriminations), 3);
waic_shuffs = NaN(length(datasets), length(discriminations), 2);
for i = 1:length(datasets)
    for j = 1:length(discriminations)
        % non-shuffled
        fitfile = [fitpath filesep datasets{i} '_' discriminations{j} '_norm.mat'];
        fprintf('Loading %s ...\n', fitfile);
        m = load(fitfile);
        linm = m.sms{1};
        limlinm = m.sms{2};
        limexpm = m.sms{3};
        assert(strcmp(linm.name, 'lin'));    % check models
        assert(strcmp(limlinm.name, 'limlin'));
        assert(strcmp(limexpm.name, 'limexp'));
        waics(i,j,1) = linm.mc.(waic_field);
        waics(i,j,2) = limlinm.mc.(waic_field);
        waics(i,j,3) = limexpm.mc.(waic_field);
        % shuffled
        fitfile = [fitpath filesep datasets{i} '_' discriminations{j} '_shuf_norm.mat'];
        fprintf('Loading %s ...\n', fitfile);
        m = load(fitfile);
        linm = m.sms{1};
        limlinm = m.sms{2};
        assert(strcmp(linm.name, 'lin'));    % check models
        assert(strcmp(limlinm.name, 'limlin'));
        waic_shuffs(i,j,1) = linm.mc.(waic_field);
        waic_shuffs(i,j,2) = limlinm.mc.(waic_field);
    end
end
warning(unreswarning.state, 'MATLAB:dispatcher:UnresolvedFunctionHandle');


%% analyze collected statistics
fprintf('\nPerforming statistic comparisons (all Wilcoxon signed rank)\n');
md = waics(:,:,1) - waics(:,:,2);
mdexp = waics(:,:,3) - waics(:,:,1);
mdshuffs = waic_shuffs(:,:,1) - waic_shuffs(:,:,2);
for i = 1:length(datasets)
    p = signrank(md(i,:));
    % does non-limiting model fit better than limiting model?
    if mean(md(i,:)) > 0.0
        if p > 0.05
            fprintf('%5s: lim > lin, but p = %f\n', datasets{i}, p);
        end
    else
        if p < 0.05
            fprintf('%5s: lim < lin, and p = %f\n', datasets{i}, p);
        else
            fprintf('%5s: lim > lin, but p = %f\n', datasets{i}, p);
        end
    end
    % for shuffled, does limiting model fit better than non-limiting model?
    p = signrank(mdshuffs(i,:));
    if mean(mdshuffs(i,:)) < 0.0
        if p > 0.05
            fprintf('%5s: shuflim < shuflin, but p = %f\n', datasets{i}, p);
        end
    else
        if p < 0.05
            fprintf('%5s: shuflim > shuflin, and p = %f\n', datasets{i}, p);
        else
            fprintf('%5s: shuflim > shuflin, but p = %f\n', datasets{i}, p);
        end
    end
    % does exponential limiting model fit better than linear limiting
    % model?
    p = signrank(mdexp(i,:));
    if mean(mdexp(i,:)) > 0.0
        if p > 0.05
            fprintf('%5s: lim > exp, but p = %f\n', datasets{i}, p);
        end
    else
        if p < 0.05
            fprintf('%5s: lim < exp, and p = %f\n', datasets{i}, p);
        else
            fprintf('%5s: lim > exp, but p = %f\n', datasets{i}, p);
        end
    end
end


%% plot collected statistics
% comparison to linear limited model
if singlefig
    figure('Name', 'Figure', 'Units', 'centimeters', 'Position', [0 0 10.7 5.1]);
    subplotcm([0.5 0.8 4.6 3.5]);  hold on;
    text(0,3.5,'a','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
plot([1 1]*prctile(reshape(md,1,[]),50), [0 80], '-', 'Color', [0.5 0.5 0.9]);
plot([1 1]*prctile(reshape(mdshuffs,1,[]),50), [0 80], '-', 'Color', [0.9 0.5 0.5]);
histogram(reshape(md, 1, []), 'BinWidth', 0.1, ...
    'EdgeColor', 'none', 'FaceColor', [0 0 0.8]);
histogram(reshape(mdshuffs, 1, []), 'BinWidth', 0.1, ...
    'EdgeColor', 'none', 'FaceColor', [0.8 0 0]);
plot([0 0], [0 80], '-', 'Color', [1 1 1]*0.5);
xlabel('WAIC_{non-lim} - WAIC_{lim}');
xlim([-1 3.5]);  ylim([0 80]);
set(gca,'YColor','none','YTick',[],'Box','off');

% comparison exponential limited vs. linear limited model
if singlefig
    subplotcm([5.6 0.8 4.6 3.5]);  hold on;
    text(0,3.5,'b','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');    
else
    figure;  hold on;
end
plot([1 1]*prctile(reshape(mdexp,1,[]),50), [0 80], '-', 'Color', [0.5 0.9 0.5]);
histogram(reshape(mdexp, 1, []), 'BinWidth', 0.1, ...
    'EdgeColor', 'none', 'FaceColor', [0 0.8 0], 'BinWidth', 0.1);
plot([0 0], [0 80], '-', 'Color', [1 1 1]*0.5);
xlabel('WAIC_{lim-exp} - WAIC_{lim}');
xlim([-1 3.5]);  ylim([0 80]);
set(gca,'YColor','none','YTick',[],'Box','off');


%% write figure to file
if singlefig
    fprintf('\nWriting figure to figS6.pdf\n');
    print(['figs' filesep 'figS6'], '-dpdf');
end
