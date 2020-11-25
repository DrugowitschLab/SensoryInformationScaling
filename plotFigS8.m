function plotFigS8
%% plots SI figure showing statistics of linear fits to 1/IN over 1/N

%% general settings
singlefig = true;   % set to true to have all planels in single figure
momentpath = ['.' filesep 'moment_cache'];
datasets = {'m25a', 'm25b', 'm26a', 'm26b', ...
    'aj42a', 'aj42b', 'aj42c', 'aj42d', 'aj42e', ...
    'aj43a', 'aj43b', 'aj43c', 'aj43d', 'aj43e', 'aj43f', 'aj43g'};
%datasets = {'m25a', 'm25b', 'm26a', 'm26b'};
discrs = {'o1-2','o2-3','o3-4','o4-5','o5-6','o6-7','o7-8','o1-8'};
discrstats = [true false true false true false true false];
sessions = length(datasets);


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
    figure('Name', 'Figure', 'Units', 'centimeters', 'Position', [0 0 9 4.8]);
end


%% collect statistics
R2adj = NaN(sessions, length(discrs));
betas = NaN(sessions, length(discrs), 6);  % [b0 (b0CI x 2) b1 (b1CI x 2)]
sigb0 = NaN(sessions, length(discrs));
for i = 1:sessions
    for j = 1:length(discrs)
        momfile = [momentpath filesep datasets{i} '_' discrs{j} '_c1.mat'];
        fprintf('Loading %s ...\n', momfile);
        d = load(momfile);
        
        % moments of 1/IN, using approx. <1/x>=1/<x>, var(1/x)=var(x)/<x>^4
        I_mu = cumsum(mean(d.Iincr_samples, 1));
        I_var = cumsum(var(d.Iincr_samples, [], 1));
        invI_mu = 1 ./ I_mu;
        invI_var = I_var ./ I_mu.^4;
        invNs = 1./(1:length(I_mu));
        
        % fit model and compute statistics
        mdl = fitlm(invNs', invI_mu', 'Weights', (1./invI_var)');
        mbetas = mdl.Coefficients.Estimate;
        mCIs = coefCI(mdl);
        betas(i,j,:) = [mbetas(1) mCIs(1,:) mbetas(2) mCIs(2,:)];
        R2adj(i,j) = mdl.Rsquared.Adjusted;
        sigb0(i,j) = mdl.Coefficients.pValue(1);
    end    
end
% compute stats for all b0
fprintf('\nAvg. R2adj = %f\n', mean(reshape(R2adj,1,[])));
fprintf('\nT-test for b0 == 0 across non-overlapping discriminations:\n');
b0test = reshape(betas(:,discrstats,1),[],1);
[~,p,~,stats] = ttest(b0test);
fprintf('<b0> = %f,  t(%d) = %f,  p = %f\n', ...
    mean(b0test), stats.df, stats.tstat, p);


%% plot coefficients against each other
if singlefig
    subplotcm([5 0.8 3.5 3.5]);  hold on;
    text(-1,4,'b','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
% plot CIs
for i = 1:sessions
    for j = 1:length(discrs)
        plot(squeeze(betas(i,j,5:6)), [1 1]*betas(i,j,1), '-', 'Color', [1 1 1]*0.5);
        plot([1 1]*betas(i,j,4), squeeze(betas(i,j,2:3)), '-', 'Color', [1 1 1]*0.5);
    end
end
% significant and non-significant b0's
b0flat = reshape(betas(:,:,1), 1, []);
b1flat = reshape(betas(:,:,4), 1, []);
sigb0flat = reshape(sigb0, 1, []) < 0.05;
plot(b1flat(~sigb0flat), b0flat(~sigb0flat), 'o', 'MarkerSize', 4, ...
    'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0]);
plot(b1flat(sigb0flat), b0flat(sigb0flat), 'o', 'MarkerSize', 4, ...
    'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
xlim([0 30]);  ylim([0 0.07]);
xlabel('\beta_1');
ylabel('\beta_0');
set(gca,'Box','off','XTick',0:10:30, 'YTick',0:0.01:0.07);


%% plot R2adj histogram
if singlefig
    subplotcm([0 0.8 3.5 3.5]);  hold on;
    text(0,4,'a','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
histogram(reshape(R2adj, 1, []), 'BinWidth', 0.00025, ...
    'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
xlabel('R^2_{adj}');  xlim([0.9945 1]);
set(gca,'Box','off','XTick',0.995:0.001:1,'YColor','none','YTick',[]);


%% write figure to file
if singlefig
    fprintf('\nWriting figure to figS8.pdf\n');
    print(['figs' filesep 'figS8'], '-dpdf');
end
