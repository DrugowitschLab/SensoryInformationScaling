function plotFigS11
%% plots panels of figure 5


%% (general) settings
singlefig = true;   % set to true to have all planels in single figure
orderpath = ['.' filesep 'reordering_cache'];
Nsamples = 1000;
dataset = 'm25b';
oricomp = [1 2];       % index within dori.oricombs
oricols = [     0    0.4470    0.7410; ...
           0.8500    0.3250    0.0980];
avgcol = [0.5000    0.8000    0.1000];
randcol = [1 1 1]*0.2;
Ifrac = 0.90;
doriid = 1;
doristats = [1 4 6 7]; % to keep in dori.oricombs for stats
Nprctiles = [50 2.5 97.5];
datasets = {'m25a', 'm25b', 'm26a', 'm26b', ...
    'aj42a', 'aj42b', 'aj42c', 'aj42d', 'aj42e', ...
    'aj43a', 'aj43b', 'aj43c', 'aj43d', 'aj43e', 'aj43f', 'aj43g'};
datasuffix = '_shuf'; % change to '_shuf' to use shuffled data

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
    figure('Name', 'Figure', 'Units', 'centimeters', 'Position', [0 0 15.5 5.1]);
end


%% load and process example data
orderfile = sprintf('%s%s%s_dori%d_c1%s.mat', ...
    orderpath, filesep, dataset, doriid, datasuffix);
fprintf('Loading %s ...\n', orderfile);
dord = load(orderfile);
N = dord.N;
% process data for first discrimination
In = dord.In_rand(:,:,oricomp(1));
In_frac = bsxfun(@rdivide, In, In(:,end));
Nrand_samples = sum(bsxfun(@times, 2:N, diff(In_frac >= Ifrac, [], 2)), 2)';
% collect info growth, scaled to max information
Infrac = [dord.In(oricomp(1),:,oricomp(1));
          dord.In(oricomp(1),:,oricomp(2));
          dord.In(oricomp(1),:,end);
          mean(In,1)];
Infrac_sd = sqrt([dord.In_var(oricomp(1),:,oricomp(1));
                  dord.In_var(oricomp(1),:,oricomp(2));
                  dord.In_var(oricomp(1),:,end);
                  var(In,[],1)]);
Infrac_sd = bsxfun(@rdivide, Infrac_sd, Infrac(1,end));
Infrac = bsxfun(@rdivide, Infrac, Infrac(1,end));
% N samples for optimized and random ordering
Nord_samples = ...
    sampleNfrac(dord.In(oricomp(1),:,end), dord.In_var(oricomp(1),:,end), ...
    Ifrac, size(In,1));
Nrand_prctiles = prctile(Nrand_samples, Nprctiles);
Nord_prctiles = prctile(Nord_samples, Nprctiles);
% plot error shares first
if singlefig
    subplotcm([1.2 0.8 4 3]);  hold on;
    text(-1,4,'a','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure; hold on;
end
alpha(errorshade(0:N, [0 Infrac(1,:)], [0 Infrac_sd(1,:)], oricols(1,:)), 0.5);
alpha(errorshade(0:N, [0 Infrac(2,:)], [0 Infrac_sd(2,:)], oricols(2,:)), 0.5);
alpha(errorshade(0:N, [0 Infrac(3,:)], [0 Infrac_sd(3,:)], avgcol), 0.5);
alpha(errorshade(0:N, [0 Infrac(4,:)], [0 Infrac_sd(4,:)], randcol), 0.1);
plot(0:N, [0 Infrac(1,:)], '-', 'LineWidth', 1, 'Color', oricols(1,:));
plot(0:N, [0 Infrac(2,:)], '-', 'LineWidth', 1, 'Color', oricols(2,:));
plot(0:N, [0 Infrac(3,:)], '-', 'LineWidth', 1, 'Color', avgcol);
plot(0:N, [0 Infrac(4,:)], '-', 'LineWidth', 1, 'Color', randcol);
plot([0 N], [1 1]*Ifrac, '-', 'Color', [1 1 1]*0.5);
plot(Nrand_prctiles(2:3), [1 1]*(Ifrac+0.005), '-', 'Color', randcol);
plot(Nord_prctiles(2:3), [1 1]*(Ifrac-0.005), '-', 'Color', avgcol);
plot(Nrand_prctiles(1), Ifrac+0.005, 'o', 'MarkerSize', 4, ...
    'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', randcol);
plot(Nord_prctiles(1), Ifrac-0.005, 'o', 'MarkerSize', 4, ...
    'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', avgcol);
xlim([0 N]);
ylim([0 1.05]);
set(gca,'Box','off','XTick',[0 100 200 300],'YTick',[0 0.25 0.5 0.75 Ifrac 1],...
    'YTickLabel',{'0%','25%','50%','75%','90%','100%'});
xlabel('Number of neurons N');
ylabel('Fraction of total population information [%]');
% process data for second discrimination
In = dord.In_rand(:,:,oricomp(2));
In_frac = bsxfun(@rdivide, In, In(:,end));
Nrand_samples = sum(bsxfun(@times, 2:N, diff(In_frac >= Ifrac, [], 2)), 2)';
% collect info growth, scaled to max information
Infrac = [dord.In(oricomp(2),:,oricomp(2));
          dord.In(oricomp(2),:,oricomp(1));
          dord.In(oricomp(2),:,end);
          mean(In,1)];
Infrac_sd = sqrt([dord.In_var(oricomp(2),:,oricomp(2));
                  dord.In_var(oricomp(2),:,oricomp(1));
                  dord.In_var(oricomp(2),:,end);
                  var(In,[],1)]);
Infrac_sd = bsxfun(@rdivide, Infrac_sd, Infrac(1,end));
Infrac = bsxfun(@rdivide, Infrac, Infrac(1,end));
% N samples for optimized and random ordering
Nord_samples = ...
    sampleNfrac(dord.In(oricomp(2),:,end), dord.In_var(oricomp(2),:,end), Ifrac, Nsamples);
Nrand_prctiles = prctile(Nrand_samples, Nprctiles);
Nord_prctiles = prctile(Nord_samples, Nprctiles);
% plot error shares first
if singlefig
    subplotcm([5.5  0.8 4 3]);  hold on;
else
    figure; hold on;
end
alpha(errorshade(0:N, [0 Infrac(1,:)], [0 Infrac_sd(1,:)], oricols(2,:)), 0.5);
alpha(errorshade(0:N, [0 Infrac(2,:)], [0 Infrac_sd(2,:)], oricols(1,:)), 0.5);
alpha(errorshade(0:N, [0 Infrac(3,:)], [0 Infrac_sd(3,:)], avgcol), 0.5);
alpha(errorshade(0:N, [0 Infrac(4,:)], [0 Infrac_sd(4,:)], randcol), 0.1);
p1 = plot(0:N, [0 Infrac(1,:)], '-', 'LineWidth', 1, 'Color', oricols(2,:));
p2 = plot(0:N, [0 Infrac(2,:)], '-', 'LineWidth', 1, 'Color', oricols(1,:));
p3 = plot(0:N, [0 Infrac(3,:)], '-', 'LineWidth', 1, 'Color', avgcol);
p4 = plot(0:N, [0 Infrac(4,:)], '-', 'LineWidth', 1, 'Color', randcol);
plot([0 N], [1 1]*Ifrac, '-', 'Color', [1 1 1]*0.5);
plot(Nrand_prctiles(2:3), [1 1]*(Ifrac+0.005), '-', 'Color', randcol);
plot(Nord_prctiles(2:3), [1 1]*(Ifrac-0.005), '-', 'Color', avgcol);
plot(Nrand_prctiles(1), Ifrac+0.005, 'o', 'MarkerSize', 4, ...
    'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', randcol);
plot(Nord_prctiles(1), Ifrac-0.005, 'o', 'MarkerSize', 4, ...
    'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', avgcol);
xlim([0 N]);
ylim([0 1.05]);
set(gca,'Box','off','XTick',[0 100 200 300],'YColor','none');
xlabel('Number of neurons N');
ylabel('Fraction of total population information [%]');
legend('boxoff');
legend([p1 p2 p3 p4], ...
    {'45\circ vs. 90\circ', '0\circ vs. 45\circ', 'best avg', 'random'},...
    'Location','southeast');


%% compute for each discrimination the pop size at which Ifrac is reached
Nrand = zeros(3, 0);
Nord = zeros(3, 0);
dostats = false(0, 0);
pdiff = [];
for i = 1:length(datasets)
    fprintf('--- Dataset %s ...\n', datasets{i});
    orderfile = sprintf('%s%s%s_dori%d_c1%s.mat', ...
        orderpath, filesep, datasets{i}, doriid, datasuffix);
    fprintf('Loading %s ...\n', orderfile);
    dord = load(orderfile);
    for j = 1:size(dord.oricombs, 2)
        In = dord.In_rand(:,:,j);
        In_frac = bsxfun(@rdivide, In, In(:,end));
        Nrand_samples = sum(bsxfun(@times, 2:dord.N, diff(In_frac >= Ifrac, [], 2)), 2)';
        % same for information increase
        Nord_samples = sampleNfrac(dord.In(j,:,end), dord.In_var(j,:,end), ...
            Ifrac, size(In,1));
        % add results to statistics
        Nrand = cat(2, Nrand, prctile(Nrand_samples, Nprctiles)');
        Nord = cat(2, Nord, prctile(Nord_samples, Nprctiles)');
        pdiff = cat(2, pdiff, mean(Nord_samples >= Nrand_samples));
        dostats = cat(2, dostats, any(j == doristats));
    end
end
% perform statistical test
fprintf('\nComparing Nrand, Nord with paired t-test\n');
[~,p,~,stats] = ttest(Nrand(1, dostats)' - Nord(1, dostats)');
fprintf('<Nrand - Nord> = %f\n', mean(Nrand(1, dostats)' - Nord(1, dostats)'));
fprintf('t(%d) = %f,  p = %f\n', stats.df, stats.tstat, p);


%% plot results
if singlefig
    subplotcm([10.5 0.8 3 3]);  hold on;
    text(-1,4,'b','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure; hold on;
end
plot([150 400], [150 400], '-', 'Color', [1 1 1]*0.5);
% error bars
for i = 1:length(pdiff)
    plot([Nrand(2,i) Nrand(3,i)], [1 1]*Nord(1,i), '-', ...
        'LineWidth', 0.25, 'Color', [1 1 1]*0.6);
    plot([1 1]*Nrand(1,i), [Nord(2,i) Nord(3,i)], '-', ...
        'LineWidth', 0.25, 'Color', [1 1 1]*0.6);
end
% centers
psig = pdiff < 0.05;
plot(Nrand(1,psig), Nord(1,psig), 'o', ...
    'MarkerSize', 3, 'MarkerFaceColor', [1 1 1] * 0.2, 'MarkerEdgeColor', 'none');
plot(Nrand(1,~psig), Nord(1,~psig), 'o', ...
    'MarkerSize', 3, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [1 1 1] * 0.2);
xlim([150 400]);  ylim([150 400]);
xlabel('population size for 90% information, random order');
ylabel('population size for 90% information, optimized order');


%% save figure
if singlefig
    fprintf('\nWriting figure to figS11.pdf\n');
    print(['figs' filesep 'figS11'], '-dpdf');
end



function Ns = sampleNfrac(In_mu, In_var, Ifrac, Nsamples)
%% returns N samples at which I = Ifrac * max(I)

% draw info scaling samples that satisfy I_mu and I_var
Iincr_mu = diff([0 In_mu(:)']);
Iincr_var = diff([0 In_var(:)']);
Iincr = bsxfun(@plus, Iincr_mu(:)', ...
    bsxfun(@times, sqrt(Iincr_var(:)'), randn(Nsamples, length(Iincr_var))));
In = cumsum(Iincr, 2);
% determine N in each sample where I = Ifrac * max(I)
In_frac = bsxfun(@rdivide, In, In(:,end));
Ns = sum(bsxfun(@times, 2:length(In_mu), diff(In_frac >= Ifrac, [], 2)), 2)';


function o = errorshade(x, ymu, ysd, c)
% plots shaded error regions ymu+/-ysd in color c over x
o = patch([x(:)' fliplr(x(:)')], [(ymu(:)'+ysd(:)') fliplr(ymu(:)'-ysd(:)')],1,...
    'FaceColor',c,'EdgeColor','none');
