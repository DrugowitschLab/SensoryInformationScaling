function plotFigS9
%% plots SI figure showing model comparison for pooled vs. indiv. fits


%% (general) settings
singlefig = true;  % format for paper if true
datasets = {'m25a', 'm25b', 'm26a', 'm26b', ...
    'aj42a', 'aj42b', 'aj42c', 'aj42d', 'aj42e', ...
    'aj43a', 'aj43b', 'aj43c', 'aj43d', 'aj43e', 'aj43f', 'aj43g'};
fitpath = ['.' filesep 'fits'];
pooledfiles = {'dori1a_c1', 'dori2a_c1', 'dori3a_c1'};
poolednames = {'45\circ','90\circ','135\circ'};
indivfiles = {{'o1-2_c1', 'o3-4_c1', 'o5-6_c1', 'o7-8_c1'}, ...
    {'o1-3_c1','o2-4_c1','o5-7_c1','o6-8_c1'}, ...
    {'o1-4_c1','o2-7_c1','o5-8_c1'}};
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


%% collect statistics across datasets and pooled fits
waics = NaN(length(datasets), length(pooledfiles), 2); % [pooled, indv]
for i = 1:length(datasets)
    for j = 1:length(pooledfiles)
        % pooled fits
        fitfile = [fitpath filesep datasets{i} '_' pooledfiles{j} '.mat'];
        fprintf('Loading %s ...\n', fitfile);
        m = load(fitfile);
        assert(strcmp(m.fits{2}.name, 'limlin'));
        waics(i,j,1) = m.fits{2}.mc.(waic_field);
        % inidividual fits
        waicindv = 0;
        for k = 1:length(indivfiles{j})
            fitfile = [fitpath filesep datasets{i} '_' indivfiles{j}{k} '_norm.mat'];
            fprintf('Loading %s ...\n', fitfile);
            m = load(fitfile);
            assert(strcmp(m.sms{2}.name, 'limlin'));
            waicindv = waicindv + m.sms{2}.mc.(waic_field);
        end
        waics(i,j,2) = waicindv;
    end
end


%% plot statistics
dx = 0.02*(0:(length(datasets)-1)) - 0.01*(length(datasets)-1);
dm = waics(:,:,2) - waics(:,:,1);
if singlefig
    figure('Name', 'Figure', 'Units', 'centimeters', 'Position', [0 0 5.1 5]);
    subplotcm([1 0.8 3.5 3.5]);  hold on;
else
    figure;  hold on;
end    
for i = 1:length(datasets)
    plot((1:length(pooledfiles)) + dx(i), dm(i,:), '-', 'Color', [1 1 1]*0.3);
end
for i = 1:length(pooledfiles)
    plot(i + dx, dm(:,i), 'o', 'MarkerSize', 4, ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0 0 0]);
end
plot([0.75 length(poolednames)+0.25], [0 0], '-', 'Color', [1 1 1]*0.5);
xlabel('drift direction difference');
ylabel('WAIC_{indv} - WAIC_{pooled}');
xlim([0.75 length(poolednames)+0.25]);  ylim([-2 8]);
set(gca,'Box','off','XTick',1:length(poolednames),'XTickLabel',poolednames);


%% write figure to file
if singlefig
    fprintf('\nWriting figure to figS9.pdf\n');
    print(['figs' filesep 'figS9'], '-dpdf');
end
