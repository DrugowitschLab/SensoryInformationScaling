function plotFigS16
%% plots statistics of the 


%% (general) settings
singlefig = true;  % format for paper if true
covpath = ['.' filesep 'cov_cache'];
N = 300;
Ms = [500 1000 5000 10000 15000];
%Mvar = 500;
Mvar = 10000;
Mcols = [27 158 119; 217 95 2; 117 112 179; 231 41 138; 102 166 30] / 255;


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


%% load data
varfrac = NaN(length(Ms), N);
for i = 1:length(Ms)
    covfile = sprintf('%s%scovsim_N%d.mat', covpath, filesep, Ms(i));
    fprintf('Loading %s ...\n', covfile);
    d = load(covfile);
    assert(d.N == N);
    avgcov = mean(d.shufcov(:,:,:), 3);
    diagvar = cumsum(diag(avgcov)');
    fullvar = arrayfun(@(j) sum(sum(avgcov(1:j,1:j))), 1:N);
    varfrac(i,:) = (diagvar-fullvar)./fullvar;

    % store variances for example M
    if Ms(i) == Mvar
        shufvar = mean(d.shufvar, 1);
        momshufvar = mean(d.momshufvar, 1);
        subvar = mean(d.subvar, 1);
        momsubvar = mean(d.momsubvar, 1);
        totalvar = mean(d.totalvar, 1);
    end
end


%% plot examples variance due to shuffling, and moment estimates
if singlefig
    figure('Name', 'Figure', 'Units', 'centimeters', 'Position', [0 0 11.5 5.1]);
    subplotcm([1 0.8 4.5 3.5]);  hold on;
    text(-1,4,'a','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
plot(1:N, shufvar, '-', 'LineWidth', 1, 'Color', [55 126 184]/255);
plot(1:N, momshufvar, '-', 'LineWidth', 1, 'Color', [228 26 28]/255);
%plot(1:N, subvar, '--', 'LineWidth', 1, 'Color', [55 126 184]/255);
%plot(1:N, momsubvar, '--', 'LineWidth', 1, 'Color', [228 26 28]/255);
plot(1:N, totalvar, '-', 'LineWidth', 1, 'Color', [1 1 1]*0.1);
plot(1:N, shufvar + momshufvar, '--', 'Color', [1 1 1]*0.5);
set(gca, 'Box', 'off');
xlabel('population size n');     xlim([1 N]);
ylabel('var(\Delta I_n)'); ylim([0 550]);
legend('boxoff');
legend({'shuffling','moment estimates','total','shuffling + moment estimates'},...
    'Location','northwest');


%% plot variance overestimate
if singlefig
    subplotcm([6.5 0.8 4.5 3.5]);  hold on;
    text(-1,4,'b','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
for i = 1:length(Ms)
    plot(1:N, varfrac(i,:), '-', 'LineWidth', 1, 'Color', Mcols(i,:));
end
set(gca, 'Box', 'off');
xlim([1 N]);  ylim([0 25]);
xlabel('population size n');
ylabel('fraction variance overestimate');
Mlegend = arrayfun(@(i) sprintf('%d', Ms(i)), 1:length(Ms), ...
    'UniformOutput', false);
legend('boxoff');
legend(Mlegend,'Location','northwest');


%% write figure to file
if singlefig
    fprintf('\nWriting figure to figS16.pdf\n');
    print(['figs' filesep 'figS16'], '-dpdf');
end
