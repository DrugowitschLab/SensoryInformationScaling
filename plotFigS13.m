function plotFigS13
%% plots N95 scaling factor


%% (general) settings
singlefig = true;   % true for paper fig, with all subplots in one figure


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


%% mapping between Na and N95
as = linspace(0.8, 0.999, 100);
if singlefig
    figure('Units','centimeters','Position',[0 0 4 3.5]);  hold on;
else
    figure;  hold on;
end
plot(as, (0.05 /0.95) * as ./ (1-as), 'k-', 'LineWidth', 1);
plot([min(as) max(as)], [1 1], '-', 'Color', [1 1 1]*0.5);
plot([0.95 0.95], [0.1 100], '-', 'Color', [1 1 1]*0.5);
text(min(as), 100, 'log', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
text(max(as), 0.1, 'log', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
xlabel('fraction of I_\inft');
ylabel('N_{95} scaling factor');
xlim([min(as) max(as)]);  ylim([0.1 100]);
set(gca,'Box','off','YScale','log','XScale','log');


%% save figure
if singlefig
    fprintf('\nWriting figure to figS13.pdf\n');
    print(['figs' filesep 'figS13'], '-dpdf');
end
