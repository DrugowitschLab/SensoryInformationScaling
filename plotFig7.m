function plotFig7
%% plots panels of figure 7


%% (general) settings
singlefig = true;   % true for paper fig, with all subplots in one figure
fitpath = ['.' filesep 'fits'];


afrac = 0.95;
obspopcol = [0.95 0.95 0.99];
plotprctiles = [5 25 50 75 95];
cprcriles = [50 25 75];
sessfiles = {...
    { {'m25a_dori1a_c1.mat',  'm25a_dori2a_c1.mat',  'm25a_dori3a_c1.mat'}, ...  % m25
      {'m25b_dori1a_c1.mat',  'm25b_dori2a_c1.mat',  'm25b_dori3a_c1.mat'}}, ...
    { {'m26a_dori1a_c1.mat',  'm26a_dori2a_c1.mat',  'm26a_dori3a_c1.mat'}, ...  % m26
      {'m26b_dori1a_c1.mat',  'm26b_dori2a_c1.mat',  'm26b_dori3a_c1.mat'}}, ...
    { {'aj42a_dori1a_c1.mat', 'aj42a_dori2a_c1.mat', 'aj42a_dori3a_c1.mat'}, ... % aj42
      {'aj42b_dori1a_c1.mat', 'aj42b_dori2a_c1.mat', 'aj42b_dori3a_c1.mat'}, ...
      {'aj42c_dori1a_c1.mat', 'aj42c_dori2a_c1.mat', 'aj42c_dori3a_c1.mat'}, ...
      {'aj42d_dori1a_c1.mat', 'aj42d_dori2a_c1.mat', 'aj42d_dori3a_c1.mat'}, ...
      {'aj42e_dori1a_c1.mat', 'aj42e_dori2a_c1.mat', 'aj42e_dori3a_c1.mat'}}, ...
    { {'aj43a_dori1a_c1.mat', 'aj43a_dori2a_c1.mat', 'aj43a_dori3a_c1.mat'}, ... % aj43
      {'aj43b_dori1a_c1.mat', 'aj43b_dori2a_c1.mat', 'aj43b_dori3a_c1.mat'}, ...
      {'aj43c_dori1a_c1.mat', 'aj43c_dori2a_c1.mat', 'aj43c_dori3a_c1.mat'}, ...
      {'aj43d_dori1a_c1.mat', 'aj43d_dori2a_c1.mat', 'aj43d_dori3a_c1.mat'}, ...
      {'aj43e_dori1a_c1.mat', 'aj43e_dori2a_c1.mat', 'aj43e_dori3a_c1.mat'}, ...
      {'aj43f_dori1a_c1.mat', 'aj43f_dori2a_c1.mat', 'aj43f_dori3a_c1.mat'}, ...
      {'aj43g_dori1a_c1.mat', 'aj43g_dori2a_c1.mat', 'aj43g_dori3a_c1.mat'}}, ...
    {{{'aj60a_dori1a_c1.mat', 'aj60a_dori2a_c1.mat', 'aj60a_dori3a_c1.mat'}, ... % aj60
      {'aj60b_dori1a_c1.mat', 'aj60b_dori2a_c1.mat', 'aj60b_dori3a_c1.mat'}, ...
      {'aj60c_dori1a_c1.mat', 'aj60c_dori2a_c1.mat', 'aj60c_dori3a_c1.mat'}, ...
      {'aj60d_dori1a_c1.mat', 'aj60d_dori2a_c1.mat', 'aj60d_dori3a_c1.mat'}}, ...
     {{'aj60a_dori1a_c2.mat', 'aj60a_dori2a_c2.mat', 'aj60a_dori3a_c2.mat'}, ...
      {'aj60b_dori1a_c2.mat', 'aj60b_dori2a_c2.mat', 'aj60b_dori3a_c2.mat'}, ...
      {'aj60c_dori1a_c2.mat', 'aj60c_dori2a_c2.mat', 'aj60c_dori3a_c2.mat'}, ...
      {'aj60d_dori1a_c2.mat', 'aj60d_dori2a_c2.mat', 'aj60d_dori3a_c2.mat'}}}, ...
    {{{'aj61a_dori1a_c1.mat', 'aj61a_dori2a_c1.mat', 'aj61a_dori3a_c1.mat'}, ... % aj61
      {'aj61b_dori1a_c1.mat', 'aj61b_dori2a_c1.mat', 'aj61b_dori3a_c1.mat'}, ...
      {'aj61c_dori1a_c1.mat', 'aj61c_dori2a_c1.mat', 'aj61c_dori3a_c1.mat'}}, ...
     {{'aj61a_dori1a_c2.mat', 'aj61a_dori2a_c2.mat', 'aj61a_dori3a_c2.mat'}, ...
      {'aj61b_dori1a_c2.mat', 'aj61b_dori2a_c2.mat', 'aj61b_dori3a_c2.mat'}, ...
      {'aj61c_dori1a_c2.mat', 'aj61c_dori2a_c2.mat', 'aj61c_dori3a_c2.mat'}}}, ...
     };
animals = length(sessfiles);
sessdorinames = {'45\circ','90\circ','135\circ'};


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
    figure('Name', 'Figure', 'Units', 'centimeters', 'Position', [0 0 16 7.5]);
end


%% N95 schematic
Ns = 0:500;
N95 = round(max(Ns) * 0.8);  % N95 at 80% of max(N)
Iinf = 0.3;
Imax = 0.5;
c = 0.95 * Iinf / (0.05 * N95); % required c for desired N95
In = 1 ./ (1 ./ (c * Ns) + 1 ./ Iinf);
if singlefig
    subplotcm([1 4 4.5 3]);  hold on;
    text(-1,3.5,'a','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
plot(Ns, In, '-', 'LineWidth', 1, 'Color', [0 0 0.8]);
plot([1 1]*N95, [0 0.95*Iinf], '-', 'Color', [0.5 0.5 0.5]);
text(N95, 0.5*0.95*Iinf, '95% I_\infty', ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
text(N95, 0, 'N_{95}', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
plot([min(Ns) max(Ns)], Iinf*[1 1], 'k--');
xlim([min(Ns) max(Ns)]);
ylim([0 Imax]);
xlabel('number of neurons N');
ylabel('Fisher information');
set(gca, 'Box', 'off', 'XTick', 0, 'YTick', 0);
text(0.99*max(Ns), Iinf, 'Asumptotic information I_\infty', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');


%% schema for different / same N95 when info increases
Ns = 0:500;
Nobs = 30;
N95s = round(max(Ns) * [0.35 0.9]); % desired N95s
Iinfs = [0.25 0.48];  % info scaling
if singlefig
    subplotcm([7 5.6 2.5 1.4]);  hold on;
    text(-1,1.9,'c','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;
    subplot('Position', [0.13 0.525 0.775 0.4]);  hold on;
end
cs = 0.95 * Iinfs ./ (0.05 * N95s(1));
patch([0 Nobs Nobs 0], [0 0 0.55 0.55], 1, ...
    'FaceColor', obspopcol, 'EdgeColor', 'none');
plot(Ns, 1 ./ (1 ./ (cs(1) * Ns) + 1 ./ Iinfs(1)), 'LineWidth', 1, 'Color', plotcols(5,[],1));
plot(Ns, 1 ./ (1 ./ (cs(2) * Ns) + 1 ./ Iinfs(2)), 'LineWidth', 1, 'Color', plotcols(5,[],2));
plot([0 max(Ns)], Iinfs(1)*[1 1], '--', 'Color', plotcols(5,[],1));
plot([0 max(Ns)], Iinfs(2)*[1 1], '--', 'Color', plotcols(5,[],2));
plot([1 1]*1.005*N95s(1), [0 0.95*Iinfs(2)], '-', 'Color', plotcols(5,[],1));
plot([1 1]*0.995*N95s(1), [0 0.95*Iinfs(1)], '-', 'Color', plotcols(5,[],2));
text(N95s(1)*1.005, 0.5*0.95*Iinfs(2), '95% I_\infty', 'Color', plotcols(5,[],2), ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
text(N95s(1)*0.995, 0.5*0.95*Iinfs(1), '95% I_\infty', 'Color', plotcols(5,[],1), ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
text(N95s(1)*1.005, 0, 'N_{95}', 'Color', plotcols(5,[],2), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
text(N95s(1)*0.995, 0, 'N_{95}', 'Color', plotcols(5,[],1), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
xlim([0 max(Ns)]);  ylim([0 0.55]);
set(gca,'Box','off', 'XTick', [], 'YTick', 0);

if singlefig
    subplotcm([7 4 2.5 1.4]);  hold on;
else
    subplot('Position', [0.13 0.11 0.775 0.4]);  hold on;
end
cs = 0.95 * Iinfs ./ (0.05 * N95s);
patch([0 Nobs Nobs 0], [0 0 0.55 0.55], 1, ...
    'FaceColor', obspopcol, 'EdgeColor', 'none');
plot(Ns, 1 ./ (1 ./ (cs(1) * Ns) + 1 ./ Iinfs(1)), 'LineWidth', 1, 'Color', plotcols(5,[],1));
plot(Ns, 1 ./ (1 ./ (cs(2) * Ns) + 1 ./ Iinfs(2)), 'LineWidth', 1, 'Color', plotcols(5,[],2));
plot([0 max(Ns)], Iinfs(1)*[1 1], '--', 'Color', plotcols(5,[],1));
plot([0 max(Ns)], Iinfs(2)*[1 1], '--', 'Color', plotcols(5,[],2));
plot([1 1]*N95s(2), [0 0.95*Iinfs(2)], '-', 'Color', plotcols(5,[],2));
plot([1 1]*N95s(1), [0 0.95*Iinfs(1)], '-', 'Color', plotcols(5,[],1));
text(N95s(2), 0.5*0.95*Iinfs(2), '95% I_\infty', 'Color', plotcols(5,[],2), ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
text(N95s(1), 0.5*0.95*Iinfs(1), '95% I_\infty', 'Color', plotcols(5,[],1), ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
text(N95s(2), 0, 'N_{95}', 'Color', plotcols(5,[],2), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
text(N95s(1), 0, 'N_{95}', 'Color', plotcols(5,[],1), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
xlim([0 max(Ns)]);  ylim([0 0.55]);
set(gca,'Box','off', 'XTick', 0, 'YTick', 0);
xlabel('Number of neurons N');
ylabel('Fisher information');


%% estimated N95 per dataset
unreswarning = warning('query', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
pw = (14.5 - 0.5*(animals-1))/animals;
ph = 2.2;
N95s = zeros(0, 3);
Iinfs = zeros(0, 3);
ccons = zeros(0, 1);
dataids = zeros(0, 3);  % [animal session contrast]
cs = zeros(0, 3);
dthetas = zeros(0, 1);
for i = 1:animals
    if singlefig
        subplotcm([(1+(i-1)*(pw+0.5)) 0.8 pw ph]);  hold on;
        if i == 1
            text(-1,2.7,'b','Units','centimeters','FontWeight','bold',...
                'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
        end
    else
        figure;  hold on;
    end
    % differentiate between single- and multi-contrast animals
    multicon = ~isempty(sessfiles{i}) && iscell(sessfiles{i}{1}{1});  % multiple contrasts?
    if multicon
        sessions = length(sessfiles{i}{1});
        cons = 2;
        asessfiles = sessfiles{i};
        sx = bsxfun(@plus, repmat(0.04*(0:(sessions-1)), 2, 1), ...
            [0.01; -0.01]-0.02*(sessions-1));
    else
        sessions = length(sessfiles{i});
        cons = 1;
        asessfiles = { sessfiles{i} };
        sx = 0.04*(0:(sessions-1)) - 0.02*(sessions-1);
    end
    % iterate over contrasts
    for j = 1:cons
        for k = 1:sessions
            meds = NaN(1, length(sessdorinames));
            for l = 1:length(sessdorinames)
                fprintf('Loading %s ...\n', [fitpath filesep asessfiles{j}{k}{l}]);
                m = load([fitpath filesep asessfiles{j}{k}{l}]);
                sms = m.fits{2};
                assert(strcmp(sms.name, 'limlin'));  % make sure to get the right model
                assert(strcmp(sms.pnames{1}, 'Iinf'));
                assert(strcmp(sms.pnames{2}, 'c'));
                Iinfss = reshape(sms.mc.ss(:,1,:),1,[]);
                css = reshape(sms.mc.ss(:,2,:),1,[]);
                N95ss = (afrac/(1-afrac)) * Iinfss ./ css;
                meds(l) = prctile(N95ss, plotprctiles(3));
                plotPrctiles(l+sx(j,k), N95ss, plotcols(i,k,j));
                % store stats for later analysis
                N95s = cat(1, N95s, prctile(N95ss, cprcriles));
                Iinfs = cat(1, Iinfs, prctile(Iinfss, cprcriles));
                ccons = cat(1, ccons, j + multicon);
                dataids = cat(1, dataids, [i k j]);
                dthetas = cat(1, dthetas, l);  % dtheta idx
                cs = cat(1, cs, prctile(css, cprcriles)); 
            end
            plot((1:length(sessdorinames))+sx(j,k), meds, '-', 'Color', plotcols(i,k,j));
        end
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


%% model Iinf vs. N95 relationship, by means of the Iinf vs. c relationship
log10Iinf = log10(Iinfs(:,1));
log10c = log10(cs(:,1));
log10N95 = log10(N95s(:,1));
tbl = table(log10Iinf, log10c, log10N95);
tbl.dtheta = categorical(dthetas);
tbl.hicon = categorical(ccons > 2);
tbl.log10diff = log10c - log10Iinf;
fprintf('\n - Relationship between Iinf and c\n');
mdl = fitlm(tbl, 'log10c ~ log10Iinf');
disp(mdl);
disp(anova(mdl));
fprintf('Testing null-hypothesis that log10Inf coefficient is 1\n');
[pVal,F,df1] = coefTest(mdl, [0 0; 0 1], [0; 1]);
fprintf('p = %f, F(%d) = %f\n', pVal, df1, F);
fprintf('\n - Relationship between Iinf and c, controlling for con & dtheta\n');
mdl = fitlm(tbl, 'log10c ~ dtheta + hicon + log10Iinf');
disp(mdl);
disp(anova(mdl));


%% plot Iinf vs. c
% re-fit initial model without controls
mdl = fitlm(tbl, 'log10c ~ log10Iinf');
b = mdl.Coefficients.Estimate;
if singlefig
    subplotcm([11 4 4.5 3]);  hold on;
    text(-1,3.5,'d','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
% plot log-linear regression lines
plot([1 1000], 10.^(b(1)+b(2).*log10([1 1000])), '-', 'Color', [0.5 0.5 0.5]);
% the error bars
for i = 1:animals
    multicon = ~isempty(sessfiles{i}) && iscell(sessfiles{i}{1}{1});
    if multicon, cons = 2; else, cons = 1; end
    for j = 1:cons
        k = find(dataids(:,1) == i & dataids(:,3) == j);
        c = plotcols(i,[],j);
        for l = 1:length(k)
            plot([Iinfs(k(l),2) Iinfs(k(l),3)], [1 1] * cs(k(l),1), '-', ...
                'LineWidth', 0.25, 'Color', c);
            plot([1 1] * Iinfs(k(l),1), [cs(k(l),2) cs(k(l),3)], '-', ...
                'LineWidth', 0.25, 'Color', c);
        end
    end
end
% then centers
for i = 1:animals
    multicon = ~isempty(sessfiles{i}) && iscell(sessfiles{i}{1}{1});
    if multicon, cons = 2; else, cons = 1; end
    for j = 1:cons
        k = find(dataids(:,1) == i & dataids(:,3) == j);
        plot(Iinfs(k,1), cs(k,1), 'o', 'MarkerSize', 3, ...
            'MarkerFaceColor', plotcols(i,[],j), 'MarkerEdgeColor', 'none');
    end
end
set(gca,'Box','off','XScale','log','YScale','log');
xlabel('asymptotic information I_\infty');  ylabel('non-limiting scaling c');
xlim([1 1000]);  ylim([0.004 0.4]);
text(1000, 0.4, ...
     sprintf('c \\approx %6.4f I_{\\infty}^{%4.2f}', 10^b(1), b(2)), ...
     'HorizontalAlignment','right','VerticalAlignment','top');


%% save figure
if singlefig
    fprintf('\nWriting figure to fig7.pdf\n');
    print(['figs' filesep 'fig7'], '-dpdf');
end
