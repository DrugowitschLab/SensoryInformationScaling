function plotFig6
%% plots panels of figure 6

%% (general) settings
singlefig = true;   % set to true to have all planels in single figure
momentpath = ['.' filesep 'moment_cache'];
fitpath = ['.' filesep 'fits'];

momentfile = 'aj60a_o1-2';
datasets = {'aj60a','aj60b','aj60c','aj60d','aj61a','aj61b','aj61c'};
d45discrs = {'o1-2','o2-3','o3-4','o4-5','o5-6','o6-7','o7-8','o1-8'};
d45dirscnonoverlap = [1 3 5 7];
dsmouse = [5 5 5 5 6 6 6];
dssession = [1 2 3 4 1 2 3];
dslabels = {'5A','5B','5C','5D','6A','6B','6C'};
obspopcol = [0.95 0.95 0.99];
cprcriles = [50 25 75];
% order of datasets needs to be the same as in 'datasets' cell array
sessfiles = {...
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
      {'aj61c_dori1a_c2.mat', 'aj61c_dori2a_c2.mat', 'aj61c_dori3a_c2.mat'}}} ...
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
    figure('Name', 'Figure', 'Units', 'centimeters', 'Position', [0 0 9.5 8.3]);
end


%% info scaling for low vs. high contrast
fprintf('Loading %s ...\n', [momentpath filesep momentfile '_c1.mat']);
dlo = load([momentpath filesep momentfile '_c1.mat']);
fprintf('Loading %s ...\n', [momentpath filesep momentfile '_c2.mat']);
dhi = load([momentpath filesep momentfile '_c2.mat']);
Iincr_mu_lo = mean(dlo.Iincr_samples, 1);
Iincr_var_lo = var(dlo.Iincr_samples, [], 1);
Iincr_mu_hi = mean(dhi.Iincr_samples, 1);
Iincr_var_hi = var(dhi.Iincr_samples, [], 1);
I_mu_lo = [0 cumsum(Iincr_mu_lo)];
I_sd_lo = [0 sqrt(cumsum(Iincr_var_lo))];
I_mu_hi = [0 cumsum(Iincr_mu_hi)];
I_sd_hi = [0 sqrt(cumsum(Iincr_var_hi))];
% plot scaling
if singlefig
    subplotcm([1 4.8 4 3]);  hold on;
    text(-0.8,3.5,'a','Units','centimeters','FontWeight','bold',...
         'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
alpha(errorshade(0:(length(I_mu_lo)-1), I_mu_lo, I_sd_lo, plotcols(5,[],1)), 0.2);
alpha(errorshade(0:(length(I_mu_hi)-1), I_mu_hi, I_sd_hi, plotcols(5,[],2)), 0.2);
plot(0:(length(I_mu_lo)-1), I_mu_lo, '-', 'Color', plotcols(5,[],1), 'LineWidth', 1);
plot(0:(length(I_mu_hi)-1), I_mu_hi, '-', 'Color', plotcols(5,[],2), 'LineWidth', 1);
xlim([0 length(I_mu_lo)]);  ylim([0 50]);
set(gca,'Box','off','XTick',[0 100 200 300],'YTick',[0 10 20 30 40 50]);
xlabel('number of neurons N');
ylabel('Fisher information [rad^{-2}]');


%% summary across all dtheta=45 comparisons
% collect info measures and perform stats comparison
Ilo = NaN(length(datasets), length(d45discrs), 2); % 3rd dim [mu var]
Ihi = NaN(length(datasets), length(d45discrs), 2);
lohip = NaN(length(datasets), length(d45discrs));
for i = 1:length(datasets)
    for j = 1:length(d45discrs)
        filebase = [momentpath filesep datasets{i} '_' d45discrs{j}];
        fprintf('Loading %s ...\n', [filebase '_c1.mat']);
        dlo = load([filebase '_c1.mat']);
        fprintf('Loading %s ...\n', [filebase '_c2.mat']);
        dhi = load([filebase '_c2.mat']);
        [Ilo(i,j,1), Ilo(i,j,2)] = fullpopinfo(dlo);
        [Ihi(i,j,1), Ihi(i,j,2)] = fullpopinfo(dhi);
        lohip(i,j) = normcdf(0, Ihi(i,j,1)-Ilo(i,j,1), ...
            sqrt(Ilo(i,j,2)+Ihi(i,j,2)));
        fprintf('p( Ihi <= Ilo ) = %f', lohip(i,j));
        if lohip(i,j) >= 0.05, fprintf(' (n.s.)\n');
        else, fprintf('\n'); end
    end
end
% datasets for which >= half of discrimination has non-sig increase
nonsigdatasets = find(sum(lohip > 0.05, 2) >= (size(lohip,2)/2));
fprintf('datasets for which >= half of discrs. have non-signif. increase:\n');
if isempty(nonsigdatasets)
    fprintf('none\n');
else
    for i = 1:(length(nonsigdatasets)-1)
        fprintf('%s, ', datasets{nonsigdatasets(i)});
    end
    fprintf('%s\n', datasets{nonsigdatasets(end)});
end
% stats comparison
fprintf('Performing paired (non-overlapping) t-test, Ilo vs. Ihi\n');
fprintf('<Ilo> = %f, <Ihi> = %f\n', ...
    mean(reshape(Ilo(:,d45dirscnonoverlap,1),1,[])), ...
    mean(reshape(Ihi(:,d45dirscnonoverlap,1),1,[])));
[~,p,~,stats] = ttest(reshape(Ihi(:,d45dirscnonoverlap,1) - ...
                              Ilo(:,d45dirscnonoverlap,1),[],1));
fprintf('t(%d) = %f, p = %f\n\n', stats.df, stats.tstat, p);
% plot data
if singlefig
    subplotcm([6 4.8 3 3]);  hold on;
    text(-0.8,3.5,'b','Units','centimeters','FontWeight','bold',...
      'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure; hold on;
end
plot([0 110], [0 110], '-', 'Color', [1 1 1]*0.5);
% plot SDs
for i = 1:length(datasets)
    for j = 1:length(d45discrs)
        c = plotcols(dsmouse(i), dssession(i))*0.5+0.5;
        plot(Ilo(i,j,1)+[-1 1]*Ilo(i,j,2), [1 1]*Ihi(i,j,1), '-', ...
            'Color', c, 'LineWidth', 0.25);
        plot([1 1]*Ilo(i,j,1), Ihi(i,j,1)+[-1 1]*Ihi(i,j,2), '-', ...
            'Color', c, 'LineWidth', 0.25);
    end
end
% plot means
plegend = NaN(1, length(datasets));
for i = 1:length(datasets)
    j = lohip(i,:) >= 0.05;
    plot(Ilo(i,j,1), Ihi(i,j,1), 'o', 'MarkerSize', 3, ...
        'MarkerFaceColor', [1 1 1], ...
        'MarkerEdgeColor', plotcols(dsmouse(i), dssession(i)));
    plegend(i) = plot(Ilo(i,~j,1), Ihi(i,~j,1), 'o', 'MarkerSize', 3, ...
        'MarkerFaceColor', plotcols(dsmouse(i), dssession(i)), ...
        'MarkerEdgeColor', 'none');    
end
xlim([0 110]);  ylim([0 110]);
xlabel('Fisher information, low contrast [rad^{-2}]');
ylabel('Fisher information, high contrast, [rad^{-2}]');
legend('boxoff');
legend(plegend, dslabels,'Location','northwest');
set(gca,'Box','off',...
    'XTick',[0 50 100],'YTick',[0 50 100]);


%% schematic about two possible ways that info can grow
Ns = 0:500;
Nobs = 30;
Iinfs = [0.3 0.45];  % info scalings
cs = [0.01 0.02];
if singlefig
    subplotcm([1 3.2 1.75 0.8]);  hold on;
    text(-1,1.3,'c','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;
    subplot('Position', [0.13 0.525 0.775 0.4]);  hold on;
end
patch([0 Nobs Nobs 0], [0 0 0.55 0.55], 1, ...
    'FaceColor', obspopcol, 'EdgeColor', 'none');
plot(Ns, 1 ./ (1 ./ (cs(1) * Ns) + 1 ./ Iinfs(1)), 'LineWidth', 1, 'Color', plotcols(5,[],1));
plot(Ns, 1 ./ (1 ./ (cs(2) * Ns) + 1 ./ Iinfs(2)), 'LineWidth', 1, 'Color', plotcols(5,[],2));
plot([0 max(Ns)], Iinfs(1)*[1 1], '--', 'Color', plotcols(5,[],1));
plot([0 max(Ns)], Iinfs(2)*[1 1], '--', 'Color', plotcols(5,[],2));
text(max(Ns), Iinfs(1), 'I_\infty', 'Color', plotcols(5,[],1), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(max(Ns), Iinfs(2), 'I_\infty', 'Color', plotcols(5,[],2), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
xlim([0 max(Ns)]);  ylim([0 0.55]);
set(gca,'Box','off', 'XTick', [], 'YTick', 0);
xlabel('Number of neurons N');
ylabel('Fisher information');

if singlefig
    subplotcm([3.25 3.2 1.75 0.8]);  hold on;
else
    subplot('Position', [0.13 0.11 0.775 0.4]);  hold on;
end
patch([0 Nobs Nobs 0], [0 0 0.55 0.55], 1, ...
    'FaceColor', obspopcol, 'EdgeColor', 'none');
plot(Ns, 1 ./ (1 ./ (cs(1) * Ns) + 1 ./ Iinfs(1)), 'LineWidth', 1, 'Color', plotcols(5,[],1));
plot(Ns, 1 ./ (1 ./ (cs(2) * Ns) + 1 ./ Iinfs(1)), 'LineWidth', 1, 'Color', plotcols(5,[],2));
plot([0 max(Ns)], Iinfs(1)*[1 1], '--', 'Color', [1 1 1]*0.5);
text(max(Ns), Iinfs(1), ...
    '{\color[rgb]{.8,0,0}I_\infty} = {\color[rgb]{0,.8,0}I_\infty}', ...
    'Color', [0 0 0], ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
xlim([0 max(Ns)]);  ylim([0 0.55]);
set(gca,'Box','off', 'XTick', 0, 'YTick', 0);


%% estimated Iinf for different datasets, dthetas
unreswarning = warning('query', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
pw = (4 - 0.5*(animals-1))/animals;
ph = 2;
totalsessions = sum(arrayfun(@(i) length(sessfiles{i}{1}), 1:animals));
assert(totalsessions == length(datasets));
Iinfs = NaN(totalsessions, 2, 3);
si = 1;
for i = 1:animals
    if singlefig
        subplotcm([(1+(i-1)*(pw+0.5)) 0.8 pw ph]);  hold on;
        if i == 1
            text(-0.8,2.5,'d','Units','centimeters','FontWeight','bold',...
                'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
        end
    else
        figure;  hold on;
    end
    % differentiate between single- and multi-contrast animals
    sessions = length(sessfiles{i}{1});
    sx = bsxfun(@plus, repmat(0.04*(0:(sessions-1)), 2, 1), ...
        [0.01; -0.01]-0.02*(sessions-1));
    % iterate over contrasts
    for j = 1:2
        for k = 1:sessions
            meds = NaN(1, length(sessdorinames));
            for l = 1:length(sessdorinames)
                fprintf('Loading %s ...\n', [fitpath filesep sessfiles{i}{j}{k}{l}]);
                m = load([fitpath filesep sessfiles{i}{j}{k}{l}]);
                sms = m.fits{2};
                assert(strcmp(sms.name, 'limlin'));  % make sure to get the right model
                assert(strcmp(sms.pnames{1}, 'Iinf'));
                assert(strcmp(sms.pnames{2}, 'c'));
                Iinfss = reshape(sms.mc.ss(:,1,:),1,[]);
                meds(l) = prctile(Iinfss, 50);
                plotPrctiles(l+sx(j,k), Iinfss, plotcols(i+4,k,j));
                % store for later, for dtheta = 45
                if l == 1
                    Iinfs(si, j, :) = prctile(Iinfss, cprcriles);
                    si = si + 1;
                    if j == 1 && k == sessions, si = si - sessions; end
                end
            end
            plot((1:length(sessdorinames))+sx(j,k), meds, '-', 'Color', plotcols(i+4,k,j));
        end
    end
    % plot formatting
    set(gca,'Box','off','YScale','log',...
        'XTick',1:length(sessdorinames),'XTickLabel',sessdorinames,....
        'YTick',[1:9 10:10:90 100:100:1000],'YTickLabel',...
        {'10^0','','','','','','','','',...
         '10^1','','','','','','','','',...
         '10^2','','','','','','','','','10^3'});
    if i == 1
        xlabel('drift direction difference');
        ylabel('Asymptotic information I_\infty [rad^{-2}]');
    else
        set(gca, 'YTickLabel',[]);
    end
    xlim([0.75 length(sessdorinames)+0.2]);  ylim([1 1000]);
    text(length(sessdorinames), 1, sprintf('mouse %d', i+4), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
warning(unreswarning.state, 'MATLAB:dispatcher:UnresolvedFunctionHandle');


%% plot hi/lo info against each other
fprintf('\n_Performing one-tailed Wilcoxon signed rank test, Iinf(lo) vs. Iinf(hi)\n');
fprintf('<Iinf(lo)> = %f, <Iinf(hi)> = %f\n', ...
    mean(Iinfs(:,1,1)), mean(Iinfs(:,2,1)));
[p,~,stats] = signrank(Iinfs(:,2,1),Iinfs(:,1,1), 'tail', 'right');
fprintf('signedrank W+ = %f, p = %f\n\n', stats.signedrank, p); 
if singlefig
    subplotcm([6 0.8 3 3]);  hold on;
    text(-1.2,3.5,'e','Units','centimeters','FontWeight','bold',...
         'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure; hold on;
end
% diagonal
plot([20 400], [20 400], '-', 'Color', [1 1 1]*0.5);
% plot error bars
for i = 1:totalsessions
    plot([Iinfs(i,1,2) Iinfs(i,1,3)], [1 1]*Iinfs(i,2,1), ...
        '-', 'Color', plotcols(dsmouse(i),dssession(i))*0.5+0.5, ...
        'LineWidth', 0.25);
    plot([1 1]*Iinfs(i,1,1), [Iinfs(i,2,2) Iinfs(i,2,3)], ...
        '-', 'Color', plotcols(dsmouse(i),dssession(i))*0.5+0.5, ...
        'LineWidth', 0.25);
end
% then means
for i = 1:totalsessions
    if any(nonsigdatasets == i), mfc = [1 1 1]; meg = plotcols(dsmouse(i),dssession(i));
    else, mfc = plotcols(dsmouse(i),dssession(i)); meg = 'none'; end
    plot(Iinfs(i,1,1), Iinfs(i,2,1), 'o', ...
        'MarkerSize', 3, 'MarkerFaceColor', mfc, 'MarkerEdgeColor', meg);
end
xlim([20 400]);  ylim([20 400]);
xlabel('Asymptotic information {\color[rgb]{0.8,0,0}I_\infty}, low contrast [rad^{-2}]');
ylabel('Asymptotic information {\color[rgb]{0,0.8,0}I_\infty}, high contrast [rad^{-2}]');
set(gca,'Box','off','XScale','log','YScale','log',...
    'XTick',[20:10:100 200 300 400],'XTickLabel',{'20','','','','','','','','100','200','300',''},...
    'YTick',[20:10:100 200 300 400],'YTickLabel',{'20','','','','','','','','100','200','300',''});


%% save figure
if singlefig
    fprintf('\nWriting figure to fig6.pdf\n');
    print(['figs' filesep 'fig6'], '-dpdf');
end



function o = errorshade(x, ymu, ysd, c)
% plots shaded error regions ymu+/-ysd in color c over x
o = patch([x(:)' fliplr(x(:)')], [(ymu(:)'+ysd(:)') fliplr(ymu(:)'-ysd(:)')],1,...
    'FaceColor',c,'EdgeColor','none');


function [mu, var] = fullpopinfo(d)
%% returns mean / variance of info in full population for given data d
T = d.T;
N = size(d.Iincr_samples, 2);
ds = d.ds;
mu = sum(d.Iincr_samples(1,:));
var = 2*mu^2 / (2*T-N-3) * (1 + 4*(2*T-3) / (T*mu*d.ds^2)+ ...
    4*N*(2*T-3) / (T^2*mu^2*ds^4)); 
