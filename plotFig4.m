function plotFig4
%% plots panels of figure 4

%% (general) settings
singlefig = true;   % set to true to have all planels in single figure
momentpath = ['.' filesep 'moment_cache'];
fitpath = ['.' filesep 'fits'];
momentfile = 'm25b_o3-4_c1.mat';
obspopcol = [0.95 0.95 0.99];
doricols = [     0    0.4470    0.7410; ...   % get(gca, 'colororder') + add 
            0.9290    0.6940    0.1250; ...
            0.4940    0.1840    0.5560; ...
            0.4660    0.6740    0.1880; ...
            0.3010    0.7450    0.9330; ...
            0.6350    0.0780    0.1840; ...
            0.5000    0.5000    0.8000; ...
            0.8500    0.3250    0.0980; ...
            0.1500    0.1500    0.1500; ...
            0.3000    0.3000    0.3000];
dorifiles = {'m25b_o1-8_c1_norm.mat' ...
    'm25b_o1-2_c1_norm.mat' ...
    'm25b_o2-3_c1_norm.mat' ...
    'm25b_o3-4_c1_norm.mat' ...
    'm25b_o4-5_c1_norm.mat' ...
    'm25b_o5-6_c1_norm.mat' ...
    'm25b_o6-7_c1_norm.mat' ...
    'm25b_o7-8_c1_norm.mat' ...
    'm25b_dori1a_c1.mat' ...
    'm25b_dori1b_c1.mat'};
dorinames = {'0\circ vs. 45\circ', '45\circ vs. 90\circ', ...
    '90\circ vs. 135\circ', '135\circ vs. 180\circ', ...
    '180\circ vs. 225\circ', '225\circ vs. 270\circ', ...
    '270\circ vs. 315\circ', '315\circ vs. 0\circ', ...
    'pooled 1', 'pooled 2'};
doriexample = 4;
fitfile = dorifiles{doriexample};
plotprctiles = [5 25 50 75 95];
prctilesat = [0.5 0.75 1 0.75 0.5];
orithreshs = [20 15 10 5 4 3 2 1];
sessfiles = {...
    {{'m25a_dori1a_c1.mat',  'm25a_dori2a_c1.mat',  'm25a_dori3a_c1.mat'}, ...  % m25
     {'m25b_dori1a_c1.mat',  'm25b_dori2a_c1.mat',  'm25b_dori3a_c1.mat'}}, ...
    {{'m26a_dori1a_c1.mat',  'm26a_dori2a_c1.mat',  'm26a_dori3a_c1.mat'}, ...  % m26
     {'m26b_dori1a_c1.mat',  'm26b_dori2a_c1.mat',  'm26b_dori3a_c1.mat'}}, ...
    {{'aj42a_dori1a_c1.mat', 'aj42a_dori2a_c1.mat', 'aj42a_dori3a_c1.mat'}, ... % aj42
     {'aj42b_dori1a_c1.mat', 'aj42b_dori2a_c1.mat', 'aj42b_dori3a_c1.mat'}, ...
     {'aj42c_dori1a_c1.mat', 'aj42c_dori2a_c1.mat', 'aj42c_dori3a_c1.mat'}, ...
     {'aj42d_dori1a_c1.mat', 'aj42d_dori2a_c1.mat', 'aj42d_dori3a_c1.mat'}, ...
     {'aj42e_dori1a_c1.mat', 'aj42e_dori2a_c1.mat', 'aj42e_dori3a_c1.mat'}}, ...
    {{'aj43a_dori1a_c1.mat', 'aj43a_dori2a_c1.mat', 'aj43a_dori3a_c1.mat'}, ... % aj43
     {'aj43b_dori1a_c1.mat', 'aj43b_dori2a_c1.mat', 'aj43b_dori3a_c1.mat'}, ...
     {'aj43c_dori1a_c1.mat', 'aj43c_dori2a_c1.mat', 'aj43c_dori3a_c1.mat'}, ...
     {'aj43d_dori1a_c1.mat', 'aj43d_dori2a_c1.mat', 'aj43d_dori3a_c1.mat'}, ...
     {'aj43e_dori1a_c1.mat', 'aj43e_dori2a_c1.mat', 'aj43e_dori3a_c1.mat'}, ...
     {'aj43f_dori1a_c1.mat', 'aj43f_dori2a_c1.mat', 'aj43f_dori3a_c1.mat'}, ...
     {'aj43g_dori1a_c1.mat', 'aj43g_dori2a_c1.mat', 'aj43g_dori3a_c1.mat'}}};
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
    figure('Name', 'Figure', 'Units', 'centimeters', 'Position', [0 0 12 11.5]);
end
unreswarning = warning('query', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle');


%% measured and fitted info scaling example
fprintf('Loading %s ...\n', [momentpath filesep momentfile]);
d = load([momentpath filesep momentfile]);
Iincr_mu = mean(d.Iincr_samples, 1);
Iincr_sd = sqrt(var(d.Iincr_samples, [], 1));
I_mu = [0 cumsum(Iincr_mu)];
I_sd = [0 sqrt(cumsum(Iincr_sd.^2))];
fprintf('Loading %s ...\n\n', [fitpath filesep fitfile]);
m = load([fitpath filesep fitfile]);
Nmax = 400;
prctsamples = 10000;

% model predictions
sms = m.sms{2};
assert(strcmp(sms.name, 'limlin'));  % make sure to get the right model
assert(strcmp(sms.pnames{1}, 'Iinf'));
Iincrs = NaN(prctsamples, Nmax);
Is = NaN(prctsamples, Nmax+1);
j = round(linspace(1, size(sms.mc.ss,1), prctsamples));
for i = 1:prctsamples
    p = sms.mc.ss(j(i),:,mod(i,size(sms.mc.ss,3))+1);
    Iincrs(i,:) = sms.Iincrfn(p, 1:Nmax);
    Is(i,:) = sms.Infn(p, 0:Nmax);
end
Iincrmode = sms.Iincrfn(sms.ml.p, 1:Nmax);
Imode = sms.Infn(sms.ml.p, 0:Nmax);

% fit statistics
assert(strcmp(m.sms{1}.name, 'lin'));
fprintf('linear model WAIC = %f\n', m.sms{1}.mc.WAIC1);
fprintf('limlin model WAIC = %f\n\n', sms.mc.WAIC1);

% Iincr data & fits
if singlefig
    subplotcm([1 10 4.5 1]);  hold on;
    text(-1,1.6,'a','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;
    subplot('Position', [0.13 0.7 0.775 0.225]);  hold on;
end
patch([0 sms.N sms.N 0], [-0.2 -0.2 0.45 0.45], 1, ...
    'FaceColor', obspopcol, 'EdgeColor', 'none');
alpha(errorshade(1:length(Iincr_mu), Iincr_mu, Iincr_sd, [1 1 1]*0.2), 0.2);
plot(1:length(Iincr_mu), Iincr_mu, '-', 'Color', [1 1 1]*0.2, 'LineWidth', 1);
for i = 1:length(plotprctiles)
    c = prctilesat(i)*doricols(doriexample,:) + (1-prctilesat(i))*[1 1 1];
    plot(1:Nmax, prctile(Iincrs, plotprctiles(i)), '-', 'Color', c);
end
plot(1:Nmax, Iincrmode, '--', 'Color', doricols(doriexample,:));
plot([1 Nmax], [0 0], '-', 'Color', [1 1 1]*0.5);
xlim([0 Nmax]);  ylim([-0.2 0.45]);
set(gca,'XColor','none','XTick',[],'Box','off','YTick',[-0.2 0 0.2 0.4]);
ylabel('Fisher information increase');

% I data & fits
if singlefig
    subplotcm([1 8 4.5 1.8]);  hold on;
else
    subplot('Position', [0.13 0.11 0.775 0.55]);  hold on;
end
patch([0 sms.N sms.N 0], [0 0 25 25], 1, ...
    'FaceColor', obspopcol, 'EdgeColor', 'none');
alpha(errorshade(0:(length(I_mu)-1), I_mu, I_sd, [1 1 1]*0.2), 0.2);
plot(0:(length(I_mu)-1), I_mu, '-', 'Color', [1 1 1]*0.2, 'LineWidth', 1);
for i = 1:length(plotprctiles)
    c = prctilesat(i)*doricols(doriexample,:) + (1-prctilesat(i))*[1 1 1];
    plot(0:Nmax, prctile(Is, plotprctiles(i)), '-', 'Color', c);
    text(Nmax*1.01, prctile(Is(:,end)', plotprctiles(i)), ...
        sprintf('%d%%', plotprctiles(i)), ...
        'HorizontalAlignment', 'left', 'Color', c);
end
plot(0:Nmax, Imode, '--', 'Color', doricols(doriexample,:));
text(Nmax*1.01, Imode(end), 'mode', ...
    'Color', doricols(doriexample,:), 'HorizontalAlignment', 'left');
xlim([0 Nmax]);  ylim([0 25]);
set(gca,'Box','off','XTick',[0 100 200 300 400],'YTick',[0 10 20],'Clipping','off');
xlabel('number of neurons N');
ylabel('Fisher information [rad^{-2}]');


%% large-N extraploation of example, with Iinf density plot
Nmax = 6000;
Nbar = 6500;
Ndensmin = 6800;
Ndensmax = 7800;
Iinfmax = 180;
% re-compute model predictions
Ns = linspace(0,Nmax,500);
Is = NaN(prctsamples, length(Ns));
Iinfs = NaN(1, prctsamples);
for i = 1:prctsamples
    p = sms.mc.ss(j(i),:,mod(i,size(sms.mc.ss,3))+1);
    Is(i,:) = sms.Infn(p, Ns);
    Iinfs(i) = p(1);
end

if singlefig
    subplotcm([1 4 4.5 3]);  hold on;
    text(-1,3.5,'c','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
% extraplocation to Nmax
patch([0 sms.N sms.N 0], [0 0 Iinfmax Iinfmax], 1, ...
    'FaceColor', obspopcol, 'EdgeColor', 'none');  hold on;
alpha(errorshade(0:(length(I_mu)-1), I_mu, I_sd, [1 1 1]*0.2), 0.2);
plot(0:(length(I_mu)-1), I_mu, '-', 'Color', [1 1 1]*0.2, 'LineWidth', 1);
Ninfprctiles = NaN(1, length(plotprctiles));
for i = 1:length(plotprctiles)
    c = prctilesat(i)*doricols(doriexample,:) + (1-prctilesat(i))*[1 1 1];
    plot(Ns, prctile(Is, plotprctiles(i)), '-', 'Color', c);
    Ninfprctiles(i) = prctile(Iinfs, plotprctiles(i));
    plot([Nmax Nbar], [prctile(Is(:,end)', plotprctiles(i)) Ninfprctiles(i)], ...
        '--', 'Color', c);
end
% bar at Nbar
c1 = prctilesat(1)*doricols(doriexample,:) + (1-prctilesat(1))*[1 1 1];
plot([Nbar Nbar], [Ninfprctiles(1) Ninfprctiles(5)], '-', ...
    'LineWidth', 0.5, 'Color', c1);
c2 = prctilesat(2)*doricols(doriexample,:) + (1-prctilesat(2))*[1 1 1];
plot([Nbar Nbar], [Ninfprctiles(2) Ninfprctiles(4)], '-', ...
    'LineWidth', 1, 'Color', c2);
c3 = prctilesat(3)*doricols(doriexample,:) + (1-prctilesat(3))*[1 1 1];
plot(Nbar, Ninfprctiles(3), 'o', 'MarkerSize', 4, ...
    'MarkerEdgeColor', 'none', 'MarkerFaceColor', c3);
% density from Ndensmin to Ndensmax
[pIinf,densIinf] = ksdensity(reshape(sms.mc.ss(:,1,:),1,[]), ...
    linspace(0,Iinfmax,1000),'Support','Positive');
pIinf = Ndensmin + (Ndensmax - Ndensmin) * pIinf ./ max(pIinf);
di = NaN(1, length(plotprctiles));
for i = 1:length(plotprctiles)
    [~,di(i)] = min(abs(densIinf - Ninfprctiles(i)));
end
alpha(patch([(ones(1,di(2)-di(1)+1)*Ndensmin) fliplr(pIinf(di(1):di(2)))], ...
    [densIinf(di(1):di(2)) fliplr(densIinf(di(1):di(2)))],1, ...
    'FaceColor', c1, 'EdgeColor', 'none'), 0.5);
alpha(patch([(ones(1,di(4)-di(2)+1)*Ndensmin) fliplr(pIinf(di(2):di(4)))], ...
    [densIinf(di(2):di(4)) fliplr(densIinf(di(2):di(4)))],1, ...
    'FaceColor', c2, 'EdgeColor', 'none'), 0.5);
alpha(patch([(ones(1,di(5)-di(4)+1)*Ndensmin) fliplr(pIinf(di(4):di(5)))], ...
    [densIinf(di(4):di(5)) fliplr(densIinf(di(4):di(5)))],1, ...
    'FaceColor', c1, 'EdgeColor', 'none'), 0.5);
plot([Ndensmin pIinf(di(3))], [1 1]*densIinf(di(3)), '-', ...
    'LineWidth', 1, 'Color', c3);
plot(pIinf, densIinf, '-', 'Color', [0.5 0.5 0.5]);
plot([Ndensmin Ndensmin], [0 Iinfmax], '-', 'Color', [0.5 0.5 0.5]);
for i = 1:length(plotprctiles)
    c = prctilesat(i)*doricols(doriexample,:) + (1-prctilesat(i))*[1 1 1];
    text(pIinf(di(i)) + 0.1 * (Ndensmax - Ndensmin), densIinf(di(i)), ...
        sprintf('%d%%', plotprctiles(i)), 'HorizontalAlignment', 'left', 'Color', c);
end
% additional settings
xlim([0 Ndensmax]); ylim([0 Iinfmax]);
set(gca,'Box','off','XTick',[0 1000 2000 3000 4000 5000 6000 Nbar],...
    'XTickLabel',{'0','1000','2000','3000','4000','5000','6000','\infty'},...
    'YTick',[0 50 100 150]);
xlabel('Number of neurons N');
ylabel('Fisher information [rad^{-2}]');


%% 1/In over 1/N
invNs = 1./(1:(length(I_mu)-1));
% data moments, using approx. <1/x>=1/<x>, var(1/x)=var(x)/<x>^4
invI_mu = 1 ./ I_mu(2:end);
invI_sd = I_sd(2:end) ./ I_mu(2:end).^2;
% best linear fit, with inverse-variance weights
mdl = fitlm(invNs', invI_mu', 'Weights', (1./invI_sd.^2)');
fprintf('Model fit for 1/In over 1/N linear fit (inv-variance-weighted):\n');
disp(mdl);
fprintf('p-value of intercept: %f\n\n', mdl.Coefficients.pValue(1));
b = mdl.Coefficients.Estimate;
% only plot from Nmin to end
Nmin = 100;
invNs = invNs(Nmin:end);
invI_mu = invI_mu(Nmin:end);
invI_sd = invI_sd(Nmin:end);
% generate plot
if singlefig
    subplotcm([7 8 4 3]);  hold on;
    text(-1,3.5,'b','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
patch(1./[Nmin sms.N sms.N Nmin], [0 0 0.15 0.15], 1, ...
    'FaceColor', obspopcol, 'EdgeColor', 'none');
alpha(errorshade(invNs, invI_mu, invI_sd, [1 1 1]*0.2), 0.2); hold on;
plot(invNs, invI_mu, '-', 'Color', [1 1 1]*0.2, 'LineWidth', 1);
plot([0 1/Nmin], [b(1) (b(1)+b(2)/Nmin)], '--', 'Color', [0.5 0.5 0.5]);
plot(0, b(1), 'o', 'MarkerSize', 4, 'MarkerFaceColor', 0.5*[1 1 1], ...
    'MarkerEdgeColor', 'none');
text(0, b(1)+0.01*0.15, 'I_\infty', 'Color', 0.5*[1 1 1],...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
xlim([0 1./Nmin]); ylim([0 0.15]);
set(gca,'Box','off');
set(gca,'XTick',[0 0.001 0.003333 0.005 0.01],...
    'XTickLabel',{'\infty','1000','300','200','100'});
set(gca,'YTick',[0 0.05 0.10 0.15],'YTickLabel',{'\infty','20','10','6.67'});
xlabel('number of neurons N (1/N scale)');
ylabel('Fisher information (1/IN scale)');


%% different densities
Iinfmin = 10;  Iinfmax = 10^4;
Iinfs = logspace(log10(Iinfmin), log10(Iinfmax), 500);
barymin = 0.8;  barydist = 0.085;
ymax = barymin+(length(dorifiles)-0.5)*barydist;
if singlefig
    subplotcm([6.5 4 4.5 3]);  hold on;
    text(-0.5,3.5,'d','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
% orientation discrimination thresholds
for i = 1:length(orithreshs)
    x = 2*norminv(0.8)^2/(orithreshs(i)*pi/180)^2;
    plot([1 1]*x, [0 ymax], '-', 'LineWidth', 0.5, 'Color', [1 1 1]*0.5);
    text(x, ymax, sprintf('%d\\circ', orithreshs(i)), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
        'Color', [1 1 1]*0.5);
end
% for each data file, plot density and bar
for i = 1:length(dorifiles)
    fprintf('Loading %s ...\n', [fitpath filesep dorifiles{i}]);
    m = load([fitpath filesep dorifiles{i}]);
    if contains(dorifiles{i}, 'dori'), sms = m.fits{2};
    else, sms = m.sms{2}; end
    assert(strcmp(sms.name, 'limlin'));  % make sure to get the right model
    assert(strcmp(sms.pnames{1}, 'Iinf'));
    ss = reshape(sms.mc.ss(:,1,:),1,[]);
    pIinf = ksdensity(ss,Iinfs,'Support','Positive');
    % re-scale density to justify log-scaling p(logI)=p(I)*exp(logI)
    plot(Iinfs, pIinf .* Iinfs, '-', 'Color', doricols(i,:));
    % summary bar plot
    bary = barymin + barydist * (i-1);
    plot(prctile(ss, plotprctiles([1 5])), [1 1]*bary, '-', ...
        'Color', prctilesat(1)*doricols(i,:) + (1-prctilesat(1))*[1 1 1]);
    plot(prctile(ss, plotprctiles([2 4])), [1 1]*bary, '-', 'LineWidth', 1.4, ...
        'Color', prctilesat(2)*doricols(i,:) + (1-prctilesat(2))*[1 1 1]);
    plot(prctile(ss, plotprctiles(3)), bary, 'o', 'MarkerSize', 3, ...
        'MarkerFaceColor', prctilesat(3)*doricols(i,:) + (1-prctilesat(3))*[1 1 1], ...
        'MarkerEdgeColor', 'None');
    text(Iinfmax, bary, dorinames{i}, 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', 'Color', doricols(i,:));
end
fprintf('\n');
xlim([Iinfmin Iinfmax]);  ylim([0 ymax]);
set(gca,'Box','off','XScale','log','YColor','none','YTick',[]);
xlabel('Fisher information');


%% per-animal across-session plots
pw = (10 - 0.5*(animals-1))/animals;
ph = 2.2;
for i = 1:animals
    sessions = length(sessfiles{i});
    sx = @(j) 0.04*(j-1)-0.02*(sessions-1);
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
            fprintf('Loading %s ...\n', [fitpath filesep sessfiles{i}{j}{k}]);
            m = load([fitpath filesep sessfiles{i}{j}{k}]);
            sms = m.fits{2};
            assert(strcmp(sms.name, 'limlin'));  % make sure to get the right model
            assert(strcmp(sms.pnames{1}, 'Iinf'));
            ss = reshape(sms.mc.ss(:,1,:),1,[]);
            meds(k) = prctile(ss, plotprctiles(3));
            plotPrctiles(k+sx(j), ss, plotcols(i,j));
        end
        plot((1:length(sessdorinames))+sx(j), meds, '-', 'Color', plotcols(i,j));
    end
    % plot formatting
    set(gca,'Box','off','YScale','log',...
        'XTick',1:length(sessdorinames),'XTickLabel',sessdorinames,...
        'YTick',[1:9 10:10:90 100:100:1000],'YTickLabel',...
        {'10^0','','','','','','','','',...
         '10^1','','','','','','','','',...
         '10^2','','','','','','','','','10^3'});
    if i == 1
        xlabel('drift direction difference');
        ylabel('Fisher information [rad^{-2}]');
    else
        set(gca, 'YTickLabel',[]);
    end
    xlim([0.75 length(sessdorinames)+0.25]);  ylim([1 1000]);
    text(length(sessdorinames), 1, sprintf('mouse %d', i), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end


%% save figure
warning(unreswarning.state, 'MATLAB:dispatcher:UnresolvedFunctionHandle');
if singlefig
    fprintf('\nWriting figure to fig4.pdf\n');
    print(['figs' filesep 'fig4'], '-dpdf');
end


function o = errorshade(x, ymu, ysd, c)
% plots shaded error regions ymu+/-ysd in color c over x
o = patch([x(:)' fliplr(x(:)')], [(ymu(:)'+ysd(:)') fliplr(ymu(:)'-ysd(:)')],1,...
    'FaceColor',c,'EdgeColor','none');
