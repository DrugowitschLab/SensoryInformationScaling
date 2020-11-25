function plotFigS7
%% plots the model fits to simulations

%% general settings
singlefig = true;  % format for paper if true
sims = {'sim1','sim2','sim3','sim4'};
Ns = [100 300 500 750 1000];
Ts = {[100 300 500 750 1000], ...
      [300 500 750 1000], ...
      [500 750 1000], ...
      [750 1000], ...
      1000};
Tnum = sum(arrayfun(@(i) length(Ts{i}), 1:length(Ts)));
Tcols1 = [208 209 230; 166 189 219; 103 169 207; 28 144 153; 1 108 89] / 255;
Tcols2 = [253 212 158; 253 187 132; 252 141 8; 227 74 51; 179 0 0] / 255;
Tplotsep = 0.1;
fitpath = ['.' filesep 'fits'];
LNPdata = ['.' filesep 'simData' filesep 'sim3.mat'];
Gaussdata = ['.' filesep 'simData', filesep 'sim1.mat'];
examplemom = 'N300_T500';
momfolder = ['.' filesep 'moment_cache'];
figpos = [1 8.7 1 6.2; 1 3.2 1 0.8; 5 8.7 5 6.2; 5 3.2 5 0.8];


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


%% get Iinf for the different datasets
fprintf('Finding Iinf for datasets\n');
fprintf('Loading %s ... \n', Gaussdata);
d = load(Gaussdata);
GaussIinf = d.par.Iinf;
fprintf('Loading %s ... \n', LNPdata);
d = load(LNPdata);
N = d.par.N;
T = d.par.T;
LNPIinf = (2*T-N-3)/(2*T-2) * (d.popmom.fp * (d.popmom.Sig \ d.popmom.fp')) - ...
    (2*N)/(T*d.par.ds^2);
simIinfs = [GaussIinf NaN LNPIinf NaN];


%% plot example scaling curves
% Gauss
if singlefig
    figure('Name', 'Figure', 'Units', 'centimeters', 'Position', [0 0 8.5 14.2]);    
    subplotcm([1 11.7 3 2]);  hold on;
    text(-1,2.5,'a','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
[Imul, Isdl] = Imoms([momfolder filesep 'sim1_' examplemom '.mat']);
[Imuu, Isdu] = Imoms([momfolder filesep 'sim2_' examplemom '.mat']);
N = length(Imul);
plot([1 N], [1 1]*GaussIinf, '--', 'Color', [1 1 1]*0.5);
alpha(errorshade(1:N, Imuu, Isdu, Tcols2(end,:)), 0.2);
alpha(errorshade(1:N, Imul, Isdl, Tcols1(end,:)), 0.2);
plot(1:N, Imuu, '-', 'Color', Tcols2(end,:), 'LineWidth', 1);
plot(1:N, Imul, '-', 'Color', Tcols1(end,:), 'LineWidth', 1);
set(gca,'Box','off','XTick',[1 100 200 300]);
xlim([1 N]);  ylim([0 45]);
xlabel('number of neurons N');
ylabel('Fisher information [rad-2]');

% LNP
if singlefig
    subplotcm([5 11.7 3 2]);  hold on;
    text(-1,2.5,'d','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
[Imul, Isdl] = Imoms([momfolder filesep 'sim3_' examplemom '.mat']);
[Imuu, Isdu] = Imoms([momfolder filesep 'sim4_' examplemom '.mat']);
N = length(Imul);
plot([1 N], [1 1]*LNPIinf, '--', 'Color', [1 1 1]*0.5);
alpha(errorshade(1:N, Imuu, Isdu, Tcols2(end,:)), 0.2);
alpha(errorshade(1:N, Imul, Isdl, Tcols1(end,:)), 0.2);
plot(1:N, Imuu, '-', 'Color', Tcols2(end,:), 'LineWidth', 1);
plot(1:N, Imul, '-', 'Color', Tcols1(end,:), 'LineWidth', 1);
set(gca,'Box','off','XTick',[1 100 200 300]);
xlim([1 N]);  ylim([0 45]);
xlabel('number of neurons N');
ylabel('Fisher information [rad-2]');


%% show fits for different models/Ts
Nstr = arrayfun(@(i) sprintf('%d', Ns(i)), 1:length(Ns), 'UniformOutput', false);
for i = 1:length(sims)
    fprintf('--- Processing %s ...\n', sims{i});
    % Fisher information
    if singlefig
        subplotcm([figpos(i,1) figpos(i,2) 3 2]);  hold on;
        text(-1,2.5,'b','Units','centimeters','FontWeight','bold',...
            'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
    else
        figure;  hold on;
    end
    plot([0.5 length(Ns)+0.5], [1 1]*simIinfs(i), '--', 'Color', [1 1 1]*0.5);
    cstore = cell(1,Tnum);
    linwin = false(1,Tnum);
    ci = 1;
    for j = 1:length(Ns)
        Tsj = Ts{j};
        Txs = j + Tplotsep * ((1:length(Tsj)) - 0.5 - length(Tsj)/2);
        for k = 1:length(Tsj)
            [Iinfs, cs, linwin(ci)] = loadFitParams(sims{i}, fitpath, Ns(j), Tsj(k));
            kc = size(Tcols1,1)-length(Tsj)+k;
            if linwin(ci), c = Tcols1(kc,:); else, c = Tcols2(kc,:); end
            plotPrctiles(Txs(k), Iinfs, c);
            cstore{ci} = cs;
            ci = ci + 1;
        end
    end
    set(gca,'Box','off','XColor','none','YScale','log');
    xlim([0.5 length(Ns)+0.5]);
    ylabel('Fisher information [rad-2]');
    
    % cs
    if singlefig
        subplotcm([figpos(i,3) figpos(i,4) 3 2]);  hold on;
    else
        figure;  hold on;
    end
    ci = 1;
    for j = 1:length(Ns)
        Tsj = Ts{j};
        Txs = j + Tplotsep * ((1:length(Tsj)) - 0.5 - length(Tsj)/2);
        for k = 1:length(Tsj)
            kc = size(Tcols1,1)-length(Tsj)+k;
            if linwin(ci), c = Tcols1(kc,:); else, c = Tcols2(kc,:); end
            plotPrctiles(Txs(k), cstore{ci}, c);
            ci = ci + 1;
        end
    end
    set(gca,'Box','off','XTick',1:length(Ns),'XTickLabel',Nstr);
    xlabel('population size N');    xlim([0.5 length(Ns)+0.5]);
    ylabel('estimated c');
end


%% write figure to file
if singlefig
    fprintf('\nWriting figure to figS7.pdf\n');
    print(['figs' filesep 'figS7'], '-dpdf');
end


function [Iinfs, cs, linwin] = loadFitParams(simname, fitpath, N, T)
%% loads fit data for limlin fit and returns posterior parameter samples

fitfile = sprintf('%s%s%s_N%d_T%d_norm.mat', fitpath, filesep, simname, N, T);
fprintf('Loading %s ...\n', fitfile);
unreswarning = warning('query', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
d = load(fitfile);
warning(unreswarning.state, 'MATLAB:dispatcher:UnresolvedFunctionHandle');
smslin = d.sms{1};
sms = d.sms{2};
assert(strcmp(smslin.name, 'lin'));
assert(strcmp(sms.name, 'limlin'));  % ensure picking the richt model
assert(strcmp(sms.pnames{1}, 'Iinf'));
assert(strcmp(sms.pnames{2}, 'c'));
Iinfs = reshape(sms.mc.ss(:,1,:),1,[]);
cs = reshape(sms.mc.ss(:,2,:),1,[]);
linwin = smslin.mc.WAIC1 > sms.mc.WAIC1;


function [Imu, Isd] = Imoms(momfile)
%% computes mean I and its SD for the given file

fprintf('Loading %s ...\n', momfile);
d = load(momfile);
Iincr_mu = mean(d.Iincr_samples,1);
Iincr_var = var(d.Iincr_samples,[],1);
Imu = cumsum(Iincr_mu);
Isd = sqrt(cumsum(Iincr_var));


function o = errorshade(x, ymu, ysd, c)
% plots shaded error regions ymu+/-ysd in color c over x
o = patch([x(:)' fliplr(x(:)')], [(ymu(:)'+ysd(:)') fliplr(ymu(:)'-ysd(:)')],1,...
    'FaceColor',c,'EdgeColor','none');
