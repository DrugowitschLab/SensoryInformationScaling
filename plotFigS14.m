function plotFigS14

%% general settings
singlefig = true;   % set to true to have all planels in single figure

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
    figure('Name', 'Figure', 'Units', 'centimeters', 'Position', [0 0 16 8.8]);
end

%% model parameters
% network parameters
P = 32;
par = struct(...
    'N', 1000, ...    % number of neurons
    'P', P, ...       % pixels in image (one side)
    'sig', P/5, ...   % gaussian envelope sd
    'lam', P/1.5, ... % pref. wavelength
    'phi', 0, ...     % pref. spatial phase
    'c', 1, ...       % Michelson contrast
    'g', 10, ...      % tuning amplitude
    'sig_a', sqrt(2), ... % tuning ampl. variability
    'sig0', 0.11, ... % input noise (prop. of pixel range)
    'ds', 15 * pi / 180, ...
    'T', 20000);
% simulation parameters
N = par.N;
% 
% plot colors and saturation function
plotIlims = [0 150];
oricol = [1 1 1]*0.2;
oricol_em = [0 1 0];

%%
pn_list = 0.1: 0.1: 0.8; 
col_pt = [27 158 119; 217 95 2; 117 112 179] / 255;

%% Info scaling for pt=0.3
if singlefig
    subplotcm([1 4.8 4 3]);  hold on;
else
    figure('Color', 'white');  hold on;
end

pt = 0; 
pn = 0; 
outfile = ['.' filesep 'moment_cache' filesep 'sim_eye_jitter_pn_' num2str(pn) '_pt_' num2str(pt) '_N_' num2str(N) '.mat'];

fprintf('Information scaling for spike count data ...\n');
d = load(outfile);

Imuori = [0 mean(cumsum(d.Iincr_samples, 2))];
plot(0:par.N, Imuori, '-', 'Color', oricol, 'LineWidth', 1);
hold on; 


%%
pt = 0.3;
for jj = 1: length(pn_list)
    pn = pn_list(jj);
    outfile = ['.' filesep 'moment_cache' filesep 'sim_eye_jitter_pn_' num2str(pn) '_pt_' num2str(pt) '_N_' num2str(N) '.mat'];

    fprintf('Information scaling for spike count data ...\n');
    d = load(outfile);

    Imu = [0 mean(cumsum(d.Iincr_samples, 2))];
    plot(0:par.N, Imu, '-', 'Color', [col_pt(1, :), (1- pn)], 'LineWidth', 0.5);
end

xlim([0 N]);  ylim(plotIlims);
title('p_t = 0.3');
set(gca,'Box','off','XTick',[0 N/2 N],'YTick',[0 50 100 150]);


%% Info scaling for pt=0.6
if singlefig
    subplotcm([6 4.8 4 3]);  hold on;
else
    figure('Color', 'white');  hold on;
end

pt = 0; 
pn = 0; 
outfile = ['.' filesep 'moment_cache' filesep 'sim_eye_jitter_pn_' num2str(pn) '_pt_' num2str(pt) '_N_' num2str(N) '.mat'];

fprintf('Information scaling for spike count data ...\n');
d = load(outfile);

Imuori = [0 mean(cumsum(d.Iincr_samples, 2))];
plot(0:par.N, Imuori, '-', 'Color', oricol, 'LineWidth', 1);
hold on; 

pt = 0.6;
for jj = 1: length(pn_list)
    pn = pn_list(jj);
    outfile = ['.' filesep 'moment_cache' filesep 'sim_eye_jitter_pn_' num2str(pn) '_pt_' num2str(pt) '_N_' num2str(N) '.mat'];

    fprintf('Information scaling for spike count data ...\n');
    d = load(outfile);

    Imu = [0 mean(cumsum(d.Iincr_samples, 2))];
    plot(0:par.N, Imu, '-', 'Color', [col_pt(2, :), (1- pn)], 'LineWidth', 0.5);
end

xlim([0 N]);  ylim(plotIlims);
title('p_t = 0.6');
set(gca,'Box','off','XTick',[0 N/2 N],'YColor','none');


%% Info scaling for pt=1
if singlefig
    subplotcm([11 4.8 4 3]);  hold on;
else
    figure('Color', 'white');  hold on;
end

pt = 0; 
pn = 0; 
outfile = ['.' filesep 'moment_cache' filesep 'sim_eye_jitter_pn_' num2str(pn) '_pt_' num2str(pt) '_N_' num2str(N) '.mat'];

fprintf('Information scaling for spike count data ...\n');
d = load(outfile);

Imuori = [0 mean(cumsum(d.Iincr_samples, 2))];
plot(0:par.N, Imuori, '-', 'Color', oricol, 'LineWidth', 1);
hold on; 

pt = 1;
for jj = 1: length(pn_list)
    pn = pn_list(jj);
    outfile = ['.' filesep 'moment_cache' filesep 'sim_eye_jitter_pn_' num2str(pn) '_pt_' num2str(pt) '_N_' num2str(N) '.mat'];

    fprintf('Information scaling for spike count data ...\n');
    d = load(outfile);

    Imu = [0 mean(cumsum(d.Iincr_samples, 2))];
    plot(0:par.N, Imu, '-', 'Color', [col_pt(3, :), (1- pn)], 'LineWidth', 0.5);
end

xlim([0 N]);  ylim(plotIlims);
title('p_t = 1');
set(gca,'Box','off','XTick',[0 N/2 N],'YColor','none');


%% fit plots for c
if singlefig
    subplotcm([1 0.8 4 3]);  hold on;
else
    figure('Color', 'white');  hold on;
end

fitpath = ['.' filesep 'fits'];
afrac = 0.95;
obspopcol = [0.95 0.95 0.99];
plotprctiles = [5 25 50 75 95];
cprcriles = [50 25 75];

pt = 0; 
pn = 0; 
fitname = ['sim_eye_jitter_pn_' num2str(pn) '_pt_' num2str(pt)  '_N_' num2str(N) '_norm'];
fprintf('Loading %s ...\n', [fitpath filesep fitname]);
m = load([fitpath filesep fitname]);
sms = m.sms{2};
assert(strcmp(sms.name, 'limlin'));  % make sure to get the right model
assert(strcmp(sms.pnames{1}, 'Iinf'));
assert(strcmp(sms.pnames{2}, 'c'));
Iinfss = reshape(sms.mc.ss(:,1,:),1,[]);
css = reshape(sms.mc.ss(:,2,:),1,[]);
N95ss = (afrac/(1-afrac)) * Iinfss ./ css;
meds = prctile(N95ss, plotprctiles(3));

hold on; plotPrctiles(0, css', [0,0,0]);
yline(median(css),'k');

pt_scaling = [0.3, 0.6, 1];
css_median = zeros(length(pn_list), 3);
Iinf_median = zeros(length(pn_list), 3);
N95ss_median = zeros(length(pn_list), 3);

for ii = 1: 3
    pt = pt_scaling(ii);
    for jj = 1: length(pn_list)
        pn = pn_list(jj);
        fitname = ['sim_eye_jitter_pn_' num2str(pn) '_pt_' num2str(pt) '_N_' num2str(N) '_norm.mat'];
        fprintf('Loading %s ...\n', [fitpath filesep fitname]);
        m = load([fitpath filesep fitname]);
        sms = m.sms{2};
        assert(strcmp(sms.name, 'limlin'));  % make sure to get the right model
        assert(strcmp(sms.pnames{1}, 'Iinf'));
        assert(strcmp(sms.pnames{2}, 'c'));
        Iinfss = reshape(sms.mc.ss(:,1,:),1,[]);
        css = reshape(sms.mc.ss(:,2,:),1,[]);
        N95ss = (afrac/(1-afrac)) * Iinfss ./ css;
        hold on; plotPrctiles(pn+ 0.01*(ii-2) , css', col_pt(ii, :));
        css_median(jj, ii) = median(css);
        Iinf_median(jj, ii) = median(Iinfss);
        N95ss_median(jj, ii) = median(N95ss);
    end
end
xlim([0 1]);  ylim([0 2]);
% title('estimated c');
set(gca,'Box','off','XTick',[0 0.5 1],'YTick',[0 1 2]);
xlabel('p_n');
ylabel('c');
b1 = plot(pn_list- 0.01*1, css_median(:, 1), 'LineWidth', 0.5, 'Color', col_pt(1, :));
b2 = plot(pn_list, css_median(:, 2), 'LineWidth', 0.5, 'Color', col_pt(2, :));
b3 = plot(pn_list+ 0.01*1, css_median(:, 3), 'LineWidth', 0.5, 'Color', col_pt(3, :));



%% fit plots for Iinf
if singlefig
    subplotcm([6 0.8 4 3]);  hold on;
else
    figure('Color', 'white');  hold on;
end


pt = 0; 
pn = 0; 
fitname = ['sim_eye_jitter_pn_' num2str(pn) '_pt_' num2str(pt)  '_N_' num2str(N) '_norm'];
fprintf('Loading %s ...\n', [fitpath filesep fitname]);
m = load([fitpath filesep fitname]);
sms = m.sms{2};
assert(strcmp(sms.name, 'limlin'));  % make sure to get the right model
assert(strcmp(sms.pnames{1}, 'Iinf'));
assert(strcmp(sms.pnames{2}, 'c'));
Iinfss = reshape(sms.mc.ss(:,1,:),1,[]);
css = reshape(sms.mc.ss(:,2,:),1,[]);
N95ss = (afrac/(1-afrac)) * Iinfss ./ css;
meds = prctile(N95ss, plotprctiles(3));

hold on; plotPrctiles(0, Iinfss', [0,0,0]);
yline(median(Iinfss),'k');

for ii = 1: 3
    pt = pt_scaling(ii);
    for jj = 1: length(pn_list)
        pn = pn_list(jj);
        fitname = ['sim_eye_jitter_pn_' num2str(pn) '_pt_' num2str(pt) '_N_' num2str(N) '_norm.mat'];
        fprintf('Loading %s ...\n', [fitpath filesep fitname]);
        m = load([fitpath filesep fitname]);
        sms = m.sms{2};
        assert(strcmp(sms.name, 'limlin'));  % make sure to get the right model
        assert(strcmp(sms.pnames{1}, 'Iinf'));
        assert(strcmp(sms.pnames{2}, 'c'));
        Iinfss = reshape(sms.mc.ss(:,1,:),1,[]);
        css = reshape(sms.mc.ss(:,2,:),1,[]);
        N95ss = (afrac/(1-afrac)) * Iinfss ./ css;
        hold on; plotPrctiles(pn+ 0.01*(ii-2) , Iinfss', col_pt(ii, :));
        css_median(jj, ii) = median(css);
        Iinf_median(jj, ii) = median(Iinfss);
        N95ss_median(jj, ii) = median(N95ss);
    end
end
xlim([0 1]);  ylim(plotIlims);
% title('estimated I_\infty');
xlabel('p_n');
ylabel('I_\infty');
b1 = plot(pn_list- 0.01*1, Iinf_median(:, 1), 'LineWidth', 0.5, 'Color', col_pt(1, :));
b2 = plot(pn_list, Iinf_median(:, 2), 'LineWidth', 0.5, 'Color', col_pt(2, :));
b3 = plot(pn_list+ 0.01*1, Iinf_median(:, 3), 'LineWidth', 0.5, 'Color', col_pt(3, :));
set(gca,'Box','off','XTick',[0 0.5 1],'YTick',[0 50 100 150]);
% legend([b1 b2 b3],'p_t = 0.3','p_t = 0.6', 'p_t = 1', 'Location', 'best');

%% fit plots for N95
if singlefig
    subplotcm([11 0.8 4 3]);  hold on;
else
    figure('Color', 'white');  hold on;
end


pt = 0; 
pn = 0; 
fitname = ['sim_eye_jitter_pn_' num2str(pn) '_pt_' num2str(pt)  '_N_' num2str(N) '_norm'];
fprintf('Loading %s ...\n', [fitpath filesep fitname]);
m = load([fitpath filesep fitname]);
sms = m.sms{2};
assert(strcmp(sms.name, 'limlin'));  % make sure to get the right model
assert(strcmp(sms.pnames{1}, 'Iinf'));
assert(strcmp(sms.pnames{2}, 'c'));
Iinfss = reshape(sms.mc.ss(:,1,:),1,[]);
css = reshape(sms.mc.ss(:,2,:),1,[]);
N95ss = (afrac/(1-afrac)) * Iinfss ./ css;
meds = prctile(N95ss, plotprctiles(3));

hold on; plotPrctiles(0, N95ss', [0,0,0]);
yline(median(N95ss),'k')
for ii = 1: 3
    pt = pt_scaling(ii);
    for jj = 1: length(pn_list)
        pn = pn_list(jj);
        fitname = ['sim_eye_jitter_pn_' num2str(pn) '_pt_' num2str(pt) '_N_' num2str(N) '_norm.mat'];
        fprintf('Loading %s ...\n', [fitpath filesep fitname]);
        m = load([fitpath filesep fitname]);
        sms = m.sms{2};
        assert(strcmp(sms.name, 'limlin'));  % make sure to get the right model
        assert(strcmp(sms.pnames{1}, 'Iinf'));
        assert(strcmp(sms.pnames{2}, 'c'));
        Iinfss = reshape(sms.mc.ss(:,1,:),1,[]);
        css = reshape(sms.mc.ss(:,2,:),1,[]);
        N95ss = (afrac/(1-afrac)) * Iinfss ./ css;
        hold on; plotPrctiles(pn+ 0.01*(ii-2) , N95ss', col_pt(ii, :));
        css_median(jj, ii) = median(css);
        Iinf_median(jj, ii) = median(Iinfss);
        N95ss_median(jj, ii) = median(N95ss);
    end
end
xlim([0 1]);  ylim([1e3 3e4]);
% title('estimated N_{95}');
xlabel('p_n');
ylabel('N_{95}');
b1 = plot(pn_list- 0.01*1, N95ss_median(:, 1), 'LineWidth', 0.5, 'Color', col_pt(1, :));
b2 = plot(pn_list, N95ss_median(:, 2), 'LineWidth', 0.5, 'Color', col_pt(2, :));
b3 = plot(pn_list+ 0.01*1, N95ss_median(:, 3), 'LineWidth', 0.5, 'Color', col_pt(3, :));
legend([b1 b2 b3],'p_t = 0.3','p_t = 0.6', 'p_t = 1', 'Location', 'best');
set(gca,'Box','off','XTick',[0 0.5 1],'YTick',[1000 5000 1e4 2e4 3e4],'YScale','log');


%% save figure
if singlefig
    fprintf('\nWriting figure to figS14.pdf\n');
    print(['figs' filesep 'figS14'], '-dpdf');
end