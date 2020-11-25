function plotFigS12
%% plots supplementary Fig. showing tuning curves and pairwise corr.

%% (general) settings
singlefig = true;   % set to true to have all planels in single figure
tcfitpath = ['.' filesep 'tuning_fits'];
dataset = 'aj60a';
datasets = {{'aj60a', 'aj60b', 'aj60c', 'aj60d'}, {'aj61a', 'aj61b', 'aj61c'}};

subsample = 'none';
% conid = 1;  % contrast to use
eventRateScale  = 30;
ylimmax = 0.2;   % limit for the y-axis for tuning plots
ymulimmax = 0.2; % limit for the y-axis for sequential effects plots
uncol = [27 158 119] / 255;
dircol = [217 95 2] / 255;
oricol = [117 112 179] / 255;
oricols = [247 244 249; 129 144 190;  12  44 132;  78  22  98; ...
           145   0  63; 161   0  50; 177   0  38; 212 122 150] / 255;
tuningrows = 5;
tuningcols = 2;
cdfsubsamples = 10000;

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
set(groot, 'defaultFigureRenderer', 'painters');
if singlefig
    figure('Name', 'Figure', 'Units', 'centimeters', 'Position', [0 0 21 11.5]);
end


%% loading data
% activity data
fprintf('Loading %s ...\n', dataset');
[sr, oris, cons] = loaddata(dataset,subsample);
uoris = unique(oris);
ucons = unique(cons);
consu = ucons;
ncons = length(ucons);


%% von-mises functions are the same for both contrast
tcfitfile_sep = [tcfitpath filesep dataset filesep 'TuningCoef_Subsample_' subsample '.mat'];
fprintf('Loading individual fits %s ...\n', tcfitfile_sep);
tcd_sep = load(tcfitfile_sep);

conid = 1;

vmfun.OriTun_vonMises   = tcd_sep.coef.OriTun_vonMises{conid};
vmfun.DirTun_vonMises   = tcd_sep.coef.DirTun_vonMises{conid};
vmfun.NullTun           = tcd_sep.coef.NullFunc{conid};

%% select neuron tuning based on high-contrast stim
conid = 2; 

nonsigneurons = find(tcd_sep.coef.alpha_OriVsNull{conid} >= 0.05);
orineurons = find(tcd_sep.coef.alpha_OriVsNull{conid} < 0.05 & ...
    tcd_sep.coef.alpha_OriVsDir{conid} > 0.05);
dirneurons = find(tcd_sep.coef.alpha_OriVsNull{conid} < 0.05 & ...
    tcd_sep.coef.alpha_OriVsDir{conid} < 0.05);

%% joint tuning curve fits
tcfitfile = [tcfitpath filesep dataset filesep 'TuningCoef_Combo_Subsample_' subsample '.mat'];
fprintf('Loading joint fits %s ...\n', tcfitfile);
tcd = load(tcfitfile);

orientationForTrials = oris;
averageSpkRateOverTrial = sr; 
nNeurons = size(averageSpkRateOverTrial, 2);


x_OriTun = nan(6, nNeurons);
for ii = 1: length(orineurons)
    x_OriTun(:, orineurons(ii))    = tcd.coef.x_combo{orineurons(ii), 1};
end

x_DirTun = nan(7, nNeurons);
for ii = 1: length(dirneurons)
    x_DirTun(:, dirneurons(ii))    = tcd.coef.x_combo{dirneurons(ii), 1};
end

for con = 1: ncons
    x_null(:, :, con)      = tcd_sep.coef.x_null{con,1};
end

%% compute R2 of fits
R2_OriTun = NaN(ncons, nNeurons);
R2_DirTun = NaN(ncons, nNeurons);
f = NaN(1, length(uoris));
f_ori = NaN(ncons, length(uoris), nNeurons);
f_dir = NaN(ncons, length(uoris), nNeurons);

for nn = 1: nNeurons
    
    if sum(dirneurons == nn)

        [mu_lo, mu_hi] = fitMultip(x_DirTun(:, nn), vmfun, uoris, uoris, 'DirTun');
        f_dir(1, :, nn) = mu_lo/ eventRateScale;
        f_dir(2, :, nn) = mu_hi/ eventRateScale;
        
    elseif sum(orineurons == nn)
        [mu_lo, mu_hi] = fitMultip(x_OriTun(:, nn), vmfun, uoris, uoris, 'OriTun');
        f_ori(1, :, nn) = mu_lo/ eventRateScale;
        f_ori(2, :, nn) = mu_hi/ eventRateScale;
        
    end
    for conid = 1: 2
        for orii = 1:length(uoris)
            f(orii) = mean(sr(oris == uoris(orii) & cons == ucons(conid), nn));
        end
        R2_OriTun(conid, nn) = 1 - sum((f_ori(conid, :, nn) - f).^2) / sum((f - mean(f)).^2);
        R2_DirTun(conid, nn) = 1 - sum((f_dir(conid, :, nn) - f).^2) / sum((f - mean(f)).^2);
    end
end


%% plot activity / tuning curves of untuned and tuned neurons for example data
k = tuningrows * tuningcols;
thetas          = linspace(0, 360, 361);
exampleids_untund = 1:k; % randsample([1: length(nonsigneurons)],k);

%%
% [~, exampleids_ori] = sort(max(R2_OriTun(:, orineurons)), 'descend');
[~, exampleids_ori] = sort(max(f_ori(2, :, orineurons)), 'descend');
exampleids_ori = exampleids_ori(:);

exampleids_ori_cont = [];
iter = 1;
for ii = 1: length(orineurons)
    nn = orineurons(exampleids_ori(ii));
%     if tcd.coef.aic_combo(nn)  < tcd.coef.aic_ind(nn)  
    if max(f_ori(2, :, nn)) > max(f_ori(1, :, nn))
        exampleids_ori_cont(iter) = exampleids_ori(ii);
        iter = iter+ 1;
    end
end
% exampleids_ori_cont = exampleids_ori;
exampleids_ori = exampleids_ori_cont(1:k); % randsample(exampleids_ori(1: k), k);

%%
% [~, exampleids_dir] = sort(max(R2_DirTun(:, dirneurons)), 'descend');
[~, exampleids_dir] = sort(max(f_dir(2, :, dirneurons)), 'descend');
exampleids_dir = exampleids_dir(:);

exampleids_dir_cont = [];
iter = 1;
for ii = 1: length(dirneurons)
    nn = dirneurons(exampleids_dir(ii));
%     if tcd.coef.aic_combo(nn)  < tcd.coef.aic_ind(nn)  
    if max(f_dir(2, :, nn)) > max(f_dir(1, :, nn))
        exampleids_dir_cont(iter) = exampleids_dir(ii);
        iter = iter+ 1;
    end
end
% exampleids_dir_cont = exampleids_dir;
exampleids_dir = exampleids_dir_cont(1:k); % randsample(exampleids_dir(1: k), k);

%%
if ~singlefig; figure; hold on; end
for jj = 1: k
    jr = mod(jj-1,tuningrows)+1;
    jc = ceil(jj/tuningrows);
    if singlefig
        subplotcm([(1+(jc-1)*2) (0.8+(tuningrows-jr+1)*1.5) 1.5 1]);  hold on;
    else, subplot(tuningrows,ceil(k/tuningrows), jj);
    end
    i = nonsigneurons(exampleids_untund(jj));
    for conid = 1: ncons
        plotntuning(sr, oris, cons, i, conid, thetas, ...
            vmfun.NullTun(tcd_sep.coef.x_null{conid}(i,:), thetas) / eventRateScale, uncol);
    end
    ylim([0 ylimmax]);
    set(gca,'Box','off','XTick',[0 max(uoris)],'XTickLabels',[],'YTick',[0 0.1]);
    if jc > 1, set(gca, 'YColor','none');  end
end
if ~singlefig; figure; hold on; end
for jj = 1: k
    jr = mod(jj-1,tuningrows)+1;
    jc = ceil(jj/tuningrows);
    if singlefig
        subplotcm([(1+(tuningcols+jc-1)*2) (0.8+(tuningrows-jr+1)*1.5) 1.5 1]);  hold on;
    else, subplot(tuningrows,ceil(k/tuningrows), jj);
    end
    i = orineurons(exampleids_ori(jj));
    [f_lo, f_hi] = fitMultip(x_OriTun(:, i), vmfun, thetas, thetas, 'OriTun');
%     f_hi = vmfun.OriTun_vonMises(tcd.coef.x_combo{i}(1: 4), thetas)/ eventRateScale;
%     f_lo = tcd.coef.x_combo{i}(5)+ tcd.coef.x_combo{i}(6)* f_hi;
    plotntuning(sr, oris, cons, i, 1, thetas, f_lo/eventRateScale, oricol);
    plotntuning(sr, oris, cons, i, 2, thetas, f_hi/eventRateScale, oricol);
    
    ylim([0 ylimmax]);
    set(gca,'Box','off','XTick',[0 max(uoris)],'XTickLabels',[],'YColor','none');
end
if ~singlefig; figure; hold on; end
for jj = 1: k
    jr = mod(jj-1,tuningrows)+1;
    jc = ceil(jj/tuningrows);
    if singlefig
        subplotcm([(1+(2*tuningcols+jc-1)*2) (0.8+(tuningrows-jr+1)*1.5) 1.5 1]);  hold on;
    else, subplot(tuningrows,ceil(k/tuningrows), jj);
    end
    i = dirneurons(exampleids_dir(jj));
    [f_lo, f_hi] = fitMultip(x_DirTun(:, i), vmfun, thetas, thetas, 'DirTun');
%     f_hi = vmfun.DirTun_vonMises(tcd.coef.x_combo{i}(1: 5), thetas)/ eventRateScale;
%     f_lo = tcd.coef.x_combo{i}(6)+ tcd.coef.x_combo{i}(7)* f_hi;
    plotntuning(sr, oris, cons, i, 1, thetas, f_lo/eventRateScale, dircol);
    plotntuning(sr, oris, cons, i, 2, thetas, f_hi/eventRateScale, dircol);
    
    ylim([0 ylimmax]);
    set(gca,'Box','off','XTick',[0 max(uoris)],'XTickLabels',[],'YColor','none');
end

%% plot adaptation
jittershift = [0, 10];
poris = [0; uoris];
contcols(:, :, 1) = oricols*0.4+[1 1 1]*0.6;
contcols(:, :, 2) = oricols;

if ~singlefig; figure; hold on; end
for jj = 1:tuningcols
    if singlefig
        subplotcm([(1+(jj-1)*2) 0.8 1.5 1]);  hold on;        
    else, subplot(2,1,jj);  hold on;
    end
    for conid = 1: ncons
        msr = oriresp(sr(cons == consu(conid), :), oris(cons == consu(conid)), nonsigneurons(exampleids_untund((jj-1)*tuningrows+1)));
        for ii = 1:(length(uoris)+1)
            if ii > length(uoris), c = [0 0 0]; lw = 1;
            else, c = contcols(ii,:, conid); lw = 0.5; end
            plot(poris+ jittershift(conid), [msr(ii,end) msr(ii,:)], '-', 'Colo', c, 'LineWidth', lw);
        end
    end
    xlim([-3 363]);  ylim([0 ymulimmax]);
    set(gca,'Box','off','XTick',[0 max(uoris)],'XTickLabels',{'0','360'},'YTick',[0 0.1]);
    if jj == 1, xlabel('drift direction');  else, set(gca, 'YColor','none'); end
end


if ~singlefig; figure; hold on; end

for jj = 1:2
    if singlefig
        subplotcm([(1+(tuningcols+jj-1)*2) 0.8 1.5 1]);  hold on;        
    else, subplot(2,1,jj);  hold on;
    end
    for conid = 1: ncons
        msr = oriresp(sr(cons == consu(conid), :), oris(cons == consu(conid)), orineurons(exampleids_ori((jj-1)*tuningrows+1)));
        for ii = 1:(length(uoris)+1)
            if ii > length(uoris), c = [0 0 0]; lw = 1;
            else, c = contcols(ii,:, conid); lw = 0.5; end
            plot(poris+ jittershift(conid), [msr(ii,end) msr(ii,:)], '-', 'Colo', c, 'LineWidth', lw);
        end
    end
    xlim([-3 363]);  ylim([0 ymulimmax]);
    set(gca,'Box','off','XTick',[0 max(uoris)],'XTickLabels',{'0','360'},'YColor','none');
end


if ~singlefig; figure; hold on; end

for jj = 1:2
    if singlefig
        subplotcm([(1+(tuningcols*2+jj-1)*2) 0.8 1.5 1]);  hold on;        
    else, subplot(2,1,jj);  hold on;
    end
    for conid = 1: ncons
        msr = oriresp(sr(cons == consu(conid), :), oris(cons == consu(conid)), dirneurons(exampleids_dir((jj-1)*tuningrows+1)));
        for ii = 1:(length(uoris)+1)
            if ii > length(uoris), c = [0 0 0]; lw = 1;
            else, c = contcols(ii,:, conid); lw = 0.5; end
            plot(poris+ jittershift(conid), [msr(ii,end) msr(ii,:)], '-', 'Colo', c, 'LineWidth', lw);
        end
    end
    xlim([-3 363]);  ylim([0 ymulimmax]);
    set(gca,'Box','off','XTick',[0 max(uoris)],'XTickLabels',{'0','360'},'YColor','none');
end



%% plot pairwise correlations and R2 across animals
animals = length(datasets);
fiwidth = 2;
ledge = 1.5+6*tuningcols;           % left edge of animals panel
fheight = 1.5*(tuningrows+1)-0.5;  % overall height across animals
fiheight = (fheight-(animals-1)*0.5)/(2*animals);  % height of per-animal panel

oricol_cont(1, :) = oricol* 0.5+ 0.5;
oricol_cont(2, :) = oricol;

dircol_cont(1, :) = dircol* 0.5+ 0.5;
dircol_cont(2, :) = dircol;

for animal = 1:animals
    sessions = length(datasets{animal});
    bedge = 0.8+fheight-animal*fiheight-(animal-1)*0.5;
    % load and process data
    c = cell(2,sessions);
    cshuf = cell(2,sessions);
    cun = cell(2,sessions);
    cori = cell(2,sessions);
    cdir = cell(2,sessions);
    R2ori = cell(2,sessions);
    R2dir = cell(2,sessions);
    clo   = cell(1,sessions);
    chi   = cell(1,sessions);
    for sessi = 1:sessions
        dsi = datasets{animal}{sessi};
        [sr, oris, cons] = loaddata(dsi,subsample);
        
        tcfitfile_sep = [tcfitpath filesep dsi filesep 'TuningCoef_Subsample_' subsample '.mat'];
        fprintf('Loading individual fits %s ...\n', tcfitfile_sep);
        tcd_sep = load(tcfitfile_sep);
        
        tcfitfile = [tcfitpath filesep dsi filesep 'TuningCoef_Combo_Subsample_' subsample '.mat'];
        uoris = unique(oris);
        ucons = unique(cons);
        fprintf('Loading %s ...\n', tcfitfile);
        tcd = load(tcfitfile);
        
        [clo{1, sessi}, chi{1, sessi}] = pairwisecorr_cont(sr, cons);
        
        conid = 1;

        vmfun.OriTun_vonMises   = tcd_sep.coef.OriTun_vonMises{conid};
        vmfun.DirTun_vonMises   = tcd_sep.coef.DirTun_vonMises{conid};
        vmfun.NullTun           = tcd_sep.coef.NullFunc{conid};

        %% select neuron tuning based on high-contrast stim
        conid = 2; 

        nonsigneurons = tcd_sep.coef.alpha_OriVsNull{conid} >= 0.05;
        orineurons = find(tcd_sep.coef.alpha_OriVsNull{conid} < 0.05 & ...
            tcd_sep.coef.alpha_OriVsDir{conid} > 0.05);
        dirneurons = find(tcd_sep.coef.alpha_OriVsNull{conid} < 0.05 & ...
            tcd_sep.coef.alpha_OriVsDir{conid} < 0.05);
        
        %%
        orientationForTrials = oris(cons == ucons(conid));
        averageSpkRateOverTrial = sr(cons == ucons(conid), :); 
        nNeurons = size(averageSpkRateOverTrial, 2);
        R2_OriTun = NaN(ncons, nNeurons);
        R2_DirTun = NaN(ncons, nNeurons);
        f = NaN(1, length(uoris));
        f_ori = NaN(ncons, length(uoris), nNeurons);
        f_dir = NaN(ncons, length(uoris), nNeurons);
        
        x_OriTun = nan(6, nNeurons);
        for ii = 1: length(orineurons)
            x_OriTun(:, orineurons(ii))    = tcd.coef.x_combo{orineurons(ii), 1};
        end

        x_DirTun = nan(7, nNeurons);
        for ii = 1: length(dirneurons)
            x_DirTun(:, dirneurons(ii))    = tcd.coef.x_combo{dirneurons(ii), 1};
        end

%         for con = 1: ncons
%             x_null(:, :, con)      = tcd_sep.coef.x_null{con,1};
%         end
        
        for nn = 1: nNeurons

            if sum(dirneurons == nn)

                [mu_lo, mu_hi] = fitMultip(x_DirTun(:, nn), vmfun, uoris, uoris, 'DirTun');
                f_dir(1, :, nn) = mu_lo/ eventRateScale;
                f_dir(2, :, nn) = mu_hi/ eventRateScale;

            elseif sum(orineurons == nn)
                [mu_lo, mu_hi] = fitMultip(x_OriTun(:, nn), vmfun, uoris, uoris, 'OriTun');
                f_ori(1, :, nn) = mu_lo/ eventRateScale;
                f_ori(2, :, nn) = mu_hi/ eventRateScale;

            end
            for conid = 1: 2
                for orii = 1:length(uoris)
                    f(orii) = mean(sr(oris == uoris(orii) & cons == ucons(conid), nn));
                end
                R2_OriTun(conid, nn) = 1 - sum((f_ori(conid, :, nn) - f).^2) / sum((f - mean(f)).^2);
                R2_DirTun(conid, nn) = 1 - sum((f_dir(conid, :, nn) - f).^2) / sum((f - mean(f)).^2);
            end
        end
        
        for conid = 1: 2
            c{conid, sessi} = pairwisecorr(sr(cons == ucons(conid), :));
            cshuf{conid, sessi} = pairwisecorr(shufsr(sr(cons == ucons(conid), :), oris(cons == ucons(conid)), cons(cons == ucons(conid))));
            cun{conid, sessi} = pairwisecorr(sr(cons == ucons(conid), nonsigneurons));
            cori{conid, sessi} = pairwisecorr(sr(cons == ucons(conid), orineurons));
            cdir{conid, sessi} = pairwisecorr(sr(cons == ucons(conid), dirneurons));
            R2ori{conid, sessi} = R2_OriTun(conid, orineurons);
            R2dir{conid, sessi} = R2_DirTun(conid, dirneurons);
        end            
    end

    % pairwise correlations, low-low & hig-high contrast
    if singlefig
        subplotcm([ledge+0.5+fiwidth bedge fiwidth fiheight]);
        hold on;
    else
        figure;  hold on;
    end
    for sessi = 1:sessions
        if strcmp(dataset, datasets{animal}{sessi}), lw = 1; else, lw = 0.5; end
        [f, x] = subsampledecdf(clo{sessi}, cdfsubsamples);
        stairs([0 x' 1],[0 f' 1],'Color',uncol,'LineWidth',lw);
        [f, x] = subsampledecdf(chi{sessi}, cdfsubsamples);
        stairs([0 x' 1],[0 f' 1],'Color',dircol,'LineWidth',lw);
        xlim([-0.2 0.4]);  ylim([0 1]);
        plot([1 1]*mean(clo{sessi}), ylim, '-', 'LineWidth', 0.5, 'Color', 0.5+0.5*uncol);
        plot([1 1]*mean(chi{sessi}), ylim, '-', 'LineWidth', 0.5, 'Color', 0.5+0.5*dircol);
    end
    set(gca,'Box','off','XTick',[-0.2 0 0.2 0.4],'YTick',[0 1]);
    if animal == animals, xlabel('correlation coefficient \rho');
    else set(gca,'XTickLabels',[]); end
    
    % R2 for fitted neurons
    if singlefig
        subplotcm([ledge+1+2*fiwidth bedge fiwidth fiheight]);
        hold on;
    else
        figure;  hold on;
    end
    for sessi = 1:sessions
        for conid = 1: 2
            if strcmp(dataset, datasets{animal}{sessi}), lw = 1; else, lw = 0.5; end
            [f,x] = ecdf(R2ori{conid, sessi});
            stairs([0 x' 1],[0 f' 1],'Color',oricol_cont(conid, :),'LineWidth',lw);
            [f,x] = ecdf(R2dir{conid, sessi});
            stairs([0 x' 1],[0 f' 1],'Color',dircol_cont(conid, :),'LineWidth',lw);
            plot([1 1]*mean(R2ori{conid, sessi}), ylim, '-', 'LineWidth', 0.5, 'Color', 0.5+0.5*oricol_cont(conid, :));
            plot([1 1]*mean(R2dir{conid, sessi}), ylim, '-', 'LineWidth', 0.5, 'Color', 0.5+0.5*dircol_cont(conid, :));
        end
    end
    xlim([0 1]);  ylim([0 1]);        
    set(gca,'Box','off','XTick',[0 0.25 0.5 0.75 1],'YTick',[0 1]);
    if animal == animals, xlabel('R^2');
    else set(gca,'XTickLabels',[]); end
end


%% write plot to file
if singlefig
    fprintf('\nWriting figure to figS12.pdf\n');
    print(['figs' filesep 'figS12'], '-dpdf');
end


function plotntuning(sr, oris, cons, neurid, conid, thetas, fthetas, col)
%% plots neural responses as well as fitted tuning curve
consu = unique(cons);
jittersd = 1;
if conid == 1 
    jittershift = -5;
    col = col*0.5+[1 1 1]*0.5;
%     col = 0.9* col; 
elseif conid == 2
    jittershift = 5;
%     shadedcol = col*0.5+[1 1 1]*0.5;
%     col = 0.5* col; 
end

shadedcol = col*0.8+[1 1 1]*0.2;

sr = sr(cons == consu(conid), :);
oris = oris(cons == consu(conid));

[sr_ori_ave, ~, sr_ori_neg, sr_ori_pos] = errbar_util(sr(:, neurid), oris);

orisu = [0 unique(oris)'];
hold on;
for i = 1:length(orisu)
    j = mod(oris, 360) == mod(orisu(i), 360);
    plot(orisu(i) + jittersd * randn(1, sum(j))+ jittershift, sr(j,neurid)','o',...
        'MarkerFaceColor',shadedcol, ...
        'MarkerEdgeColor','None','MarkerSize',0.5);  alpha(0.1);
    if i == 1, j = length(orisu) - 1; else j = i-1; end
    plot((orisu(i)+ jittershift)*[1 1], [sr_ori_neg(j) sr_ori_pos(j)], '-', ...
        'Color', col, 'LineWidth', 1.4);
end

plot(orisu+ jittershift, [sr_ori_ave(end); sr_ori_ave], 'o', 'MarkerSize', 3, ...
    'MarkerFaceColor', col, 'MarkerEdgeColor', 'None');
plot(thetas+ jittershift, fthetas, '-', 'LineWidth', 0.5, 'Color', shadedcol);
xlim([(min(thetas)-6*jittersd) (max(thetas)+6*jittersd)]);

    
function [sr_ori_ave, sr_ori_std, sr_ori_neg, sr_ori_pos] = errbar_util(sr, ori)
uoris           = sort(unique(ori), 'ascend');
orin            = length(uoris);
sr_ori_ave      = nan(orin, 1);
sr_ori_std      = nan(orin, 1);
sr_ori_neg      = nan(orin, 1);
sr_ori_pos      = nan(orin, 1);    
for ii = 1: orin
    sr_ori_tmp  = sr(ori == uoris(ii));
    sr_ori_ave(ii)  = mean(sr_ori_tmp);
    sr_ori_std(ii)  = std(sr_ori_tmp);
    sr_ori_neg(ii)  = prctile(sr_ori_tmp, 25);
    sr_ori_pos(ii)  = prctile(sr_ori_tmp, 75);
end


function sr = shufsr(sr, oris, cons)
%% trial-shuffle data within each condition
N = size(sr, 2);
ucons = unique(cons);  conn = length(ucons);
uoris = unique(oris);  orin = length(uoris);
for coni = 1:conn
    for orii = 1:orin
        condtrials = find((oris == uoris(orii)) & (cons == ucons(coni)));
        T = length(condtrials);
        for n = 2:N
            sr(condtrials,n) = sr(condtrials(randperm(T)),n);
        end
    end
end


function [f,x] = subsampledecdf(x, subsamples)
if length(x) <= subsamples
    [f, x] = ecdf(x);
else
    [f, x] = ecdf(x(randperm(length(x), subsamples)));
end


function c = pairwisecorr(sr)
c = corr(sr); 
c = c(logical(triu(ones(size(sr, 2)), 1)));

function [clo, chi] = pairwisecorr_cont(sr, cons)
conus = unique(cons);
sr_lo = sr(cons == conus(1), :);
sr_hi = sr(cons == conus(2), :);

corr_lo = corr(sr_lo); 
mask = ~~triu(ones(size(corr_lo)), 1); 
clo = corr_lo(mask);

corr_hi = corr(sr_hi); 
chi = corr_hi(mask);


function msr = oriresp(sr, oris, n)
%% computes mean responses for neuron n conditional on previous orientation

uoris = unique(oris);  orin = length(uoris);
msr = NaN(orin + 1, orin);
prevoris = [-1; oris(1:(end-1))];
for i = 1:(orin+1)  % prev. orientation
    if i > orin, previ = true(size(sr,1),1);  % last i across all trials
    else, previ = prevoris == uoris(i);  end
    for j = 1:orin  % curr. orientation
        msr(i,j) = mean(sr(oris == uoris(j) & previ, n));
    end
end

function [mu_lo, mu_hi] = fitMultip(P, vmfun, theta_lo, theta_hi, Tuning)
           % model with 4 or 5 parameters
if strcmp(Tuning, 'DirTun')
    mu_hi = vmfun.DirTun_vonMises(P(1: 5), theta_hi);
    mu_lo = P(6)+ P(7)* vmfun.DirTun_vonMises(P(1: 5), theta_lo);
elseif strcmp(Tuning, 'OriTun')
    mu_hi = vmfun.OriTun_vonMises(P(1: 4), theta_hi);
    mu_lo = P(5)+ P(6)* vmfun.OriTun_vonMises(P(1: 4), theta_lo);
elseif strcmp(Tuning, 'null')
    mu_hi = vmfun.NullTun(P(1), theta_hi);
    mu_lo = P(2)+ P(3)* vmfun.NullTun(P(12), theta_lo);          
end