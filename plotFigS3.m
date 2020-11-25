function plotFigS3
%% plots SI figure showing per-neuron tuning and response time-course

momentpath = ['.' filesep 'moment_cache'];

%% (general) settings
singlefig = true;   % set to true to have all planels in single figure
tcfitpath = ['.' filesep 'tuning_fits'];
datapath = ['.' filesep 'data'];

eventRateScale = 30;
uncol = [27 158 119] / 255;
dircol = [217 95 2] / 255;
oricol = [117 112 179] / 255;
plotxdist = 6.5;
plotydist = 5.5;
plotns = 3;

% AJ high contrast
datafile = [datapath filesep 'AJ060_190902.mat'];
dataset = 'aj60a';
conid = 2;
maxdFs = [1 0.75 1; ... % direction-tuned
          2 4 2; ... % orientation-tuned
          2 2 2];    % untuned
exampleids_dir = [47 14 43];
exampleids_ori = [128 18 144];
Fscale = 1;


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
    figure('Name', 'Figure', 'Units', 'centimeters', 'Position', ...
           [0 0 plotns*plotxdist 4.2*plotydist]);
end


%% loading data
% activity data
fprintf('Loading %s ...\n', datafile);
d = load(datafile);
if length(fields(d)) == 1
    f = fields(d);
    d = d.(f{1});
    datatype = 2;
else
    datatype = 1;
end
if exist('datafile_dF','var')   % replace df with data in file
    fprintf('Loading %s ...\n', datafile_dF);
    d_dF = load(datafile_dF);
    f = fields(d_dF);
    if length(f) == 1 && ~strcmp(f{1}, 'dF'), d_dF = d_dF.(f{1}); end
    d.stimExpt.dF = d_dF.dF;
    clear d_dF
end
fprintf('Loading %s ...\n', dataset');
[sr, oris, cons] = loaddata(dataset);
uoris = unique(oris);
ucons = unique(cons);
% tuning curve fits
tcfitfile = [tcfitpath filesep dataset filesep 'TuningCoef_Subsample_none.mat'];
fprintf('Loading %s ...\n', tcfitfile);
tcd = load(tcfitfile);
OriTun_vonMises = tcd.coef.OriTun_vonMises{conid};
DirTun_vonMises = tcd.coef.DirTun_vonMises{conid};
NullTun         = tcd.coef.NullFunc{conid};
x_OriTun    = tcd.coef.x_OriTun{conid,1};
x_DirTun    = tcd.coef.x_DirTun{conid,1};
orientationForTrials = oris;
averageSpkRateOverTrial = sr; 
nNeurons = size(averageSpkRateOverTrial, 2);
% compute R2 of fits
R2_OriTun = NaN(1, nNeurons);
R2_DirTun = NaN(1, nNeurons);
f = NaN(1, length(uoris));
for nn = 1: nNeurons
    for orii = 1:length(uoris)
        f(orii) = mean(sr(oris == uoris(orii) & cons == ucons(conid), nn));
    end
    f_ori = OriTun_vonMises(x_OriTun(nn, :), uoris)' / eventRateScale;
    f_dir = DirTun_vonMises(x_DirTun(nn, :), uoris)' / eventRateScale;
    R2_OriTun(nn) = 1 - sum((f_ori - f).^2) / sum((f - mean(f)).^2);
    R2_DirTun(nn) = 1 - sum((f_dir - f).^2) / sum((f - mean(f)).^2);
    R2_OriTun(nn)   = corr(OriTun_vonMises(x_OriTun(nn, :), orientationForTrials),  averageSpkRateOverTrial(:, nn))^2;
    R2_DirTun(nn)   = corr(DirTun_vonMises(x_DirTun(nn, :), orientationForTrials), averageSpkRateOverTrial(:, nn))^2;
end
% separate neurons into categories by fit significance
nonsigneurons = find(tcd.coef.alpha_OriVsNull{conid} >= 0.05);
orineurons = find(tcd.coef.alpha_OriVsNull{conid} < 0.05 & ...
    tcd.coef.alpha_OriVsDir{conid} > 0.05);
dirneurons = find(tcd.coef.alpha_OriVsNull{conid} < 0.05 & ...
    tcd.coef.alpha_OriVsDir{conid} < 0.05);


%% plot activity / tuning curves of untuned and tuned neurons for example data
if ~exist('exampleids_untund','var')
    exampleids_untund = 1:plotns; % randsample([1: length(nonsigneurons)],k);
end
if ~exist('exampleids_ori','var')
    [~, exampleids_ori] = sort(R2_OriTun(orineurons), 'descend');
    exampleids_ori = exampleids_ori(1:plotns); % randsample(exampleids_ori(1: k), k);
end
if ~exist('exampleids_dir','var')
    [~, exampleids_dir] = sort(R2_DirTun(dirneurons), 'descend');
    exampleids_dir = exampleids_dir(1:plotns); % randsample(exampleids_dir(1: k), k);
end


%% find optimal decoder directions for dori=45 pairs
optdec = NaN(nNeurons, 4);
fprintf('Computing optimal decoders...\n');
for orii = 1:4
    ori1 = orii*2-1;
    ori2 = orii*2;
    momentfile = sprintf('%s_o%d-%d_c1.mat', dataset, ori1, ori2);
    fprintf('Loading %s ...\n', [momentpath filesep momentfile]);
    dp = load([momentpath filesep momentfile]);
    fsig = pinv(dp.S)* dp.mu'/ (dp.mu* pinv(dp.S)* dp.mu');
    optdec(:,orii) = fsig/norm(fsig); 
end

%% direction-tuned neurons
[stimFrames, oris, cons] = getStimFrames(d.stimExpt, datatype);
fps = 30;
plotframes = round(-0.25*fps):round(1.5*fps);

for ni = 1:plotns
    nid = dirneurons(exampleids_dir(ni));
    dF = d.stimExpt.dF(nid,:);
    plotRaw([plotxdist*(ni-0.5) plotydist*1.5], ...
        dF, oris, cons, conid, stimFrames, maxdFs(1, ni), Fscale, dircol);    
    plotProj([plotxdist*(ni-0.5) plotydist*0.5], ...
        d.stimExpt.dF, nid, optdec, oris, cons, conid, stimFrames, maxdFs(1, ni), Fscale, dircol);
end

%% orientation-tuned neurons
for ni = 1:plotns
    nid = orineurons(exampleids_ori(ni));
    dF = d.stimExpt.dF(nid,:);
    plotRaw([plotxdist*(ni-0.5) plotydist*3.7], ...
        dF, oris, cons, conid, stimFrames, maxdFs(2, ni), Fscale, oricol);
    plotProj([plotxdist*(ni-0.5) plotydist*2.7], ...
        d.stimExpt.dF, nid, optdec, oris, cons, conid, stimFrames, maxdFs(2, ni), Fscale, oricol);
end

%% untuned neurons
% for ni = 1:plotns
%     nid = nonsigneurons(exampleids_untund(ni));
%     dF = d.stimExpt.dF(nid,:);
%     plotRaw([plotxdist*(ni-0.5) plotydist*5.9], ...
%         dF, oris, cons, conid, stimFrames, maxdFs(3, ni), uncol);
%     plotProj([plotxdist*(ni-0.5) plotydist*4.], ...
%         d.stimExpt.dF, nid, optdec, oris, cons, conid, stimFrames, maxdFs(3, ni), uncol);
% end

% save figure
if singlefig
    fprintf('\nWriting figure to figS3.pdf\n');
    print(['figs' filesep 'figS3'], '-dpdf');
end


function plotRaw(pos, dF, oris, cons, conid, stimFrames, maxdF, Fscale, col)

fps = 30;
Ksig = 0.025;
nsamples = 200;
figh = 1;  figw = 1.5;
% 45 is direction northwest, 90 is north, 135 is northeast, etc.
centers = figh * [-1.6 1.1; 0 1.9; 1.6 1.1; 2.2 0; 1.6 -1.1; 0 -1.9; -1.6 -1.1; -2.2 0]; 

plotframes = round(-0.25*fps):round(1.5*fps);
% gaussian smoothing kerenel
K = exp(- bsxfun(@minus, 1:length(plotframes), (1:length(plotframes))').^2 ...
    / (2*(Ksig*fps)^2));
K = bsxfun(@rdivide, K, sum(K, 2));

uoris = unique(oris);
ucons = unique(cons);

for orii = 1:length(uoris)
    % collect activity for particular orientation
    oritrials = oris == uoris(orii) & cons == ucons(conid);
    oriframes = bsxfun(@plus, stimFrames(oritrials), plotframes);
    oridF = dF(oriframes);
    subplotcm([pos(1)+ centers(orii, 1)-figw/2  ...
               pos(2)+centers(orii, 2)-figh/2 figw figh]);  hold on;
    patch([0 0 0.5 0.5], [0 maxdF maxdF 0], [1 1 1]*0.9, 'EdgeColor', 'none');
    sampletrials = round(linspace(1,size(oridF,1),nsamples));
    dFavg   = mean(oridF,1);
    l2nrm = sum((oridF- dFavg).^2, 2);
    l2nrm = l2nrm./std(l2nrm);
    l2nrm(l2nrm > 0.99) = 0.99;
    for j = sampletrials
        plot(plotframes / fps, K * oridF(j,:)', '-', 'Color', (1- l2nrm(j))*[col 0.08]);
    end
    plot(plotframes / fps, K * dFavg', '-', 'LineWidth', 1, 'Color', col);
    xlim([plotframes(1) plotframes(end)] / fps);
    ylim([0 maxdF]); xlim([-0.25, 1.5]);
    set(gca,'Box','off','XColor','none','YColor','none','XTick',[],'YTick',[]);
    if orii == 8, plot(-[1 1]*0.25, [0 Fscale], 'k-'); end  % scale bar
end


function plotProj(pos, dF, nid, optdec, oris, cons, conid, stimFrames, maxdF, Fscale, col)

fps = 30;
Ksig = 0.025;
nsamples = 200;
figh = 1;  figw = 1.5;
% 45 is direction northwest, 90 is north, 135 is northeast, etc.
centers = figh * [-1.6 1.1; 0 1.9; 1.6 1.1; 2.2 0; 1.6 -1.1; 0 -1.9; -1.6 -1.1; -2.2 0]; 

plotframes = round(-0.25*fps):round(1.5*fps);
% gaussian smoothing kerenel
K = exp(- bsxfun(@minus, 1:length(plotframes), (1:length(plotframes))').^2 ...
    / (2*(Ksig*fps)^2));
K = bsxfun(@rdivide, K, sum(K, 2));

uoris = unique(oris);
ucons = unique(cons);

for orii = 1:length(uoris)
    % collect activity for particular orientation
    oritrials = oris == uoris(orii) & cons == ucons(conid);
    oriframes = bsxfun(@plus, stimFrames(oritrials), plotframes);
    oridF = NaN(size(dF,1), size(oriframes,2), size(oriframes,1));
    for i = 1:size(dF,1)
        dFi = dF(i,:);
        oridF(i,:,:) = dFi(oriframes)';  % neuron x frame x trial
    end
    avgdF = mean(oridF, 3);  % neuron x frame
    w = optdec(:,floor((orii+1)/2));
    % distance along decoder, d = w^T (dF - <dF>)
    projdF = sum(bsxfun(@times, oridF - avgdF, w), 1); % 1 x frame x trial
    % map back into dF space, <dF> + w d
    projdF = avgdF + bsxfun(@times, w, projdF);
    % pick out neuron of interest
    projdF = squeeze(projdF(nid,:,:))';  % trial x frame
    avgdF = avgdF(nid,:);                % 1 x frame
    
    subplotcm([pos(1)+ centers(orii, 1)-figw/2  ...
               pos(2)+centers(orii, 2)-figh/2 figw figh]);  hold on;
    patch([0 0 0.5 0.5], [0 maxdF maxdF 0], [1 1 1]*0.9, 'EdgeColor', 'none');
    
    sampletrials = round(linspace(1,size(projdF,1),nsamples));
    l2nrm = sum((projdF - avgdF).^2, 2);
    l2nrm = l2nrm./std(l2nrm);
    l2nrm(l2nrm > 0.99) = 0.99;
    for j = sampletrials
        plot(plotframes / fps, K * projdF(j,:)', '-', 'Color', (1- l2nrm(j))*[col 0.08]);
    end
    plot(plotframes / fps, K * avgdF', '-', 'LineWidth', 1, 'Color', col);
    xlim([plotframes(1) plotframes(end)] / fps);
    ylim([0 maxdF]); xlim([-0.25, 1.5]);
    set(gca,'Box','off','XColor','none','YColor','none','XTick',[],'YTick',[]);
    if orii == 8, plot(-[1 1]*0.25, [0 Fscale], 'k-'); end  % scale bar
end


function [stimFrames, oris, cons] = getStimFrames(stimExpt, datatype)
%% returns stimulus frames, associated orientations and contrasts
% datatype = 1 for mxx data (default), 2 for ajxx data

if nargin < 2, datatype = 1; end
oris = [];  cons = [];  stimFrames = [];
for nBlock = find(stimExpt.stimBlocks)
    if datatype == 1
        % mxx data
        blockOffsetFrame = length(cat(1,stimExpt.frameTimes{1:nBlock-1}));
        lastTrial = find(~isnan(stimExpt.psych2frame{nBlock}),1,'last')-1;
        thisStimFrames = stimExpt.psych2frame{nBlock}(1:lastTrial) + blockOffsetFrame;
        thisOris = stimExpt.stimInfo{nBlock}((1:3:lastTrial*3)+4); % stimSynch grating angle
        thisCons = stimExpt.stimInfo{nBlock}((2:3:lastTrial*3)+4); % stimSynch Contrast
    else
        % ajxx data
        blockOffsetFrame = length(cat(2,stimExpt.frameTimes{1:nBlock-1}));
        lastTrial = length(stimExpt.psych2frame{nBlock});
        thisStimFrames = (stimExpt.psych2frame{nBlock}(1:lastTrial) + blockOffsetFrame)';
        if size(stimExpt.stimInfo, 1) == 2  % 2 constrast data
            thisOris = stimExpt.stimInfo{1,nBlock};
            thisCons = stimExpt.stimInfo{2,nBlock};
        else                                % single contrast data
            thisOris = stimExpt.stimInfo{nBlock};
            thisCons = zeros(lastTrial,1) + 0.1;
        end
    end
   
    stimFrames = cat(1, stimFrames, thisStimFrames);
    oris = cat(1, oris, thisOris);
    cons = cat(1, cons, thisCons);
end
