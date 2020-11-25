function plotFig2
%% plots part of Fig. 2

%% (general) settings
singlefig = true;   % set to true to have all planels in single figure
dataset = 'm25b';
datafile = 'm25_170523';  % data file corresponding to dataset
conid = 1;  % contrast to use
tcfitpath = ['.' filesep 'tuning_fits'];
datapath = ['.' filesep 'data'];
bkgminc = 0;     % background image intensity
bkgmaxc = 1;
bkgpowc = 0.5;   % [0 1] -> [0 1]^bkgpowc
imagesize = 680;  % side with of image in micrometers
tuningn = 20;
exampleids = [1 3 6];  % example neurons [non-tuned ori dir]
eventRateScale  = 30;
uncol = [27 158 119] / 255;
dircol = [217 95 2] / 255;
oricol = [117 112 179] / 255;


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
    figure('Name', 'Figure', 'Units', 'centimeters', 'Position', [0 0 10.8 5.1]);
end


%% plot image with color-labeled cells
neurfile = [datapath filesep datafile '.mat'];
bkgfile = [datapath filesep datafile '_refIm.mat'];
tcfitfile = [tcfitpath filesep dataset filesep 'TuningCoef_Subsample_none.mat'];
fprintf('Loading %s ...\n', neurfile);
neurd = load(neurfile);
fprintf('Loading %s ...\n', bkgfile);
bkgd = load(bkgfile);
i = fields(bkgd);
bkgd = bkgd.(i{1});
fprintf('Loading %s ...\n', tcfitfile);
tcd = load(tcfitfile);

% normalize background to [bkgminc bkcmaxc]
bkgsize = size(bkgd);
bkgd = reshape(bkgd,[],1);  % turn into vector
bkgmin = min(bkgd);
bkgmax = max(bkgd);
nimage = repmat(...
    bkgminc+((bkgd-bkgmin)/(bkgmax-bkgmin)).^bkgpowc * (bkgmaxc-bkgminc), ...
    1, 3);  % add RBG dimension

% normalize cell shapes
nshapes = neurd.stimExpt.cellFilts;
shapemax = max(max(nshapes));

% add individual neurons
signeurons = find(tcd.coef.alpha_OriVsNull{conid} < 0.05);
oripref = mod(tcd.coef.x_OriTun{conid}(signeurons, 4), 180);
for i = 1:length(signeurons)
    ic = hsv2rgb([(oripref(i)/180) 1 1]);  % ori = hue in HSV space
    nshape = nshapes(:,signeurons(i)) / shapemax;
    nimage(nshape > 0,:) = bsxfun(@times, nimage(nshape > 0,:), ic);
end
% back from vector x RGB to 3D matrix
nimage = reshape(nimage, bkgsize(1), bkgsize(2),3);

if singlefig
    subplotcm([1 0.8 3.5 3.5]);  hold on;
else
    figure; hold on;
end
image(linspace(0,1,bkgsize(1)), linspace(0,1,bkgsize(2)), nimage);
xlim([0 1]);  ylim([0 1]);
% bar corresponding to 100 micrometers
plot(0.95-[0 (100/imagesize)], [1 1]*0.05, '-', 'LineWidth', 1, 'Color', [1 1 1]);
text(0.95-50/imagesize, 0.05, '100 \mu m', 'Color', [1 1 1], ...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
% colormap
hsvmap = [linspace(0,1,100)' ones(100,1) ones(100,1)];
colormap(hsv2rgb(hsvmap));
colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',...
    {'0^\circ','45^\circ','90^\circ','145^\circ','180^\circ'});
set(gca,'Box','off','XColor','none','YColor','none','DataAspectRatio', [1 1 1]);


%% plot tuning curves
if singlefig
    subplotcm([5.5 0.8 3.5 2.9]);  hold on;
else
    figure; hold on;
end
thetas = linspace(0, 360, 361);
OriTun_vonMises = tcd.coef.OriTun_vonMises{conid};
DirTun_vonMises = tcd.coef.DirTun_vonMises{conid};
NullTun = tcd.coef.NullFunc{conid};
plotneurons = signeurons(round(linspace(1,length(signeurons),tuningn)));
for i = 1:tuningn
    j = plotneurons(i);
    if tcd.coef.alpha_OriVsDir{conid}(j) < 0.05
        % direction-tuned neuron
        plot(thetas, ...
            DirTun_vonMises(tcd.coef.x_DirTun{conid}(j,:), thetas) / eventRateScale, ...
            '-', 'Color', [1 1 1]*0.2, 'LineWidth', 0.25);
    else
        % orientation-tuned neuron
        plot(thetas, ...
            OriTun_vonMises(tcd.coef.x_OriTun{conid}(j,:), thetas) / eventRateScale, ...
            '-', 'Color', [1 1 1]*0.2, 'LineWidth', 0.25);
    end
end
xlim([0 360]);  ylim([0 0.25]);
set(gca,'Box','off','XTick',[0 90 180 270 360],...
    'XTickLabel',{'0^\circ','90^\circ','180^\circ','270^\circ','360^\circ'},...
    'YTick',[0 0.1 0.2]);
xlabel('drift direction');
ylabel('\Delta F / F');


%% plot examples for each tuning
nonsigneurons = find(tcd.coef.alpha_OriVsNull{conid} >= 0.05);
orineurons = find(tcd.coef.alpha_OriVsNull{conid} < 0.05 & ...
    tcd.coef.alpha_OriVsDir{conid} > 0.05);
dirneurons = find(tcd.coef.alpha_OriVsNull{conid} < 0.05 & ...
    tcd.coef.alpha_OriVsDir{conid} < 0.05);

if singlefig
    subplotcm([9.2 2.8 1.1 0.9]);  hold on;
else
    figure; hold on;
end
i = nonsigneurons(exampleids(1));
plotntuning(neurd, i, conid, thetas, ...
    NullTun(tcd.coef.x_null{conid}(i,:), thetas) / eventRateScale, uncol);
ylim([0 0.3]);

if singlefig
    subplotcm([9.2 1.8 1.1 0.9]);  hold on;
else
    figure; hold on;
end
i = orineurons(exampleids(2));
plotntuning(neurd, i, conid, thetas, ...
    OriTun_vonMises(tcd.coef.x_OriTun{conid}(i,:), thetas) / eventRateScale, oricol);
ylim([0 0.3]);

if singlefig
    subplotcm([9.2 0.8 1.1 0.9]);  hold on;
else
    figure; hold on;
end
i = dirneurons(exampleids(3));
plotntuning(neurd, i, conid, thetas, ...
    DirTun_vonMises(tcd.coef.x_DirTun{conid}(i,:), thetas) / eventRateScale, dircol);
ylim([0 0.3]);


%% save figure
if singlefig
    fprintf('\nWriting figure to fig2.pdf\n');
    print(['figs' filesep 'fig2'], '-dpdf');
end



function plotntuning(d, neurid, conid, thetas, fthetas, col)
%% plots neural responses as well as fitted tuning curve

jittersd = 2.5;
shadecol = col*0.8 + [1 1 1]*0.2;

cons = unique(d.visCon);
oris = [0 unique(d.visOri)'];
f = NaN(length(oris), 3);
for i = 1:length(oris)
    j = d.visCon == cons(conid) & mod(d.visOri, 360) == mod(oris(i), 360);
    plot(oris(i) + jittersd * randn(1, sum(j)), d.deResp(j,neurid)','o',...
        'MarkerFaceColor',shadecol, 'MarkerEdgeColor','none','MarkerSize',0.5); alpha(0.1);
    f(i,:) = [mean(d.deResp(j,neurid)) ...
        prctile(d.deResp(j,neurid),25) prctile(d.deResp(j,neurid),75)];
    plot(oris(i)*[1 1], f(i,2:3), '-', 'Color', col, 'LineWidth', 1);
end
plot(thetas, fthetas, '-', 'LineWidth', 0.5, 'Color', shadecol);
plot(oris, f(:,1), 'o', 'MarkerFaceColor', col, 'MarkerEdgeColor', 'none', ...
    'MarkerSize', 2);
xlim([(min(thetas)-3*jittersd) (max(thetas)+3*jittersd)]);
set(gca,'Box','off','YColor','none','XTick',[0 90 180 270 360],'XTickLabel',[]);
