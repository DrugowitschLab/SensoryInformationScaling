function plotFig8
%% plots panels of figure 8

%% (general) settings
singlefig = true;   % set to true to have all planels in single figure
dataset = 'm25b';
datasets = {'m25a', 'm25b', 'm26a', 'm26b', ...
    'aj42a', 'aj42b', 'aj42c', 'aj42d', 'aj42e', ...
    'aj43a', 'aj43b', 'aj43c', 'aj43d', 'aj43e', 'aj43f', 'aj43g'};
%datasets = {'m25a', 'm25b', 'm26a', 'm26b'};
doriid = 1;  % focus on 45deg
doriexample = 4;  % 135 vs. 180 deg, as in asym info scaling figure
varcol = [0.12 0.47 0.71];
aligncol = [0.2 0.63 0.17];
infcol = [0 0 0];
dot90col = [0.89 0.1 0.11];
dot90size = 2;
dorisnonoverlap = [1 3 5 7];
cvn = 10;


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
    figure('Name', 'Figure', 'Units', 'centimeters', 'Position', [0 0 8.5 8.3]);
end


%% load example data and collect statistics
[sr, ori, con] = loaddata(dataset, 'none');
cons = unique(con);
assert(length(cons) == 1);  %% only deal with single-constrat data
N = size(sr, 2);
s = dataInfo(dataset);
oricomb = s.oricomb(1:2, s.oricomb(3,:) == s.doris(doriid));
noricomb = size(oricomb,2);
% compute various statistics
varn = NaN(noricomb, N);
varfracn = NaN(noricomb, N);
alignn = NaN(noricomb, N);
Ifracn = NaN(noricomb, N);
fprintf('Cross-validated statistics');
for combi = 1:noricomb
    fprintf(' %d', combi);
    ori1 = s.oris(oricomb(1,combi));
    ori2 = s.oris(oricomb(2,combi));
    [varni, alignni, Ini] = cvstats(sr, ori, con, cons(1), ori1, ori2, cvn);
    varn(combi, :) = mean(varni, 1);
    varfracn(combi, :) = mean(...
        bsxfun(@rdivide, cumsum(varni, 2), sum(varni, 2)), 1);
    alignn(combi, :)= mean(alignni, 1);
    Ifracn(combi, :) = mean(bsxfun(@rdivide, Ini, Ini(:,end)), 1);
end
fprintf('\n');


%% plot example variance per dimension
if singlefig
    subplotcm([1 4.8 3 3]);  hold on;
    text(-1,3.5,'b','Units','centimeters','FontWeight','bold',...
         'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
var90 = NaN(2, noricomb);
for combi = 1:noricomb
    i = find(varfracn(combi,:) > 0.9, 1, 'first');
    var90(:,combi) = [i; varn(combi,i)];
    if combi == doriexample, c = varcol;
    else, c = varcol*0.5 +0.5; end
    plot(1:N, varn(combi,:), '-', 'Color', c);
end
% add 90% variance dots
plot(var90(1,(1:noricomb) ~= doriexample), var90(2,(1:noricomb) ~= doriexample), ...
    'o', 'MarkerSize', dot90size, ...
    'MarkerFaceColor', dot90col*0.5+0.5, 'MarkerEdgeColor', 'none');
plot(var90(1,doriexample), var90(2,doriexample), 'o', 'MarkerSize', dot90size, ...
    'MarkerFaceColor', dot90col, 'MarkerEdgeColor', 'none');
xlim([1 N]);  ylim([1e-4 1]);
set(gca,'Box','off','YScale','log','XTick',[1 100 200 300]);
xlabel('principal dimension');
ylabel('variance');


%% plot example alignment of f' to principal axes
if singlefig
    subplotcm([5 4.8 3 3]);  hold on;
    text(-1,3.5,'c','Units','centimeters','FontWeight','bold',...
         'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
align90 = NaN(2, noricomb);
for combi = 1:noricomb
    i = find(alignn(combi,:) >= 0.9, 1, 'first');
    align90(:,combi) = [i alignn(combi,i)];
    if combi == doriexample, c = aligncol;
    else, c = aligncol*0.5 + 0.5; end
    plot(1:N, alignn(combi,:), '-', 'Color', c);
end
% add 90% alignment dots
plot(align90(1,(1:noricomb) ~= doriexample), align90(2,(1:noricomb) ~= doriexample), ...
    'o', 'MarkerSize', dot90size, ...
    'MarkerFaceColor', dot90col*0.5+0.5, 'MarkerEdgeColor', 'none');
plot(align90(1,doriexample), align90(2,doriexample), 'o', 'MarkerSize', dot90size, ...
    'MarkerFaceColor', dot90col, 'MarkerEdgeColor', 'none');

xlim([1 N]);  ylim([0 1]);
set(gca,'Box','off','XTick',[1 100 200 300],'YTick',[0 1]);
xlabel('principal dimension');
ylabel('cos^2(\alpha_n)'); 


%% plot example fraction total information / variance over principal axes
if singlefig
    subplotcm([1 0.8 3 3]);  hold on;
    text(-1,3.5,'d','Units','centimeters','FontWeight','bold',...
         'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
inf90 = NaN(2, noricomb);
for combi = 1:noricomb
    i = find(Ifracn(combi,:) >= 0.9, 1, 'first');
    inf90(:,combi) = [i; Ifracn(combi,i)];
    if combi == doriexample, c = varcol;
    else, c = varcol*0.5 +0.5; end
    plot(1:N, varfracn(combi,:), '-', 'Color', c);
    if combi == doriexample, c = infcol;
    else, c = infcol*0.5 + 0.5; end
    plot(1:N, Ifracn(combi,:), '-', 'Color', c);
end
% add 90% information dots
plot(inf90(1,(1:noricomb) ~= doriexample), inf90(2,(1:noricomb) ~= doriexample), ...
    'o', 'MarkerSize', dot90size, ...
    'MarkerFaceColor', dot90col*0.5+0.5, 'MarkerEdgeColor', 'none');
plot(inf90(1,doriexample), inf90(2,doriexample), 'o', 'MarkerSize', dot90size, ...
    'MarkerFaceColor', dot90col, 'MarkerEdgeColor', 'none');

xlim([1 N]);  ylim([0 1]);
set(gca,'Box','off','XTick',[1 100 200 300],'YTick',[0 1]);
xlabel('principal dimension');
ylabel('fraction variance / information');


%% per-dataset statistics
var90 = NaN(length(datasets), noricomb);
align90 = NaN(length(datasets), noricomb);
inf90 = NaN(length(datasets), noricomb);
inffrac = NaN(length(datasets), noricomb);
for i = 1:length(datasets)
    % load data
    [sr, ori, con] = loaddata(datasets{i}, 'none');
    cons = unique(con);
    assert(length(cons) == 1);  %% only deal with single-constrat data
    N = size(sr, 2);
    s = dataInfo(datasets{i});
    oricomb = s.oricomb(1:2, s.oricomb(3,:) == s.doris(doriid));
    assert(noricomb == size(oricomb,2));
    % iterate across ori combinations
    fprintf('Cross-validated statistics');
    for combi = 1:noricomb
        fprintf(' %d', combi);
        ori1 = s.oris(oricomb(1,combi));
        ori2 = s.oris(oricomb(2,combi));
        [varni, alignni, Ini] = cvstats(sr, ori, con, cons(1), ori1, ori2, cvn);
        varfracn = mean(...
            bsxfun(@rdivide, cumsum(varni, 2), sum(varni, 2)), 1);
        alignn = mean(alignni, 1);
        Ifracn = mean(bsxfun(@rdivide, Ini, Ini(:,end)), 1);
        j = find(varfracn > 0.9, 1, 'first');
        var90(i,combi) = j / N;
        inffrac(i,combi) = Ifracn(j);  % fraction of information at 90% var
        j = find(alignn >= 0.9, 1, 'first');
        align90(i,combi) = j / N;
        j = find(Ifracn >= 0.9, 1, 'first');
        inf90(i,combi) = j / N;
    end
    fprintf('\n');
end
% some stats
fprintf('\nAverage fraction of principal directions required:\n');
fprintf('90%% variance    : %5.3f +/- %5.3f (mean +/- 1SD)\n', ...
    mean(var90(:)), sqrt(var(var90(:))));
fprintf('90%% f'' alignment: %5.3f +/- %5.3f\n', ...
    mean(align90(:)), sqrt(var(align90(:))));
fprintf('90%% information : %5.3f +/- %5.3f\n', ...
    mean(inf90(:)), sqrt(var(inf90(:))));
fprintf('inf at 90%% vari : %5.3f +/- %5.3f\n\n', ...
    mean(inffrac(:)), sqrt(var(inffrac(:))));
x = reshape(inf90(:,dorisnonoverlap)-var90(:,dorisnonoverlap),[],1);
[~,p,~,stats] = ttest(x);
fprintf('Paired t-test on inf90 - var90\n');
fprintf('<inf90>-<var90> = %f +/- %f (mean +/- 1SEM)\n', ...
    mean(x), sqrt(var(x) / length(x)));
fprintf('t(%d) = %f,  p = %f\n', stats.df, stats.tstat, p); 


%% plot cumulative distributions
if singlefig
    subplotcm([5 0.8 3 3]);  hold on;
    text(-1,3.5,'e','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
alpha(histogram(reshape(var90,1,[]),'Normalization','probability',...
    'FaceColor',varcol,'EdgeColor','none','BinWidth',0.02),0.5);
alpha(histogram(reshape(align90,1,[]),'Normalization','probability',...
    'FaceColor',aligncol,'EdgeColor','none','BinWidth',0.02),0.5);
alpha(histogram(reshape(inf90,1,[]),'Normalization','probability',...
    'FaceColor',infcol,'EdgeColor','none','BinWidth',0.02),0.5);
[f,x] = ecdf(reshape(var90,1,[]));
stairs([0 x' 1],[0 f' 1],'Color',varcol);
[f,x] = ecdf(reshape(align90,1,[]));
stairs([0 x' 1],[0 f' 1],'Color',aligncol);
[f,x] = ecdf(reshape(inf90,1,[]));
stairs([0 x' 1],[0 f' 1],'Color',infcol);
xlim([0 1]);  ylim([0 1]);
set(gca,'Box','off','XTick',[0 1],'YTick',[0 1]);
xlabel('fraction principal dimensions');
ylabel('cumulative count');


%% save figure
if singlefig
    fprintf('\nWriting figure to fig8.pdf\n');
    print(['figs' filesep 'fig8'], '-dpdf');
end



function [S, mu, T] = datastats(sr, ori, con, con1, ori1, ori2)
%% returns activity statistics for selected contrast, and orientation combo

sr1 = sr(ori == ori1 & con == con1, :);
sr2 = sr(ori == ori2 & con == con1, :);
T = min(size(sr1, 1), size(sr2, 1));
sr1 = sr1(1:T, :);
sr2 = sr2(1:T, :);
ds = abs(ori1 - ori2);
ds = min(ds,abs(ds-360)) * pi / 180; % angular distance
mu = (mean(sr1, 1) - mean(sr2, 1)) / ds;
S = 0.5 * cov(sr1) + 0.5 * cov(sr2);


function [varn, alignn, In] = cvstats(sr, ori, con, con1, ori1, ori2, cvn)
%% performs cross-validated variance, f' alignment, and info scaling
%
% sr, ori, con provide dataset stats. ori1 and ori2 specify the
% orientations to discriminate. cvn determine the number of
% cross-validations. The returns variance per dimension, f' alignment, and
% information scaling, across the different cross-validations.

ds = abs(ori1 - ori2);
ds = min(ds,abs(ds-360)) * pi / 180; % angular distance
trials1 = find(ori == ori1 & con == con1);
trials2 = find(ori == ori2 & con == con1);
N = size(sr, 2);
T = min(length(trials1), length(trials2));
Tcv = floor(T / 2);
varn = NaN(cvn, N);
alignn = NaN(cvn, N);
In = NaN(cvn, N);
parfor cvi = 1:cvn
    % randomly pick Tcv train and test trials
    train1 = randperm(length(trials1), Tcv);  % indicies into trials1
    train2 = randperm(length(trials2), Tcv);  % indicies into trials2
    test1 = 1:length(trials1);
    test2 = 1:length(trials2);
    test1(train1) = [];  % testx now all non-test trials in trials1
    test2(train2) = [];
    test1 = trials1(test1(randperm(length(test1), Tcv)));
    test2 = trials2(test2(randperm(length(test2), Tcv)));
    train1 = trials1(train1);
    train2 = trials2(train2);
    % train set statistics
    sr1 = sr(train1,:); %#ok<PFBNS>
    sr2 = sr(train2,:);
    Strain = 0.5*(cov(sr1) + cov(sr2));
    [Qtrain,Dtrain] = eig(Strain);
    Dtrain = diag(Dtrain);
    [~,j] = sort(Dtrain,'descend');
    Qtrain = Qtrain(:,j);
    % test set statistics
    sr1 = sr(test1,:);
    sr2 = sr(test2,:);
    mutest = (mean(sr1,1) - mean(sr2,1)) / ds;
    Stest = 0.5*(cov(sr1) + cov(sr2));
    
    % variance in Stest for different Qtrain directions
    % we do that by computing Qtrain(:,n)' * Stest * Qtrain(:,n) for each n
    varn(cvi,:) = sum((Stest * Qtrain) .* Qtrain, 1);
    
    % cumulative alignment of ftest' to Qtrain
    mutest2 = sum(mutest.^2);
    cos2theta = (mutest * Qtrain).^2 ./ mutest2;
    alignn(cvi,:) = cumsum(cos2theta);
    
    % information scaling in test set for different Qtrain directions
    for n = 1:N
        Qnmutest = Qtrain(:,1:n)' * mutest';
        QnStest = Qtrain(:,1:n)' * Stest *  Qtrain(:,1:n);
        % when inverting QnStest, lower-bound eigenvalues to avoid
        % numerical instabilities
        [Qn,Dn] = eig(QnStest);
        Dninv = 1./diag(Dn);
        Dninv(diag(Dn) < 1e-15) = 0;
        In(cvi, n) = Qnmutest' * bsxfun(@times, Qn, Dninv') * Qn' * Qnmutest; 
    end
end
