function plotFigS15
%% plots figure showing impact of saturating dF/F on simulated data


%% general settings
singlefig = true;   % set to true to have all planels in single figure
datafile = ['.' filesep 'data' filesep 'm25_170512.mat'];
ecdfsamples = 1000;

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
    figure('Name', 'Figure', 'Units', 'centimeters', 'Position', [0 0 12 12.8]);
end


%% model parameters
% network parameters
P = 32;
par = struct(...
    'N', 300, ...    % number of neurons
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
theta_sim = 0;
ds = 5 * pi / 180;  % in rad

% plot colors and saturation function
plotIlims = [0 150];
oricol = [1 1 1]*0.2;
lincol = [77 175 73] / 255;
quadcol = [55 126 184] / 255;
sat1col = [228 26 28] / 255;
sat2col = [236 102 103] / 255;
sat3col = [245 178 180] / 255;
sat1fn = @(sc) 0.5 * tanh((0.05*sc).^2).^0.5;
sat2fn = @(sc) 0.24 * tanh((0.1*sc).^2).^0.5;
sat3fn = @(sc) 0.05 * tanh((0.45*sc).^2).^0.5;


%% generate population activity
par = OrientedPopulation(par);
fprintf('Simulating %d trials...\n', par.T);
[r1, r2] = simulateTrials(par, theta_sim, ds, 0.5);

% load datafile and fit mapping
fprintf('Loading %s ...\n', datafile); 
d = load(datafile);
pcs = 1:99;
r1prctls = prctile(r1(:),pcs);
dFFprctls = prctile(d.deResp(:),pcs);
linm = fitlm(r1prctls, dFFprctls);
quadm = fitlm(r1prctls, dFFprctls, 'quadratic');

% plot relationship between spike count and dF/F
if singlefig
    subplotcm([1 8.3 4.5 3.5]);  hold on;
    text(-1,3.5,'a','Units','centimeters','FontWeight','bold',...
         'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure('Color', 'white');  hold on;
end
plot(r1prctls, dFFprctls, 'Color', oricol, 'LineWidth', 1);
r1s = linspace(min(r1prctls), max(r1prctls), 100);
plot(r1s, predict(linm, r1s'), 'Color', lincol, 'LineWidth', 0.5);
plot(r1s, predict(quadm, r1s'), 'Color', quadcol, 'LineWidth', 0.5);
plot(r1s, sat1fn(r1s), 'Color', sat1col, 'LineWidth', 0.5);
plot(r1s, sat2fn(r1s), 'Color', sat2col, 'LineWidth', 0.5);
plot(r1s, sat3fn(r1s), 'Color', sat2col, 'LineWidth', 0.5);
xlim([0 25]);  ylim([-0.05 0.8]);
xlabel('model spike count');
ylabel('dF/F');
set(gca,'Box','off','XTick',[0 10 20],'YTick',[0 0.2 0.4 0.6 0.8]);


% plot dF/F CDF for different models
if singlefig
    subplotcm([6.5 8.3 4.5 3.5]);  hold on;
    text(-1,3.5,'b','Units','centimeters','FontWeight','bold',...
         'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure('Color', 'white');  hold on;
end
[f,x] = ecdf(d.deResp(round(linspace(1,numel(d.deResp),ecdfsamples))));
stairs([0 x'], [0 f'], 'Color', oricol, 'LineWidth', 1);
[f,x] = ecdf(r1(round(linspace(1,numel(r1),ecdfsamples))));
stairs(predict(linm, [0; x]), [0 f'], 'Color', lincol, 'LineWidth', 0.5);
stairs(predict(quadm, [0; x]), [0 f'], 'Color', quadcol, 'LineWidth', 0.5);
stairs(sat1fn([0 x']), [0 f'], 'Color', sat1col, 'LineWidth', 0.5);
stairs(sat2fn([0 x']), [0 f'], 'Color', sat2col, 'LineWidth', 0.5);
stairs(sat3fn([0 x']), [0 f'], 'Color', sat3col, 'LineWidth', 0.5);
xlim([-0.05 0.5]);  ylim([0 1]);
xlabel('dF/F');  ylabel('cumulative fraction');
set(gca,'Box','off','XTick',[0 0.2 0.4],'YTick',[0 1]);


% compute info scaling for different types of saturation
if singlefig
    subplotcm([1 4.3 3 3]);  hold on;
    text(-1,3.5,'c','Units','centimeters','FontWeight','bold',...
         'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure('Color', 'white');  hold on;
end
fprintf('Information scaling for spike count data ...\n');
[Imuori, Isdori] = avgInfoScaling(r1, r2, ds);
[Imuorishuf, ~] = avgInfoScaling(r1, r2, ds, true);
alpha(patch([0:N fliplr(0:N)], [(Imuori+Isdori) fliplr(Imuori-Isdori)], 1, ...
        'FaceColor', oricol, 'EdgeColor', 'none'), 0.2);
plot(0:par.N, Imuori, '-', 'Color', oricol, 'LineWidth', 1);
plot(0:par.N, Imuorishuf, '--', 'Color', oricol, 'LineWidth', 1);
xlim([0 N]);  ylim(plotIlims);
title('spike counts');
set(gca,'Box','off','XColor','none','YTick',[0 50 100 150]);


if singlefig
    subplotcm([4.5 4.3 3 3]);  hold on;
else
    figure('Color', 'white');  hold on;
end
fprintf('Information scaling for linear model ...\n');
r1m = reshape(predict(linm, reshape(r1, [], 1)), par.T, par.N);
r2m = reshape(predict(linm, reshape(r2, [], 1)), par.T, par.N);
[Imu, Isd] = avgInfoScaling(r1m, r2m, ds);
[Imushuf, ~] = avgInfoScaling(r1m, r2m, ds, true);
plot(0:par.N, Imuori, '-', 'Color', 0.5*oricol+0.5, 'LineWidth', 0.5);
plot(0:par.N, Imuorishuf, '--', 'Color', 0.5*oricol+0.5, 'LineWidth', 0.5);
alpha(patch([0:N fliplr(0:N)], [(Imu+Isd) fliplr(Imu-Isd)], 1, ...
        'FaceColor', lincol, 'EdgeColor', 'none'), 0.2);
plot(0:par.N, Imu, '-', 'Color', lincol, 'LineWidth', 1);
plot(0:par.N, Imushuf, '--', 'Color', lincol, 'LineWidth', 1);
xlim([0 N]);  ylim(plotIlims);
title('linear model');
set(gca,'Box','off','XColor','none','YColor','none');


if singlefig
    subplotcm([8 4.3 3 3]);  hold on;
else
    figure('Color', 'white');  hold on;
end
fprintf('Information scaling for quadratic model ...\n');
r1m = reshape(predict(quadm, reshape(r1, [], 1)), par.T, par.N);
r2m = reshape(predict(quadm, reshape(r2, [], 1)), par.T, par.N);
[Imu, Isd] = avgInfoScaling(r1m, r2m, ds);
[Imushuf, ~] = avgInfoScaling(r1m, r2m, ds, true);
plot(0:par.N, Imuori, '-', 'Color', 0.5*oricol+0.5, 'LineWidth', 0.5);
plot(0:par.N, Imuorishuf, '--', 'Color', 0.5*oricol+0.5, 'LineWidth', 0.5);
alpha(patch([0:N fliplr(0:N)], [(Imu+Isd) fliplr(Imu-Isd)], 1, ...
        'FaceColor', quadcol, 'EdgeColor', 'none'), 0.2);
plot(0:par.N, Imu, '-', 'Color', quadcol, 'LineWidth', 1);
plot(0:par.N, Imushuf, '--', 'Color', quadcol, 'LineWidth', 1);
xlim([0 N]);  ylim(plotIlims);
title('quadratic model');
set(gca,'Box','off','XColor','none','YColor','none');


if singlefig
    subplotcm([1 0.8 3 3]);  hold on;
else
    figure('Color', 'white');  hold on;
end
fprintf('Information scaling for saturating model 1 ...\n');
r1m = sat1fn(r1);
r2m = sat1fn(r2);
[Imu, Isd] = avgInfoScaling(r1m, r2m, ds);
[Imushuf, ~] = avgInfoScaling(r1m, r2m, ds, true);
plot(0:par.N, Imuori, '-', 'Color', 0.5*oricol+0.5, 'LineWidth', 0.5);
plot(0:par.N, Imuorishuf, '--', 'Color', 0.5*oricol+0.5, 'LineWidth', 0.5);
alpha(patch([0:N fliplr(0:N)], [(Imu+Isd) fliplr(Imu-Isd)], 1, ...
        'FaceColor', sat1col, 'EdgeColor', 'none'), 0.2);
plot(0:par.N, Imu, '-', 'Color', sat1col, 'LineWidth', 1);
plot(0:par.N, Imushuf, '--', 'Color', sat1col, 'LineWidth', 1);
xlim([0 N]);  ylim(plotIlims);
title('saturating model 1');
set(gca,'Box','off','XTick',[0 100 200 300],'YTick',[0 50 100 150]);


if singlefig
    subplotcm([4.5 0.8 3 3]);  hold on;
else
    figure('Color', 'white');  hold on;
end
fprintf('Information scaling for saturating model 2 ...\n');
r1m = sat2fn(r1);
r2m = sat2fn(r2);
[Imu, Isd] = avgInfoScaling(r1m, r2m, ds);
[Imushuf, ~] = avgInfoScaling(r1m, r2m, ds, true);
plot(0:par.N, Imuori, '-', 'Color', 0.5*oricol+0.5, 'LineWidth', 0.5);
plot(0:par.N, Imuorishuf, '--', 'Color', 0.5*oricol+0.5, 'LineWidth', 0.5);
alpha(patch([0:N fliplr(0:N)], [(Imu+Isd) fliplr(Imu-Isd)], 1, ...
        'FaceColor', sat2col, 'EdgeColor', 'none'), 0.2);
plot(0:par.N, Imu, '-', 'Color', sat2col, 'LineWidth', 1);
plot(0:par.N, Imushuf, '--', 'Color', sat2col, 'LineWidth', 1);
xlim([0 N]);  ylim(plotIlims);
title('saturating model 2');
set(gca,'Box','off','XTick',[0 100 200 300],'YColor','none');


if singlefig
    subplotcm([8 0.8 3 3]);  hold on;
else
    figure('Color', 'white');  hold on;
end
fprintf('Information scaling for saturating model 3 ...\n');
r1m = sat3fn(r1);
r2m = sat3fn(r2);
[Imu, Isd] = avgInfoScaling(r1m, r2m, ds);
[Imushuf, ~] = avgInfoScaling(r1m, r2m, ds, true);
plot(0:par.N, Imuori, '-', 'Color', 0.5*oricol+0.5, 'LineWidth', 0.5);
plot(0:par.N, Imuorishuf, '--', 'Color', 0.5*oricol+0.5, 'LineWidth', 0.5);
alpha(patch([0:N fliplr(0:N)], [(Imu+Isd) fliplr(Imu-Isd)], 1, ...
        'FaceColor', sat3col, 'EdgeColor', 'none'), 0.2);
plot(0:par.N, Imu, '-', 'Color', sat3col, 'LineWidth', 1);
plot(0:par.N, Imushuf, '--', 'Color', sat3col, 'LineWidth', 1);
xlim([0 N]);  ylim(plotIlims);
title('saturating model 3');
set(gca,'Box','off','XTick',[0 100 200 300],'YColor','none');


%% save figure
if singlefig
    fprintf('\nWriting figure to figS15.pdf\n');
    print(['figs' filesep 'figS15'], '-dpdf');
end



function par = OrientedPopulation(par)
% generates a population of neurons with parameters par

N = par.N; 
P = par.P;    
alpha = par.sig_a^2 + 1;
par.an = par.g * lognrnd(-log(sqrt(alpha)), sqrt(log(alpha)), N, 1);
ThetaPref = linspace(-pi, pi* (N - 1)/ N, N);
J = NaN(P, P, N);
for ii = 1: N
    J(:, :, ii) = pixelizedGarbor(par, ThetaPref(ii), [0 0]);
end
par.F = standardizevector(reshape(J, P^2, [])); % compute neural filters


function J = pixelizedGarbor(par, theta, mu)
% gernates P x P pixelated Garbor of orientation theta with parameters par,
% centered on mu = [x y] (in pixels);

sig2 = -0.5/ par.sig^2;
P = par.P;
cosTheta = 2*pi *cos(theta)/ par.lam;    
sinTheta = 2* pi *sin(theta)/ par.lam;    
x = linspace(-P/2 + 0.5, P/2 - 0.5, P) + mu(1); 
y = linspace(-P/2 + 0.5, P/2 - 0.5, P) + mu(2);
J = NaN(P, P);
for jj = 1: P  
    J(:, jj) = par.c * ...
        exp(sig2 * (x.^2 + y(jj)^2)) .* ...
        cos(x * cosTheta + y(jj) * sinTheta + par.phi);    
end


function [r1, r2] = simulateTrials(par, theta, ds, dur)
%% simulates a set of trials of two orientations

r1 = NaN(par.T, par.N);
r2 = NaN(par.T, par.N);
theta1 = theta - ds/2;
theta2 = theta + ds/2;
P2 = par.P^2;
sig0 = par.sig0;
ReLU = @(x) (x >= 0) .* x;
for i = 1:par.T
    % response to theta1
    Jtheta = standardizevector(...
        reshape(pixelizedGarbor(par, theta1, [0 0]), P2, []));
    r1(i,:) = poissrnd(ReLU(dur * par.an .* ...
                            (par.F' * (Jtheta + sig0 * randn(P2, 1)))));
    % response to theta2
    Jtheta = standardizevector(...
        reshape(pixelizedGarbor(par, theta2, [0 0]), P2, []));
    r2(i,:) = poissrnd(ReLU(dur * par.an .* ...
                            (par.F' * (Jtheta + sig0 * randn(P2, 1)))));
end


function ftheta = popresp(par, theta, mu)
% computes mean population response of population par to stim theta

P = sqrt(size(par.F, 1));
Jtheta = pixelizedGarbor(par, theta, mu);
JthetaNorm = standardizevector(reshape(Jtheta, P^2, []));
ftheta = par.an.*(par.F'* JthetaNorm);
ftheta(ftheta < 0) = 0;


function SigTheta = poprespcov(par, theta)
% computes population response covariance to stim theta

P = sqrt(size(par.F, 1));    
Jtheta = pixelizedGarbor(par, theta, [0 0]);
JthetaNorm = standardizevector(reshape(Jtheta, P^2, []));
sig02 = par.sig0^2;
ftheta = par.an'.*(JthetaNorm'* par.F);
FF = par.F'*par.F;
ftheta(ftheta < 0) = 0;
FF(FF < 0) = 0;
SigTheta = sig02 * (par.an*par.an') .* FF + diag(ftheta);            

       
function Y = standardizevector(X)
% standardize vector X or columns of matrix X
[n, m] = size(X);
if n==1 || m ==1
    Y = X- mean(X);
    Y = Y/norm(Y);
else
    Y = NaN(size(X));
    for jj = 1: m
        Y(:, jj) = standardizevector(X(:, jj));
    end
end   


function [Imu, Isd] = avgInfoScaling(r1, r2, ds, shuf)
%% computes average info scaling across different orders
%
% performs trial-identify shuffling if [optimal] shuf is given and true
shuf = nargin >= 4 && logical(shuf);
ordn = 50;
[T, N] = size(r1);

% shuffle trial-order per neuron, if requested
if shuf
    for n = 2:N
        r1(:,n) = r1(randperm(T),n);
        r2(:,n) = r2(randperm(T),n);
    end
end

% pop statistics
fp = (mean(r1,1) - mean(r2,1)) / ds;
Sig = 0.5*(cov(r1) + cov(r2));

% info scaling across different ordering
In = NaN(ordn, N);
for ordi = 1:ordn
    if mod(ordi, 10) == 0, fprintf('%d ', ordi); end
    In(ordi,:) = empInfscaling(fp, Sig, randperm(N), T, ds);
end
fprintf('\n');

% statistics
Iincr = diff([zeros(ordn,1) In], 1, 2);
Imu = [0 mean(cumsum(Iincr, 2))];
Isd = [0 sqrt(cumsum(var(Iincr, [], 1)))];
