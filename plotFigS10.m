function plotFigS10
% plots how information scaling changes with dtheta, using rot. pop. code


%% model parameters
singlefig = true;  % format for paper if true
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
    'T', 2500);
% simulation parameters
N = par.N;
theta_sim = 0;
dss = [1 15 30 45 60 75 90];  % in deg
dsN = length(dss);
popsamples = 50;
dscol = [241 238 246; 208 209 230; 166 189 219; 116 169 207; ...
    54 144 192; 5 112 176; 3 78 123] / 255;
assert(dsN <= size(dscol, 1));
plotIlims = [0 130];


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


%% generate population
par = OrientedPopulation(par);


%% compute info scaling for different dss
Imu = NaN(dsN, N);
Isd = NaN(dsN, N);
In = NaN(popsamples, N);
for dsi = 1:dsN
    fprintf('%d ', dss(dsi));
    % get population response moments
    ds = dss(dsi) * pi / 180;  % ds in rad
    f1 = popresp(par, theta_sim - ds/2);
    f2 = popresp(par, theta_sim + ds/2);
    Sigma = 0.5 * (poprespcov(par, theta_sim - ds/2) + poprespcov(par, theta_sim + ds/2));  
    muq = (f2 - f1) / ds;
    f = muq(1:N)';
    Sig = Sigma(1:N,1:N);
    % estimate information based on those
    parfor popi = 1:popsamples
        In(popi, :) = empInfscaling(f,Sig,randperm(N));
    end
    Imu(dsi,:) = mean(In,1);
    Isd(dsi,:) = sqrt(var(In,[],1));
end
fprintf('\n');


%% plot info scaling
if singlefig
    figure('Name', 'Figure', 'Units', 'centimeters', 'Position', [0 0 8.5 4.5]);    
    subplotcm([1 0.8 4 3]);  hold on;
    text(-1,3.5,'b','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
% shaded error bars
for dsi = 1:dsN
    patch([1:N fliplr(1:N)], ...
      [(Imu(dsi,:)+Isd(dsi,:)) fliplr(Imu(dsi,:)-Isd(dsi,:))],1,...
      'FaceColor', 0.7+0.3*dscol(dsi,:), 'EdgeColor', 'none');
end
% actual scaling + final point
for dsi = 1:dsN
    plot(1:N, Imu(dsi,:), '-', 'Color', dscol(dsi,:));
end
for dsi = 1:dsN
    plot([N N], Imu(dsi,N)+Isd(dsi,N)*[-1 1], '-', 'Color', dscol(dsi,:));
    plot(N, Imu(dsi,N), 'o', 'MarkerSize', 3, 'MarkerFaceColor', dscol(dsi,:), ...
        'MarkerEdgeColor', 'none');
end
ylim(plotIlims);
xlabel('N');
ylabel('Fisher information');
set(gca,'Box','off','YTick',[0 100],'XTick',[1 500 1000]);

% show final info measure
if singlefig
    subplotcm([5.5 0.8 2.5 3]);  hold on;
    text(-1,3.5,'c','Units','centimeters','FontWeight','bold',...
        'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
else
    figure;  hold on;
end
plot(dss, Imu(:,N), '-', 'Color', [1 1 1]*0.1, 'LineWidth', 1);
for dsi = 1:dsN
    plot([1 1]*dss(dsi), Imu(dsi,N)+Isd(dsi,N)*[-1 1], '-', 'Color', dscol(dsi,:));
    plot(dss(dsi), Imu(dsi,N), 'o', 'MarkerSize', 3, 'MarkerFaceColor', dscol(dsi,:), ...
        'MarkerEdgeColor', 'none');
end
ylim(plotIlims);
xlabel('d\theta');
set(gca,'XTick',dss,'YColor','none');


%% write figure to file
if singlefig
    fprintf('\nWriting figure to figS10.pdf\n');
    print(['figs' filesep 'figS10'], '-dpdf');
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
    J(:, :, ii) = pixelizedGarbor(P, ThetaPref(ii), par);
end
par.F = standardizevector(reshape(J, P^2, [])); % compute neural filters


function J = pixelizedGarbor(P, theta, par)
% gernates P x P pixelated Garbor of orientation theta with parameters par

sig2 = -0.5/ par.sig^2;    
cosTheta = 2*pi *cos(theta)/ par.lam;    
sinTheta = 2* pi *sin(theta)/ par.lam;    
x = linspace(-P/2 + 0.5, P/2 - 0.5, P); 
y = linspace(-P/2 + 0.5, P/2 - 0.5, P);
J = NaN(P, P);
for jj = 1: P  
    J(:, jj) = par.c * ...
        exp(sig2 * (x.^2 + y(jj)^2)) .* ...
        cos(x * cosTheta + y(jj) * sinTheta + par.phi);    
end


function ftheta = popresp(par, theta)
% computes mean population response of population par to stim theta

P = sqrt(size(par.F, 1));
Jtheta = pixelizedGarbor(P, theta, par);
JthetaNorm = standardizevector(reshape(Jtheta, P^2, []));
ftheta = par.an.*(par.F'* JthetaNorm);
ftheta(ftheta < 0) = 0;


function SigTheta = poprespcov(par, theta)
% computes population response covariance to stim theta

P = sqrt(size(par.F, 1));    
Jtheta = pixelizedGarbor(P, theta, par);
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
