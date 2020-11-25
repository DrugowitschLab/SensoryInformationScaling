function simPopActivity(simid)
%% simulates population activity and writes data to simData
%
% simid is the integer simulation id. Data is written to
% simData/sim[simid].mat
%
% The script supports the following simid's:
% 1 - Gaussian population with limited information
% 2 - Gaussian population with unlimited information
% 3 - LNP population with limited information
% 4 - LNP population with unlimited information

%% settings
% parameters for different simid's
Gausspar = struct(...
    'N',     [1000 1000], ...
    'Iinf',  [20 Inf], ...
    'g',     [20 20], ...
    'sigb2', [1 1], ...
    'sig02', [0.001 0.001], ...
    'Npow',  [0.1 0.1], ...
    'ds',    [15 15] * pi / 180, ...
    'T',     [1000 1000]);
P = 32;
LNPparlim = struct(...
    'N', 2500, ...    % number of neurons
    'P', P, ...       % pixels in image (one side)
    'sig', P/5, ...   % gaussian envelope sd
    'lam', P/1.5, ... % pref. wavelength
    'phi', 0, ...     % pref. spatial phase
    'c', 1, ...       % Michelson contrast
    'g', 20, ...      % tuning amplitude
    'sig_a', sqrt(2), ... % tuning ampl. variability
    'sig0', 0.25, ... % input noise (prop. of pixel range)
    'ds', 15 * pi / 180, ...
    'T', 2500);
LNPparunlim = struct(...
    'N', 2500, ...    % number of neurons
    'P', P, ...       % pixels in image (one side)
    'sig', P/5, ...   % gaussian envelope sd
    'lam', P/1.5, ... % pref. wavelength
    'phi', 0, ...     % pref. spatial phase
    'c', 1, ...       % Michelson contrast
    'g', 20, ...      % tuning amplitude
    'sig_a', sqrt(2), ... % tuning ampl. variability
    'sig0', 0, ...    % input noise (prop. of pixel range)
    'ds', 15 * pi / 180, ...
    'T', 2500);
%popsamples = 10;
popsamples = 1;
outfile = sprintf('simData%ssim%d.mat', filesep, simid);


switch simid
    case {1,2}
        %% Gaussian population model
        plotsubpop = 500;
        N = Gausspar.N(simid);
        Iinf = Gausspar.Iinf(simid);
        g = Gausspar.g(simid);
        sigb2 = Gausspar.sigb2(simid);
        sig02 = Gausspar.sig02(simid);
        Npow = Gausspar.Npow(simid);
        ds = Gausspar.ds(simid);
        T = Gausspar.T(simid);

        if exist(outfile, 'file') == 2
            fprintf('Found file %s; Loading...\n', outfile);
            d = load(outfile);
            N = d.simp.N;
            T = d.simp.T;
            Iinf = d.simp.Iinf;
            ds = d.simp.ds;
            fp = d.popmom.fp;
            Sig = d.popmom.Sig;
            Sig0 = d.popmom.Sig0;
            X1 = d.deResp(1:T,:);
            X2 = d.deResp((T+1):end,:);
        else
            fprintf('Generating population activity samples...\n');
            [fp, Sig, Sig0] = popMoments(g, sig02, sigb2, Iinf, N, Npow);
            % assume f1 and f2 to be centered around zero
            f1 = -ds * fp/2;
            f2 = -f1;
            % generate population samples X1, X2
            X1 = mvnrnd(f1, Sig, T);
            X2 = mvnrnd(f2, Sig, T);

            %% write data to file
            outfile = sprintf('simData%ssim%d.mat', filesep, simid);
            fprintf('Writing data to %s...\n', outfile);
            visOri = cat(1, zeros(T,1), (ds*180/pi)*ones(T,1));
            deResp = cat(1, X1, X2);
            visCon = ones(2*T,1);
            modelType = 'Gauss';
            par = struct('N', N, 'T', T, ...
                'Iinf', Iinf, 'g', g, 'sigb2', sigb2, 'sig02', sig02, 'ds', ds);
            popmom = struct('fp', fp, 'Sig', Sig, 'Sig0', Sig0);
            save(outfile, 'visOri', 'deResp', 'visCon', ...
                'par', 'popmom', 'modelType');
        end
        
        %% estimate info scaling
        fprintf('Estimating empirical info scaling...\n');
        S = 0.5*(cov(X1) + cov(X2));
        Ssub = 0.5*(cov(X1(:,1:plotsubpop)) + cov(X2(:,1:plotsubpop)));
        mu = (mean(X2)-mean(X1))./ds;
        Itrue = empInfscaling(fp, Sig, 1:N);
        Iemp = NaN(popsamples, N);
        for i = 1:popsamples
            if mod(i, 10) == 0, fprintf('%d ', i); end
            Iemp(i,:) = empInfscaling(mu, S, randperm(N,N), T, ds);
        end
        fprintf('\n');

        %% plot covariance spectra, info scaling
        Seig = sort(eig(S), 'descend');
        Ssubeig = sort(eig(Ssub), 'descend');
        Sigeig = sort(eig(Sig), 'descend');
        Sig0eig = sort(eig(Sig0), 'descend');
        figure('Color', 'white');
        subplot(2,1,1);  hold on;
        plot(1:N, Sig0eig, 'b-', 'LineWidth', 1);
        plot(1:N, Sigeig, 'r-', 'LineWidth', 1);
        plot(1:N, Seig, 'k-', 'LineWidth', 1);
        plot(1:plotsubpop, sort(Ssubeig, 'descend'), 'k--', 'LineWidth', 1);
        ylabel('Sig spectrum');
        legend('I0', 'I', 'I data', 'I subdata');
        subplot(2,1,2);
        semilogy(1:N, Sig0eig, 'b-', 'LineWidth', 1);  hold on;
        semilogy(1:N, Sigeig, 'r-', 'LineWidth', 1);
        semilogy(1:N, Seig, 'k-', 'LineWidth', 1);
        semilogy(1:plotsubpop, Ssubeig, 'k--', 'LineWidth', 1);
        xlabel('N');
        ylabel('Sig spectrum (log scale)');

        % info scaling
        Iincr = diff(cat(2, zeros(popsamples, 1), Iemp), 1, 2);
        Iemp_mu = mean(Iemp, 1);
        Iemp_sd = sqrt(var(Iemp, [], 1));
        Iincr_mu = mean(Iincr, 1);
        Iincr_sd = sqrt(var(Iincr, [], 1));
        figure('Color', 'white');
        subplot(2,1,1);  hold on;
        patch([1:N fliplr(1:N)], ...
            [(Iemp_mu(1,:)+Iemp_sd(1,:)) fliplr(Iemp_mu(1,:)-Iemp_sd(1,:))],1,...
            'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
        p1 = plot(1:N, Iemp_mu, 'k-', 'LineWidth', 2);
        p2 = plot(1:N, Itrue, 'r-', 'LineWidth', 2);
        p3 = plot([1 N], Itrue(end)*[1 1], 'r--');
        p4 = plot([1 N], Iinf*[1 1], 'b--');
        legend([p1 p2 p3 p4], 'from data', 'true', sprintf("I%d", N), 'Iinf');
        ylabel('Fisher info');
        subplot(2,1,2);  hold on;
        patch([1:N fliplr(1:N)], ...
            [(Iincr_mu(1,:)+Iincr_sd(1,:)) fliplr(Iincr_mu(1,:)-Iincr_sd(1,:))],1,...
            'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
        plot(1:N, Iincr_mu, 'k-', 'LineWidth', 2);
        plot(1:N, diff([0 Itrue]), 'r-', 'LineWidth', 2);
        xlabel('N');
        ylabel('Fisher info increase');
      
    case {3,4}
        %% LNP population model        
        %% load from file if already exists, otherwise generate
        if exist(outfile, 'file') == 2
            fprintf('Found file %s; Loading...\n', outfile);
            d = load(outfile);
            par = d.par;
            fp = d.popmom.fp;
            Sig = d.popmom.Sig;
            X1 = d.deResp(1:par.T,:);
            X2 = d.deResp((par.T+1):(2*par.T),:);
        else
            if simid == 3, par = LNPparlim;  % bounded or unbounded
            else, par = LNPparunlim; end
        
            %% generate neural population and compute population activity moments
            par = LNPorientedPopulation(par);
            % shuffle filter/neural order (for population subsampling)
            norder = randperm(par.N);
            par.an = par.an(norder);
            par.F = par.F(:,norder);
            % generate population activity
            X1 = LNPpopulationActivity(-par.ds/2, par.T, par);
            X2 = LNPpopulationActivity(par.ds/2, par.T, par);
            fp = (mean(X2) - mean(X1)) / par.ds;
            Sig = 0.5*(cov(X1) + cov(X2));
            % concatenate into trial sequence
            visOri = cat(1, zeros(par.T,1), (par.ds*180/pi)*ones(par.T,1));
            deResp = cat(1, X1, X2);
            visCon = ones(2*par.T,1);
            modelType = 'LNP';
            popmom = struct('fp', fp, 'Sig', Sig);
            save(outfile, 'visOri', 'deResp', 'visCon', ...
                'par', 'popmom', 'modelType');
        end
        
        fprintf('Estimating empirical info scaling...\n');
        N = par.N;
        T = par.T;
        ds = par.ds;
        Iemp = NaN(popsamples, N);
        for i = 1:popsamples
            if mod(i, 10) == 0, fprintf('%d ', i); end
            Iemp(i,:) = empInfscaling(fp, Sig, randperm(N), T, ds);
        end
        fprintf('\n');
        
        figure('Color', 'white');
        Iemp_mu = mean(Iemp, 1);
        Iemp_sd = sqrt(var(Iemp, [], 1));
        figure('Color', 'white');
        subplot(2,1,1);  hold on;
        patch([1:N fliplr(1:N)], ...
            [(Iemp_mu(1,:)+Iemp_sd(1,:)) fliplr(Iemp_mu(1,:)-Iemp_sd(1,:))],1,...
            'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
        p1 = plot(1:N, Iemp_mu, 'k-', 'LineWidth', 2);
        p2 = plot([1 N], Iemp_mu(end)*[1 1], 'r--');
        legend([p1 p2], 'from data', sprintf("I(%d)", N));
        ylabel('Fisher info');
        
    otherwise
        error('Unknown simid');
end





function X = LNPpopulationActivity(theta, T, par)
%% returns T trials of generated activity for stimulus theta

P = sqrt(size(par.F, 1));
Jtheta = pixelizedGarbor(P, theta, par);
JthetaNorm = standardizevector(reshape(Jtheta, P^2, []));
% neural activity (non-thresholded) from T noisy images
X = bsxfun(@plus, JthetaNorm', par.sig0 * randn(T, P^2)) * par.F;
X(X < 0) = 0;    % threshold-linear activation function
X = poissrnd(X); % neural activity


function par = LNPorientedPopulation(par)
%% generates a population of neurons with parameters par

N = par.N; 
P = par.P;    
% this alpha yields a logN draw with mean 0, variance sig_a^2
alpha = par.sig_a^2 + 1;
par.an = par.g * lognrnd(-log(sqrt(alpha)), sqrt(log(alpha)), N, 1);
ThetaPref = linspace(-pi, pi* (N - 1)/ N, N);
J = NaN(P, P, N);
for ii = 1: N
    J(:, :, ii) = pixelizedGarbor(P, ThetaPref(ii), par);
end
par.F = standardizevector(reshape(J, P^2, [])); % compute neural filters


function J = pixelizedGarbor(P, theta, par)
%% gernates P x P pixelated Garbor of orientation theta with parameters par

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

function Y = standardizevector(X)
%% standardize vector X or columns of matrix X
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
