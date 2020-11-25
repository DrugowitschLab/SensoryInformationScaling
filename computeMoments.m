function computeMoments(dataset, ori1id, ori2id, conid, popsamples, subsample, subdims)
addpath('shared')
%% computes the fisher info increase moments for given dataset/ori/con
%
% The function can be called as
%
% computeMoments(dataset)
%
% Outputs dataset information.
%
% computeMoments(dataset, ori1id, ori2id, conid, popsamples)
% computeMoments(dataset, ori1id, ori2id, conid, popsamples, subsample)
% computeMoments(dataset, ori1id, ori2id, conid, popsamples, subsample, subdims)
%
% Computes the fisher info increase moments for the given empirical
% dataset, combination of provided orientations, and given contrast, using
% popsample bootstrap samples. Aruments are
% - dataset: name of the empirical dataset
% - ori1id, ori2id: indices of the first/second orientation
% - conid: contrast index
% - popsamples: number of bootstrap samples
% - subsample (optional, defaults to 'none'): supports loading trials
%   subsamples, and can take the same values as for loaddata(.).
% - subdims (optimal): if given, data will be dimensionality-reduced to
%   subdims dimensions, and information scaling is only computed up to
%   those dimensions.
% The moments are written to
% moment_cache/[dataset]_o[ori1id]-[ori2id]_c[conid]_[subsample]_[subdims].mat
% The 'subsample' part is only added if not 'none', or if 'subdims' is
% provided.
%
% computeMoments(dataset, doriid, conid, popsamples)
% computeMoments(dataset, doriid, conid, popsamples, subsample)
% computeMoments(dataset, doriid, conid, popsamples, subsample, subdims)
%
% Computes the fisher info increase moments for the given empirical
% dataset, providing proxy moments across all possible orientation
% combinations whose orientation difference if given by the doriid. The
% arguments are the same as before, except for
% - doriid: index of orientation difference
% The moments are written to
% moment_cache/[dataset]_dori[doriid]_c[conid]_[subsample]_[subdims].mat
% The 'subsample' part is only added if not 'none', or if 'subdims' is
% provided.
% 
% computeMoments(dataset, N, T, popsamples)
%
% Computed the fisher info increase for simulated datasets. The arguments
% are
% - dataset: simulated dataset name, with signature 'simx'
% - N: number of neurons to consider
% - T: number of trials to consider
% - popsamples: number of bootstrap samples
% The moments are written to
% moment_cache/[dataset]_N[N]_T[T].mat


%% settings
doricols = [     0    0.4470    0.7410; ...   % get(gca, 'colororder') + add 
            0.8500    0.3250    0.0980; ...
            0.9290    0.6940    0.1250; ...
            0.4940    0.1840    0.5560; ...
            0.4660    0.6740    0.1880; ...
            0.3010    0.7450    0.9330; ...
            0.6350    0.0780    0.1840; ...
            0.5000    0.5000    0.8000; ...
                 0    0.3803    0.0470];


%% only dataset - provide dataset info
if nargin == 1
    [sr, ori, con] = loaddata(dataset);
    if startsWith(dataset,'sim')
        fprintf('Simulated dataset with\n- %d neurons\n- %d trials\n', ...
            size(sr,2), size(sr,1)/2);
    else
        oris = unique(ori);
        cons = unique(con);
        fprintf('Orientations for dataset %s\n', dataset);
        for i = 1:length(oris)
            fprintf('%2d: %5.1f\n', i, oris(i));
        end
        fprintf('Contrasts for dataset %s\n', dataset);
        for i = 1:length(cons)
            fprintf('%2d: %5.3f\n', i, cons(i));
        end
    end
    return
end


%% identify how function was called
if startsWith(dataset,'sim')
    % simulated dataset
    fmode = 'sim';
    if nargin > 4
        error('Only accept 4 arguments for simulated datasets');
    end
    N = ori1id;
    T = ori2id;
    Tb = T; 
    popsamples = conid;
    subsample = 'none';
    ori1id = 1;
    ori2id = 2;
    conid = 1;
elseif nargin < 5 || ischar(popsamples)
    % empirical dataset, dori mode
    fmode = 'dori';
    if nargin > 6
        error('Only accept 6 arguments for empirical data, dori mode');
    end
    if nargin >= 6, subdims = subsample; else, subdims = NaN; end
    if nargin >= 5, subsample = popsamples; else, subsample = 'none'; end
    popsamples = conid;
    conid = ori2id;
    doriid = ori1id;
else
    % empirical dataset, oricomb mode
    fmode = 'oricomb';
    if nargin > 7
        error('Only accept 7 arguments for empirical data, oricomb mode');
    end
    if nargin < 7, subdims = NaN; end
    if nargin < 6, subsample = 'none'; end
end


%% load data and check according to mode
[sr, ori, con] = loaddata(dataset, subsample);
oris = unique(ori);
cons = unique(con);
Tb   = inf; 

if any(strcmp(subsample,{'lospdb'})) 
    [sralt, orialt, conalt] = loaddata(dataset, 'hispdb');
    sr1alt = sralt(orialt == oris(ori1id) & conalt == cons(conid), :);
    sr2alt = sralt(orialt == oris(ori2id) & conalt == cons(conid), :);
    Tb = min([size(sr1alt, 1) size(sr2alt, 1)]);
elseif any(strcmp(subsample,{'hispdb'})) 
    [sralt, orialt, conalt] = loaddata(dataset, 'lospdb');
    sr1alt = sralt(orialt == oris(ori1id) & conalt == cons(conid), :);
    sr2alt = sralt(orialt == oris(ori2id) & conalt == cons(conid), :);
    Tb = min([size(sr1alt, 1) size(sr2alt, 1)]);
end



if strcmp(fmode, 'sim')
    % sim mode
    if N > size(sr, 2)
        error('Demanded neuron# larger than available (%d > %d)', ...
            N, size(sr, 2));
    end
    T1 = sum(ori == oris(ori1id) & con == cons(conid));
    T2 = sum(ori == oris(ori2id) & con == cons(conid));
    if T > min(T1, T2)
        error('Demanded trials larger than available (%d > %d)', ...
            T, min(T1, T2));
    end
    fprintf('Using N = %d, T = %d from dataset %s\n', N, T, dataset);
    outfile = sprintf('moment_cache%s%s_N%d_T%d.mat', ...
        filesep, dataset, N, T);
    sr = sr(:,1:N);  % ensure that we only use N neurons
    subdims = N;
else
    % dori and oricomb modes
    d = dataInfo(dataset);
    N = size(sr, 2);
    T = Inf;
    if isnan(subdims), subdims = N; else, subdims = min(subdims, N); end
    if subdims == N, subdim_str = '';
    else, subdim_str = sprintf('_%d', subdims); end
    if strcmp(subsample, 'none'), subsample_str = '';
    else, subsample_str = ['_' subsample]; end
    if strcmp(fmode, 'dori')
        fprintf('Using dori=%d, con=%4.2f from dataset %s\n', ...
            d.doris(doriid), cons(conid), dataset);
        outfile = sprintf('moment_cache%s%s_dori%d_c%d%s%s.mat', filesep, ...
            dataset, doriid, conid, subsample_str, subdim_str);
    else
        fprintf('Using ori1=%d, ori2=%d, con=%4.2f from dataset %s\n', ...
            oris(ori1id), oris(ori2id), cons(conid), dataset);
        outfile = sprintf('moment_cache%s%s_o%d-%d_c%d%s%s.mat', filesep, ...
            dataset, ori1id, ori2id, conid, subsample_str, subdim_str);
    end
end
outexists = (exist(outfile, 'file') == 2);


%% process data
switch fmode
    case {'sim', 'oricomb'}
        [mu, S, ds, T, Nmax] = datamoments(sr, ori, con, ...
            oris(ori1id), oris(ori2id), cons(conid), subdims, min([T Tb]));
        
        
        % compute moments
        fprintf('Computing information estimates (N=%d, Nmax=%d) ...\n', N, Nmax);
        [Iincr_mu, Iincr_var, Iincr_samples] = ...
            estIincrMoments(mu, S, T, ds, popsamples, Nmax);
        
        % plot moments
        figure('Color', 'white');
        subplot(2, 1, 1);  hold on;
        patch([1:Nmax fliplr(1:Nmax)], ...
            [(Iincr_mu+sqrt(Iincr_var)) fliplr(Iincr_mu-sqrt(Iincr_var))],1,...
            'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
        plot(1:Nmax, Iincr_mu, 'k-', 'LineWidth', 2);
        plot([1 Nmax], [0 0], 'k--');
        ylabel('Info increase estimate');

        subplot(2, 1, 2);  hold on;
        I_var = cumsum(Iincr_var);
        I_mu = cumsum(Iincr_mu);
        patch([1:Nmax fliplr(1:Nmax)], ...
            [(I_mu+sqrt(I_var)) fliplr(I_mu-sqrt(I_var))],1,...
            'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
        plot(1:Nmax, I_mu, 'k-', 'LineWidth', 2);
        xlabel('N');
        ylabel('Info estimate');
        
    case 'dori'
        % iterate over different orientation combinations
        oricomb = d.oricomb(1:2, d.oricomb(3,:) == d.doris(doriid));
        oricombn = size(oricomb, 2);
        % find required T
        T = Inf;
        contrials = con == cons(conid);
        for i = 1:oricombn
            T = min(T, sum(ori == oris(oricomb(1,i)) & contrials));
            T = min(T, sum(ori == oris(oricomb(2,i)) & contrials));            
        end
        Nmax = min([N (2*T-3) subdims]);
        if outexists
            % load moments
            fprintf('Found moments file %s - loading...\n', outfile);
            m = load(outfile);
            assert(Nmax == length(m.Iincr_mu));
            Iincr_mu = m.Iincr_mu;
            Iincr_var = m.Iincr_var;
            Iincr_samples = m.Iincr_samples;
        else
            % compute moments
            rngstate = rng();
            fprintf('%d orientation pair(s) with dori = %d\n', oricombn, d.doris(doriid));
            Iincr_samples = NaN(popsamples, oricombn, Nmax);
            for i = 1:oricombn
                ori1 = oris(oricomb(1,i));
                ori2 = oris(oricomb(2,i));
                fprintf('Computing information estimates %d vs. %d (N=%d, Nmax=%d) ...\n', ...
                    ori1, ori2, N, Nmax);
                [mu, S, ds, T, Nmax] = datamoments(sr, ori, con, ...
                    ori1, ori2, cons(conid), subdims, T);
                rng(rngstate);  % ensure the same Norder permutation sequence
                [~, ~, Iincr_samples(:,i,:)] = ...
                    estIincrMoments(mu, S, T, ds, popsamples, Nmax);
            end
        
            % combine moments across orientation pairs
            Iincr_mu = NaN(1, Nmax);
            Iincr_var = NaN(1, Nmax);
            for n = 1:Nmax
                mu = mean(Iincr_samples(:,:,n), 1);
                S = cov(Iincr_samples(:,:,n));
                Sinv1 = S \ ones(oricombn, 1);
                Iincr_var(n) = 1 / sum(Sinv1);
                Iincr_mu(n) = Iincr_var(n) * (mu * Sinv1);
            end
        end
        
        % plot moments
        figure('Color', 'white');
        subplot(2, 1, 1);  hold on;
        patch([1:Nmax fliplr(1:Nmax)], ...
            [(Iincr_mu+sqrt(Iincr_var)) fliplr(Iincr_mu-sqrt(Iincr_var))],1,...
            'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
        plot(1:Nmax, Iincr_mu, 'k-', 'LineWidth', 2);
        for i = 1:oricombn
            plot(1:Nmax, reshape(mean(Iincr_samples(:,i,:),1),1,[]), ...
                '-', 'Color', doricols(i,:));
        end
        plot([1 Nmax], [0 0], 'k--');
        ylabel('Info increase estimate');

        subplot(2, 1, 2);  hold on;
        I_var = cumsum(Iincr_var);
        I_mu = cumsum(Iincr_mu);
        patch([1:Nmax fliplr(1:Nmax)], ...
            [(I_mu+sqrt(I_var)) fliplr(I_mu-sqrt(I_var))],1,...
            'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
        plot(1:Nmax, I_mu, 'k-', 'LineWidth', 2);
        for i = 1:oricombn
            plot(1:Nmax, cumsum(reshape(mean(Iincr_samples(:,i,:),1),1,[])), ...
                '-', 'Color', doricols(i,:));
        end
        xlabel('N');
        ylabel('Info estimate');
        
        % plot average covariance
        avgcov = zeros(oricombn,oricombn);
        for n = 1:Nmax
            avgcov = avgcov + cov(Iincr_samples(:,:,n));
        end
        avgcov = avgcov / Nmax;
        avgcorr = avgcov ./ sqrt(diag(avgcov) * diag(avgcov)');
        comblabels = cell(1, oricombn);
        for i = 1:oricombn
            comblabels{i} = sprintf('%d-%d', oris(oricomb(1,i)), oris(oricomb(2,i)));
        end
        figure('Color', 'white');
%        b = bar3(avgcorr);
        b = bar3(avgcov);
        for k = 1:length(b)
            zdata = b(k).ZData;
            b(k).CData = zdata;
            b(k).FaceColor = 'interp';
        end
        zlabel('correlation');
        set(gca,'XTick',1:oricombn,'XTickLabel',comblabels,...
            'YTick',1:oricombn,'YTickLabel',comblabels);
end


%% write data to file, if didn't exist before
if outexists
    fprintf('Data not saved, as %s already exists\n', outfile);
else
    fprintf('Writing data to %s\n', outfile);
    save(outfile, 'Iincr_mu', 'Iincr_var', 'Iincr_samples', ...
        'T', 'mu', 'S', 'ds', 'subdims');
end



function [mu, S, ds, T, Nmax] = datamoments(...
    sr, ori, con, ori1, ori2, con1, subdims, T)
%% returns mean, covariance, and orientation difference

% moments for requested trials
if nargin < 8, T = Inf; end
sr1 = sr(ori == ori1 & con == con1, :);
sr2 = sr(ori == ori2 & con == con1, :);
T = min([T size(sr1, 1) size(sr2, 1)]);
sr1 = sr1(1:T, :);
sr2 = sr2(1:T, :);
ds = abs(ori1 - ori2);
ds = min(ds,abs(ds-360)) * pi / 180; % angular distance
mu = (mean(sr1, 1) - mean(sr2, 1)) / ds;
S = 0.5 * cov(sr1) + 0.5 * cov(sr2);
N = length(mu);
Nmax = min(2*T-3, N);  % avoid infinite var estimator
% subdimensions, if requested
if subdims < N
    % order eigenvectors by decreasing order of eigenvalues
    [Q,Lam] = eig(S);
    [~,i] = sort(diag(Lam),'descend');
    Q = Q(:,i);
    % map into lower-dimensional subspace
    QM = Q(:,1:subdims) * Q(:,1:subdims)';
    mu = mu * QM;
    S = QM * S * QM;
    Nmax = min(subdims, Nmax);
end
