function [spikerate, orientation, contrast, speed, F] = loaddata(dataset, subsample)
%% loads given dataset
%
% For possible real datasets, see the 'datasets' cell array in the source code.
% The m* datasets are with, and the remaining datasets without
% microstimulation.
%
% Possible simulated datasets are 'sim[x]' with x being the simid. simid's
% >= 10 are LNP simulations.
%
% If subsample is given and not none, the following trial subsampling is
% performed:
% - 'none'  : no trial subsampling
% - 'cons'  : ensure that distribution over speeds is about the same across
%             orientations across different contrasts
% - 'spd'   : ensure that distribution over speeds is about the same across
%             orientations, but not necessarity across contrasts
% - 'lospd' : perform median split on speed and return only low speed
%             trials. Otherwise the same as 'spd'
% - 'hispd' : perform median split on speed and return only high speed
%             trials. Otherwise the same as 'spd' 
% - 'shuf'  : does not perform subsampling, but instead for each ori/cons
%             combination shuffles the trial order independently across
%             neurons, breaking correlations.
% - 'loshuf': combines lospd and shuf.
% - 'hishuf': combined hispd and shuf.

%% EDITED BY MK 08/07/2020
% - 'lospdb' : perform balanced split on speed and return only low speed
%             trials. Otherwise the same as 'spd'
% - 'hispdb' : perform balance split on speed and return only high speed
%             trials. Otherwise the same as 'spd' 

%% processing arguments
if nargin < 2, subsample = 'none'; end
if any(strcmp(dataset, {'sim1', 'sim2', 'sim3', 'sim4'}))
    filename = sprintf('simData%s%s.mat', filesep, dataset);
    filetype = 1;
    F = [];
elseif strcmp(dataset, 'sim10')
    filename = ['simData' filesep ...
        'surrogate_num_neuro_1000_num_trials_1000_diff_angle_0.26_epsilon_0.0100_file_0.mat'];
    filetype = 2;
    F = [];
else
    % load datasets from files
    datafolder = ['.' filesep 'data'];
%     datafolder = '/Users/mehdi/Dropbox (HMS)/201809-Icurve_fits/data';
    
    datasets = {'j1a', 'j1b', 'm25a', 'm25b', 'm26a', 'm26b', ...
        'aj30a', 'aj30b', 'aj31a', 'aj31b', 'aj31c', ...
        'aj42a', 'aj42b', 'aj42c', 'aj42d', 'aj42e', ...
        'aj43a', 'aj43b', 'aj43c', 'aj43d', 'aj43e', 'aj43f', 'aj43g', ...
        'aj60a', 'aj60b', 'aj60c', 'aj60d', ...
        'aj61a', 'aj61b', 'aj61c'};
    datafiles = {'j1_171018.mat', 'j1_171019.mat', 'm25_170512.mat', ...
        'm25_170523.mat', 'm26_170511.mat', 'm26_170524.mat', ...
        'AJ030_190324.mat', 'AJ030_190326.mat', ...
        'AJ031_190321.mat', 'AJ031_190322.mat', 'AJ031_190325.mat', ...
        'AJ042_190522.mat', 'AJ042_190530.mat', ...
        'AJ042_190625.mat', 'AJ042_190626.mat', 'AJ042_190627.mat', ...
        'AJ043_190525.mat', 'AJ043_190530.mat', 'AJ043_190531.mat', ...
        'AJ043_190626.mat', 'AJ043_190628.mat', ...
        'AJ043_190629.mat', 'AJ043_190630.mat', ...
        'AJ060_190902.mat', 'AJ060_190904.mat', ...
        'AJ060_190905.mat', 'AJ060_190906.mat', ...
        'AJ061_190902.mat', 'AJ061_190904.mat', 'AJ061_190905', ...
        };
    Ffiles = {'', '', ...  % j1
        'm25_170512_F.mat', 'm25_170523_F.mat', ...  % m25
        'm26_170511_F.mat', 'm26_170524_F.mat', ...  % m26
        '', '', ... % aj30
        '', '', '', ...  % aj31
        'AJ042_190522_F.mat', 'AJ042_190530_F.mat', ...
        'AJ042_190625_F.mat', 'AJ042_190626_F.mat', 'AJ042_190627_F.mat', ...
        'AJ043_190525_F.mat', 'AJ043_190530_F.mat', 'AJ043_190531_F.mat', ...
        'AJ043_190626_F.mat', 'AJ043_190628_F.mat', ...
        'AJ043_190629_F.mat', 'AJ043_190630_F.mat', ...
        '', '', ...  % aj60
        '', '', ...
        '', '', '', ...  % aj61
        };
    assert(length(datasets) == length(datafiles));
    i = strcmp(dataset, datasets);
    i = sum(i .* (1:length(i)));
    if i == 0
        error('Unknown dataset %s', dataset);
    end
    filename = [datafolder filesep datafiles{i}];
    filetype = 1;
    if nargout >= 5
        Fd = load([datafolder filesep Ffiles{i}]);
        f = fields(Fd);
        F = Fd.(f{1}).F;
    end
end


%% loading and extracting data
fprintf('Loading %s...\n', filename);
d = load(filename);
% data is one 'layer' down in AJ data
spdField = 'mvSpd';
if contains(dataset, 'aj')
    f = fields(d);
    if length(f) ~= 1
        error('Unexpected number of fields in %s dataset', dataset);
    end
    d = d.(f{1});
    spdField = 'mvSpdActual';
end

if filetype == 1
    % Gaussian simulations & mouse data
    spikerate = d.deResp;
    orientation = d.visOri;
    contrast = d.visCon;
elseif filetype == 2
    % LNP simulations
    spikerate = d.data;
    orientation = round(d.stimulus*180/pi)';
    contrast = ones(size(orientation)) * 0.2;
end

% add movement speed, if present
if isfield(d, spdField)
    speed = d.(spdField);
else
    speed = zeros(size(orientation));
end
if any(any(isnan(spikerate)))
    warning('Found NaNs in spikerate matrix');
end


%% performing subsampling/shuffling
switch subsample
    case 'none'
        fprintf('No trial subsampling\n');
        return
    case 'cons'
        % equalize speed distribution across all con/ori combinations
        ucons = unique(contrast);  conn = length(ucons);
        uoris = unique(orientation);  orin = length(uoris);
        condspds = cell(1, conn*orin);
        condtrials = cell(1, conn*orin);
        for coni = 1:conn
            for orii = 1:orin
                condi = (coni-1)*orin + orii;
                condtrials{condi} = find((orientation == uoris(orii)) ...
                    & (contrast == ucons(coni)));
                condspds{condi} = speed(condtrials{condi});
            end
        end
        rmtrials = equalizedists(condspds);
        keeptrials = [];
        for coni = 1:conn
            for orii = 1:orin
                condi = (coni-1)*orin + orii;
                keeptrials = cat(1, keeptrials, ...
                    condtrials{condi}(~rmtrials{condi}));
            end
        end
    case {'spd','lospd','hispd','loshuf','hishuf', 'lospdb', 'hispdb'}
        % perform median split by speed, if required
        if any(strcmp(subsample,{'lospd','loshuf', 'lospdb'}))
            spdtrials = speed < prctile(speed,50);
        elseif any(strcmp(subsample,{'hispd','hishuf', 'hispdb'}))
            spdtrials = speed >= prctile(speed,50);            

        else
            spdtrials = true(size(speed));
        end
        
        % equalize speed distribution across ori's, separately for cons
        ucons = unique(contrast);  conn = length(ucons);
        uoris = unique(orientation);  orin = length(uoris);
        keeptrials = [];
        for coni = 1:conn
            condspds = cell(1, orin);
            condtrials = cell(1, orin);
            for orii = 1:orin
                condtrials{orii} = find((orientation == uoris(orii)) ...
                    & (contrast == ucons(coni)) & spdtrials);
                condspds{orii} = speed(condtrials{orii});
            end
            rmtrials = equalizedists(condspds);
            for orii = 1:orin
                keeptrials = cat(1, keeptrials, ...
                    condtrials{orii}(~rmtrials{orii}));
            end
        end
        if any(strcmp(subsample,{'loshuf','hishuf'}))
            fprintf('Shuffling trials within each condition\n');
            [T, N] = size(spikerate);
            keepi = false(T, 1);
            keepi(keeptrials) = true;
            ucons = unique(contrast);  conn = length(ucons);
            uoris = unique(orientation);  orin = length(uoris);
            for coni = 1:conn
                for orii = 1:orin
                    condtrials = find((orientation == uoris(orii)) & ...
                        (contrast == ucons(coni)) & keepi);
                    T = length(condtrials);
                    for n = 2:N
                        spikerate(condtrials,n) = ...
                            spikerate(condtrials(randperm(T)),n);
                    end
                end
            end
        end
    case 'shuf'
        % shuffle trial order for each neuron/condition independently
        fprintf('Shuffling trials within each condition\n');
        N = size(spikerate,2);
        ucons = unique(contrast);  conn = length(ucons);
        uoris = unique(orientation);  orin = length(uoris);
        for coni = 1:conn
            for orii = 1:orin
                condtrials = find((orientation == uoris(orii)) & ...
                    (contrast == ucons(coni)));
                T = length(condtrials);
                for n = 2:N
                    spikerate(condtrials,n) = ...
                        spikerate(condtrials(randperm(T)),n);
                end
            end
        end
        return
    otherwise
        error('Unknown subsampling type %s', subsample)
end
fprintf('Subsampling %s, removing %4.1f%% trials (%d of %d)\n',...
    subsample, 100*(length(speed)-length(keeptrials))/length(speed), ...
    length(speed)-length(keeptrials), length(speed));
keeptrials = sort(keeptrials);
spikerate = spikerate(keeptrials,:);
orientation = orientation(keeptrials);
contrast = contrast(keeptrials);
speed = speed(keeptrials);


function rmtrials = equalizedists(x)
%% takes the cell array x of sample vectors returns equalization vector
%
% x is a cell array that contains sample vectors. These will be binned,
% after which the number of trials in each bin will be adjusted to the
% largest share by all sample vectors. The trials to be removed (sampled in
% regular order) are indicated for each sample vector in the cell array
% rmtrials.

% collect data stats to determine bin number
K = length(x);
xmin = Inf;
xmax = -Inf;
binnumk = NaN(1, K);
for k = 1:K
    xmin = min(xmin, min(x{k}));
    xmax = max(xmax, max(x{k}));
    binnumk(k) = fdbinnum(x{k});
end
binnum = ceil(mean(binnumk));

% determine bin edges and find bin distributions
c = NaN(K, binnum);
y = linspace(xmin, xmax, binnum+1);
y(end) = Inf;  % avoid boundary
for k = 1:K
    for bini = 1:binnum
        c(k, bini) = sum((x{k} >= y(bini)) & (x{k} < y(bini+1)));
    end
end

% determine which trials to remove for bin equalization
crm = bsxfun(@minus, c, min(c,[],1));  % number of bins to remove
rmtrials = cell(1, K);
for k = 1:K
    tk = false(1, length(x{k}));
    for bini = 1:binnum
        % find which trials to remove for bini in kth sample vector
        ti = false(1, c(k, bini));
        if crm(k, bini) == 1
            ti(round((length(ti)-1)/2)+1) = true;
        elseif crm(k, bini) > 1
            ti(round(linspace(1, length(ti), crm(k,bini)))) = true;
        end
        % assign these trials to complete sample vector
        tk((x{k} >= y(bini)) & (x{k} < y(bini+1))) = ti;
    end
    rmtrials{k} = tk;
end


function b = fdbinnum(x)
%% optimal number of bins by Freedman-Diaconis rule
b = ceil(range(x)*length(x)^(1/3) / (2*iqr(x)));