function computeOrderings(dataset, doriid, conid, shuffle, removeF)
%% computes optimal orderings for a given dataset, ori-diff, and contrast
%
% The shuffle and removeF arguments are optional. If shuffle is provided
% and true, then the data is trial-shuffled before analysis. If removeF is
% given, then it specifies the number of low-F neurons to remove before
% performing the analysis.
%
% The resulting ordering are written to
% reordering_cache/dataset_dori[dori]_c[conid][_shuf][_Fxx].mat
% The '_shuf' suffix is only added if shuffling was requested. The '_Fxx'
% is only added if removeF is provided and positive.

%% settings
randsamples = 1000;
shuffle = nargin >= 4 & logical(shuffle);
if nargin < 5, removeF = 0; end

%% load and pre-process data
d = dataInfo(dataset);
[spks, oris, cons, ~, F] = loaddata(dataset);
uoris = unique(oris);
ucons = unique(cons);
conn = length(ucons);
assert(length(uoris) == length(d.oris));
assert(conn == length(d.cons));
assert(sum(ucons == d.cons(conid)) == 1);
N = size(spks,2);
contrials = cons == d.cons(conid);
% pick desired orientation combinations
oricombs = d.oricomb(1:2, d.oricomb(3,:) == d.doris(doriid));
oricombn = size(oricombs,2);


%% pre-compute population statistics for different discriminations
highF = true(1, N);
if removeF > 0
    fprintf('Removing %d lowest-F neurons ...\n', removeF);
    [~, lowF] = sort(F);
    highF(lowF(1:removeF)) = false;
    N = N - removeF;
end
fp = NaN(oricombn, N);
Sig = NaN(N, N, oricombn);
ds = NaN(1, oricombn);
Ts = NaN(1, oricombn);
fprintf('Pre-computing population statistics ...\n');
for dorii = 1:oricombn
    ori1 = d.oris(oricombs(1,dorii));
    ori2 = d.oris(oricombs(2,dorii));
    % collect activity data for discrimination, ensuring same trial number
    ds(dorii) = min(abs(ori1-ori2),abs(ori2-ori1)) * pi/180;
    spk1 = spks((oris == ori1) & contrials, :);
    spk2 = spks((oris == ori2) & contrials, :);
    % trial-shuffle per stimulus, if requested
    if shuffle
        T1 = size(spk1,1);  T2 = size(spk2,1);
        for ni = 2:size(spk1,2)
            spk1(:,ni) = spk1(randperm(T1),ni);
            spk2(:,ni) = spk2(randperm(T2),ni);
        end
    end
    Ts(dorii) = min(size(spk1,1), size(spk2,1));
    spk1 = spk1(1:Ts(dorii),:);
    spk2 = spk2(1:Ts(dorii),:);
    fpfull = (mean(spk1,1) - mean(spk2,1)) / ds(dorii);
    Sigfull = 0.5 * (cov(spk1) + cov(spk2));
    fp(dorii,:) = fpfull(highF);
    Sig(:,:,dorii) = Sigfull(highF, highF);    
end


%% iterate over orientation combinations and compute optimal orders
In2var = @(In, T, ds) ...
    (2*In.^2 + 8*(2*T-3)/(T*ds^2)*In + 8*(2*T-3)/(T*ds^2)^2*(1:N)) ./ (2*T-(1:N)-3);
norder = NaN(oricombn+1, N);
In = NaN(oricombn, N, oricombn+1);     % discr, ~, ordering
In_var = NaN(oricombn, N, oricombn+1);  
In_ind = NaN(oricombn, N, oricombn+1);
for dorii = 1:oricombn
    ori1 = d.oris(oricombs(1,dorii));
    ori2 = d.oris(oricombs(2,dorii));
    fprintf('Computing optimal neuron order for %d vs. %d ...\n', ori1, ori2);
    % determine optimal population order
    [norder(dorii,:), In(dorii,:,dorii)] = ...
        greedyPopOrder(fp(dorii,:), Sig(:,:,dorii));
    In_var(dorii,:,dorii) = In2var(In(dorii,:,dorii), Ts(dorii), ds(dorii));
end
fprintf('Computing optimal neural order for average discrimination ...\n');
norder(oricombn+1,:) = greedyPopOrder(fp, Sig);


%% fill in information scaling for other orderings across discriminations
fprintf('Computing other info scalings for ordering ');
for dorii = 1:(oricombn+1)
    fprintf(' %d', dorii);
    % compute information scaling for other discriminations
    for dorij = 1:oricombn
        % independent information with dorii ordering
        In_ind(dorij,:,dorij) = fp(dorij, norder(dorii,:)).^2 ./ ...
            diag(Sig(norder(dorii,:), norder(dorii,:), dorij))';
        if dorii == dorij, continue; end
        In(dorij,:,dorii) = empInfscaling(fp(dorij,:), Sig(:,:,dorij), norder(dorii,:));
        In_var(dorij,:,dorii) = In2var(In(dorij,:,dorii), Ts(dorij), ds(dorij));
    end
end
fprintf('\n');


%% compute samples for random orderings for different discriminations
In_rand = NaN(randsamples, N, oricombn);
for dorii = 1:oricombn
    ori1 = d.oris(oricombs(1,dorii));
    ori2 = d.oris(oricombs(2,dorii));
    fprintf('Computing %d In samples optimal neuron order for %d vs. %d ...', ...
        randsamples, ori1, ori2);
    % for ill-condinited problems, reduce dimensionality
    [Z,Lam] = eig(Sig(:,:,dorii));
    nhi = diag(Lam)' > 1e-10 * max(diag(Lam)); % dims inclusion criteria
    Lam = Lam(nhi,nhi);
    Z = Z(:,nhi);

    for j = 1:randsamples
        if mod(j,100) == 0, fprintf(' %d', j); end
        % find empirical moments with T trials
        % resample S through eigenvalues, for numerical stability
        S = Z * wishrnd(Lam, 2*(Ts(dorii)-1)) * Z' / (2*Ts(dorii)-2);
        mup = mvnrnd(fp(dorii,:), (2/(Ts(dorii)*ds(dorii)^2)) * Sig(:,:,dorii));
        % find info scaling for random order,
        % use bias correction to get desired <In>
        In_rand(j,:,dorii) = robustInfscaling(mup, S, randperm(N,N), Ts(dorii), ds(dorii));
    end
    fprintf('\n');
end


%% writing results to file
outsuffix = '';
if shuffle, outsuffix = '_shuf'; end
if removeF > 0, outsuffix = sprintf('%s_F%d', outsuffix, removeF); end
outfile = sprintf('reordering_cache%s%s_dori%d_c%d%s.mat', ...
        filesep, dataset, doriid, conid, outsuffix);
fprintf('Writing data to %s ....\n', outfile);
save(outfile, 'N', 'oricombs', 'fp', 'Sig', 'Ts', 'ds', ...
    'norder', 'In', 'In_var', 'In_ind', 'In_rand');


function In = robustInfscaling(fp, Sig, norder, T, ds)
%% version of empInfscaling that tries to handle badly scaled Sig
%
% The arguments are the same as for empInfScaling()

if nargin < 3
    N = size(Sig, 1);
    norder = 1:N;
else
    N = length(norder);
end
In = NaN(1, N);
invSig = NaN(N, N);  % to be filled incrementally

% upper precision bound for robustness
cmax = 10^10 / max(diag(Sig));  

% reorder elements in fp and Sig
fp = fp(norder);
Sig = Sig(norder, norder); 

% start with n=1 case
invSig(1,1) = 1 / Sig(1, 1);
In(1) = invSig(1, 1) * fp(1)^2;
if isinf(invSig(1,1))
    invSig(1,1) = 0;
    In(1) = 0;
end

cStore = NaN(1, N);

% recursively handle n=2..N cases
for n = 2:N
    % Block matrix pseudo-inv. (Rohde, 1995) to add n'th element to invSig
    a = Sig(1:(n-1), n);
    b = Sig(n, n);
    invSiga = invSig(1:(n-1), 1:(n-1)) * a;
    c = 1 / (b - a' * invSiga);
    cStore(n) = c;
    if c > cmax
        % by Rohde (1995), use pseudo-inv. c = 0 for b - a' invSiga < 1/cmax
        invSig(1:(n-1), 1:(n-1)) = invSig(1:(n-1), 1:(n-1));
        invSig(1:n, n) = 0;
        invSig(n, 1:(n-1)) = 0;
        In(n) = In(n-1);
    else
        d = - c * invSiga;
        invSig(1:(n-1), 1:(n-1)) = invSig(1:(n-1), 1:(n-1)) + c * (invSiga * invSiga');
        invSig(1:(n-1), n) = d;
        invSig(n, 1:(n-1)) = d';
        invSig(n, n) = c;
        % update information measure by increment
        res = fp(1:(n-1)) * invSiga - fp(n);
        In(n) = In(n-1) + c * res^2;
    end
end

% apply bias correction, if requested
if nargin == 5
    In = (2*T-(1:N)-3)/(2*T-2) .* In - (2/(T*ds^2)) * (1:N);
end
