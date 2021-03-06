function d = dataInfo(dataset)
%% returns a structure containing information about the dataset

NperD = containers.Map(...
    {  'j1a',  'j1b', 'm25a', 'm25b', 'm26a', 'm26b',...
     'aj30a','aj30b','aj31a','aj31b','aj31c',...
     'aj42a','aj42b','aj42c','aj42d','aj42e'...
     'aj43a','aj43b','aj43c','aj43d','aj43e','aj43f','aj43g',...
     'aj60a','aj60b','aj60c','aj60d',...
     'aj61a','aj61b','aj61c'},...
    [   503     493     344     339     344     386 ...
        208     303     223     225     155 ...
        280     328     346     337     282 ...
        366     362     336     351     273     304     279 ...
        387     349     388     300 ...
        392     299     338]);

switch dataset
    case {'j1a','j1b'}
        d = struct(...
            'oris', [0 15 30 45 60 75 90],...
            'cons', [0.05 0.2]);
    case {'m25a','m25b','m26a','m26b',...
          'aj30a','aj30b','aj31a','aj31b','aj31c',...
          'aj42a','aj42b','aj42c','aj42d','aj42e',...
          'aj43a','aj43b','aj43c','aj43d','aj43e','aj43f','aj43g'}
        d = struct(...
            'oris', [45 90 135 180 225 270 315 360],...
            'cons', 0.1);
    case {'aj60a','aj60b','aj60c','aj60d','aj61a','aj61b','aj61c'}
        d = struct(...
            'oris', [45 90 135 180 225 270 315 360],...
            'cons', [0.1 0.25]);
    otherwise
        error('Unknown dataset %s');
end
d.N = NperD(dataset);
% differences in orientation
doris = abs(bsxfun(@minus, d.oris, d.oris'));
doris = min(doris,abs(doris-360)); % angular distance
doris = unique(doris(:))';
doris = doris(2:end);  % remove zero difference
d.doris = doris;
% possible orientation combinations
orin = length(d.oris);
oricomb = zeros(3,0);
for i = 1:(orin-1)
    od = abs(d.oris(i) - d.oris((i+1):end));
    od = min(od,abs(od-360));
    oricomb = cat(2, oricomb, [(ones(1,orin-i)*i); (i+1):orin; od]);
end
d.oricomb = oricomb;