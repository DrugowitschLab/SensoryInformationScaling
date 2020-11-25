function s = subplotcm(pos, u, axisonly)
% creates new subplot with pos = [left bottom width height] in cm's (or u's)
if nargin < 3, axisonly = false; end
if nargin < 2, u = 'centimeters'; end
uprev = get(gcf, 'Units');
set(gcf, 'Units', u);
fpos = get(gcf, 'Position');
set(gcf, 'Units', uprev);
if axisonly
    s = axes('Position', ...
        [(pos(1)/fpos(3)) (pos(2)/fpos(4)) (pos(3)/fpos(3)) (pos(4)/fpos(4))]);
else
    s = subplot('Position', ...
        [(pos(1)/fpos(3)) (pos(2)/fpos(4)) (pos(3)/fpos(3)) (pos(4)/fpos(4))]);
end