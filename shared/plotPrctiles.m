function plotPrctiles(x, ss, col, cp, thickp, thinp)
%% plots single percentile line plot of samples ss at x
%
% The function plots a thin line covering the thinp credible interval
% (defaults to 90), a thicker line covering the thickp credible interval
% (defaults to 50), and a dot at the cp percentile (defaults to 50), with
% color col (defaults to black).

%% process arguments
if nargin < 6
    thinp = 90;
    if nargin < 5
        thickp = 50;
        if nargin < 4
            cp = 50;
            if nargin < 3
                col = [0 0 0];
            end
        end
    end
end

%% plot lines
plot([x x], [prctile(ss,(100-thinp)/2) prctile(ss,100-(100-thinp)/2)], ...
    '-', 'LineWidth', 0.5, 'Color', col);
plot([x x], [prctile(ss,(100-thickp)/2) prctile(ss,100-(100-thickp)/2)], ...
    '-', 'LineWidth', 2, 'Color', col);
plot(x, prctile(ss,cp), 'o', 'MarkerSize', 6, ...
    'MarkerFaceColor', col, 'MarkerEdgeColor', 'none');
