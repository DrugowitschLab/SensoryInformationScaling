function c = plotcols(m,s,c)
%% returns the plot color for mouse m, sessions s, and contrast c
%
% s and c are optional. plotcols(m,[],c) returns the default color for
% mouse m and contrast c.

% blue  31 120 180

% low-contrast
% orange  255 127 0
% red     227 26 28

% green   51 160 44
% magenta 152  78 163

if nargin == 1 || isempty(s)
    % default, session-independent
    if m > 4 && nargin > 2
        % contrast-specific colors for mice 5 & 6
        if c == 1
            
            cs = [255 127 0; 227 26 28];
            c = cs(m-4,:) / 255;
        elseif c == 2
            cs = [51 160 44; 152 78 163];
            c = cs(m-4,:) / 255;
        else
            error('Unknown contrast %d', c);
        end
    else
        cs = [ 31 120 180;
               31 120 180;
               31 120 180;
               31 120 180;
              255 127   0;
              227  26  28];
        c = cs(m,:) / 255;
    end
else
    % session-dependent colors
    if m == 1
        cs = [31 120 180; 165 201 225];
        c = cs(s,:) / 255;
    elseif m == 2
        cs = [31 120 180; 165 201 225];
        c = cs(s,:) / 255;
    elseif m == 3
        w = (s-1)/4;
        c = ((1-w)*[31 120 180] + w*[165 201 225]) / 255;
    elseif m == 4
        w = (s-1)/6;
        c = ((1-w)*[31 120 180] + w*[165 201 225]) / 255;
    elseif m == 5
        w = 0.6*(s-1)/3;
        if nargin > 2
            cs = [255 127 0; 51 160 44];
            c = (1-w)*cs(c,:) / 255 + w;
        else
            c = (1-w)*[255 127 0] / 255 + w;
        end
    elseif m == 6
        w = 0.6*(s-1)/2;
        if nargin > 2
            cs = [227 26 28; 152 78 163];
            c = (1-w)*cs(c,:) / 255 + w;
        else
            c = (1-w)*[227 26 28] / 255 + w;
        end
    else
        error('Unknown mouse %d', m);
    end
end
