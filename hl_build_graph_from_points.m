function G = hl_build_graph_from_points(pts, slits)
% HL_BUILD_GRAPH_FROM_POINTS  Build graph from grid points, respecting slits.
%
%   G = hl_build_graph_from_points(pts, slits)
%
%   Connects each node to its 8 grid neighbors by Euclidean distance,
%   but rejects any edge whose segment crosses a slit obstacle.
%
%   Inputs:
%     pts   - n x 2 array of node coordinates
%     slits - struct array with fields x, y0, y1 (from hl_draw_slits_iter)
%
%   Output:
%     G     - MATLAB graph object with edge weights = Euclidean distances

n = size(pts, 1);

% Build coordinate lookup table
key = containers.Map;
for i = 1:n
    k = sprintf('%.10f_%.10f', pts(i,1), pts(i,2));
    key(k) = i;
end

% Grid spacing
dx = unique(diff(unique(pts(:,1))));
dy = unique(diff(unique(pts(:,2))));
h  = min([dx; dy]);

% 8-neighbor offsets
nbr = [-1 -1; -1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0; 1 1];

% Precompute slit data for fast crossing checks
slit_x  = [slits.x]';
slit_y0 = [slits.y0]';
slit_y1 = [slits.y1]';

S = []; T = []; W = [];

for i = 1:n
    x = pts(i,1);
    y = pts(i,2);

    for k = 1:8
        xn = x + nbr(k,1) * h;
        yn = y + nbr(k,2) * h;
        keystr = sprintf('%.10f_%.10f', xn, yn);

        if ~isKey(key, keystr), continue; end
        j = key(keystr);
        if j <= i, continue; end

        % Reject edge if it crosses any slit
        if edge_crosses_slit(x, y, xn, yn, slit_x, slit_y0, slit_y1)
            continue;
        end

        S(end+1) = i;
        T(end+1) = j;
        W(end+1) = hypot(x - xn, y - yn);
    end
end

G = graph(S, T, W, n);
fprintf('Graph edges: %d\n', numedges(G));
end

% -------------------------------------------------------
function crosses = edge_crosses_slit(x1, y1, x2, y2, slit_x, slit_y0, slit_y1)
% Returns true if the segment (x1,y1)-(x2,y2) crosses any slit.
% A slit is a vertical segment at x = slit_x(k), from y0 to y1.
% Crossing occurs when the edge passes through the slit's x-coordinate
% and the intersection y-value falls within the slit's y-range.

crosses = false;

for k = 1:length(slit_x)
    sx  = slit_x(k);
    sy0 = slit_y0(k);
    sy1 = slit_y1(k);

    % Edge must straddle the slit x-coordinate
    if (x1 < sx && x2 > sx) || (x1 > sx && x2 < sx)
        % Interpolate y at x = sx
        t  = (sx - x1) / (x2 - x1);
        yi = y1 + t * (y2 - y1);

        % Check if intersection is within slit y-range
        if yi >= sy0 && yi <= sy1
            crosses = true;
            return;
        end
    end
end
end