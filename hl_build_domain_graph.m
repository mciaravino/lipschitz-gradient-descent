function [G_D, pts_D] = hl_build_domain_graph(MaxN)
% HL_BUILD_DOMAIN_GRAPH  Build full uniform grid graph over [0,1]^2.
%
%   [G_D, pts_D] = hl_build_domain_graph(MaxN)
%
%   No nodes removed, no edges rejected. This is the graph for the
%   flat domain D = [0,1]^2 with no obstacles.
%   Distances on G_D equal Euclidean distances.
%
%   Input:
%     MaxN  - grid resolution parameter (same as hl_overlay_grid_with_slit_holes)
%             grid is (2^(MaxN+1) + 1) x (2^(MaxN+1) + 1) points
%
%   Outputs:
%     G_D   - MATLAB graph object, full grid
%     pts_D - n x 2 node coordinates

L      = 2^(MaxN + 1);
coords = (0:L)' / L;
[X, Y] = meshgrid(coords, coords);
pts_D  = [X(:), Y(:)];
n      = size(pts_D, 1);

% Build coordinate lookup
key = containers.Map;
for i = 1:n
    k = sprintf('%.10f_%.10f', pts_D(i,1), pts_D(i,2));
    key(k) = i;
end

h   = 1 / L;
nbr = [-1 -1; -1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0; 1 1];

S = []; T = []; W = [];
for i = 1:n
    x = pts_D(i,1);
    y = pts_D(i,2);
    for k = 1:8
        xn = x + nbr(k,1) * h;
        yn = y + nbr(k,2) * h;
        keystr = sprintf('%.10f_%.10f', xn, yn);
        if ~isKey(key, keystr), continue; end
        j = key(keystr);
        if j <= i, continue; end
        S(end+1) = i;
        T(end+1) = j;
        W(end+1) = hypot(x - xn, y - yn);
    end
end

G_D = graph(S, T, W, n);
fprintf('Domain graph: %d nodes, %d edges\n', n, numedges(G_D));
end