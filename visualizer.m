%% visualize_filling.m
% Shows the vertical filling map c_bar as a before/after plot.
% Left:  5 labeled points in the domain (clean unit square, no slits)
% Right: where c_bar sends those points in the slit domain X
%
% Requires in same directory:
%   hl_draw_slits_iter.m
%   hl_overlay_grid_with_slit_holes.m
%   hl_build_graph_from_points.m
%   gd_get_row.m

clear; clc;

%% Build domain
N    = 2;
MaxN = 6;
r    = 1./(1:(N+2));

slits = hl_draw_slits_iter(r, N, 'DrawGrid', false);
pts   = hl_overlay_grid_with_slit_holes(slits, N, 'MaxN', MaxN);
G     = hl_build_graph_from_points(pts, slits);
n     = size(pts, 1);
close all;

%% Boundary and interior
tol = 1e-10;
is_boundary = (abs(pts(:,1)) < tol) | (abs(pts(:,1)-1) < tol) | ...
              (abs(pts(:,2)) < tol) | (abs(pts(:,2)-1) < tol);
boundary_idx = find(is_boundary);
interior_idx = find(~is_boundary);
top_bnd    = boundary_idx(abs(pts(boundary_idx,2) - 1) < tol);
bottom_bnd = boundary_idx(abs(pts(boundary_idx,2))     < tol);

%% Build vertical filling
dist_cache = containers.Map('KeyType','int32','ValueType','any');

bt_pairs = zeros(length(interior_idx), 2);
for i = 1:length(interior_idx)
    u  = interior_idx(i);
    xu = pts(u,1);
    [~, ti] = min(abs(pts(top_bnd,1)    - xu));
    [~, bi] = min(abs(pts(bottom_bnd,1) - xu));
    bt_pairs(i,:) = [bottom_bnd(bi), top_bnd(ti)];
end

unique_bnodes = unique(bt_pairs(:,1));
for k = 1:length(unique_bnodes)
    gd_get_row(unique_bnodes(k), G, dist_cache);
end

path_cache = containers.Map('KeyType','char','ValueType','any');
F = zeros(n, 1);
for i = 1:length(boundary_idx)
    F(boundary_idx(i)) = boundary_idx(i);
end
for i = 1:length(interior_idx)
    u      = interior_idx(i);
    yu     = pts(u,2);
    b_node = bt_pairs(i,1);
    t_node = bt_pairs(i,2);
    pkey   = sprintf('%d_%d', b_node, t_node);

    if ~path_cache.isKey(pkey)
        pnodes = shortestpath(G, b_node, t_node);
        b_row  = dist_cache(int32(b_node));
        plen   = b_row(t_node);
        path_cache(pkey)          = pnodes;
        path_cache([pkey '_len']) = plen;
    end

    path_nodes  = path_cache(pkey);
    path_len    = path_cache([pkey '_len']);
    target_dist = yu * path_len;

    cum_dist = 0;
    chosen   = path_nodes(1);
    for k = 2:length(path_nodes)
        cum_dist = cum_dist + norm(pts(path_nodes(k),:) - pts(path_nodes(k-1),:));
        if cum_dist >= target_dist
            chosen = path_nodes(k);
            break
        end
    end
    F(u) = chosen;
end

%% Pick 5 random interior nodes
rng(42);
sample = interior_idx(randperm(length(interior_idx), 5));
labels = 'ABCDE';

%% Plot
figure('Position', [100 100 1000 460]); 

% --- Left: clean unit square, no slits ---
subplot(1,2,1); hold on; axis equal;
xlim([0 1]); ylim([0 1]);
title('Domain Q (pre-image)');
xlabel('x'); ylabel('y');
rectangle('Position', [0 0 1 1], 'EdgeColor', 'k', 'LineWidth', 1.5);

for i = 1:5
    u = sample(i);
    plot(pts(u,1), pts(u,2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    text(pts(u,1)+0.02, pts(u,2)+0.02, labels(i), ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
end
hold off;

% --- Right: slit domain X ---
subplot(1,2,2); hold on; axis equal;
xlim([0 1]); ylim([0 1]);
title(sprintf('Slit domain X: image \\bar{c}(u) (N=%d)', N));
xlabel('x'); ylabel('y');

for k = 1:numel(slits)
    plot([slits(k).x slits(k).x], [slits(k).y0 slits(k).y1], ...
        'k-', 'LineWidth', 2);
end

for i = 1:5
    u  = sample(i);
    fu = F(u);
    plot(pts(fu,1), pts(fu,2), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
    text(pts(fu,1)+0.02, pts(fu,2)+0.02, labels(i), ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');
end
hold off;