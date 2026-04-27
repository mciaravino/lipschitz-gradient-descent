%% gd_top.m
% Coordinate descent to minimize Lip(c_bar_N) starting from top filling.
% Captures a video frame at each iteration where Lip improves.
%
% Requires:
%   hl_draw_slits_iter.m
%   hl_overlay_grid_with_slit_holes.m
%   hl_build_domain_graph.m
%   hl_build_graph_from_points.m
%   hl_top_filling.m
%   gd_get_row.m
%   gd_dlookup.m
%   global_lip.m
%   draw_colored_arrows.m

clear; clc;

%% ---- Parameters ----
N         = 3;
MaxN      = 6;
max_iters = 1000;
%% -------------------

r   = 1./(1:(N+2));
tol = 1e-10;

%% Build geometry
slits = hl_draw_slits_iter(r, N, 'DrawGrid', false);
close all;

fprintf('Building domain graph G_D...\n');
[~, pts_D] = hl_build_domain_graph(MaxN);

fprintf('Building slit domain graph G_X...\n');
pts_X = hl_overlay_grid_with_slit_holes(slits, N, 'MaxN', MaxN);
G_X   = hl_build_graph_from_points(pts_X, slits);
close all;

fprintf('Domain nodes: %d  |  X_N nodes: %d\n', size(pts_D,1), size(pts_X,1));

%% Build top filling
fprintf('Building top filling...\n');
%[F, ~] = hl_top_filling(pts_D, slits, G_X, pts_X);
[F, ~] = hl_lane_filling(pts_D, slits, G_X, pts_X);
F_init = F;

%% Node sets
is_bnd  = (abs(pts_D(:,1)) < tol) | (abs(pts_D(:,1)-1) < tol) | ...
          (abs(pts_D(:,2)) < tol) | (abs(pts_D(:,2)-1) < tol);
all_idx = (1:size(pts_D,1))';
n_all   = length(all_idx);
int_D   = find(~is_bnd);

% Only check pairs within r_max — sort by x for early break
ux_D  = unique(pts_D(:,1));
h_D   = min(diff(ux_D));
r_max = 3 * h_D;

% Sort all_idx by x for early break in pair loops
[~, ord] = sort(pts_D(all_idx, 1));
all_idx  = all_idx(ord);
n_all    = length(all_idx);

% Pre-build radius neighbor list for each node (used in descent loops)
fprintf('Building radius neighbor lists...\n');
radius_nbrs = cell(size(pts_D,1), 1);
for a = 1:n_all
    ua = all_idx(a);
    nbrs = [];
    for b = 1:n_all
        vb = all_idx(b);
        if vb == ua, continue; end
        if pts_D(vb,1) - pts_D(ua,1) > r_max, break; end
        if abs(pts_D(ua,2) - pts_D(vb,2)) > r_max, continue; end
        d_dom = norm(pts_D(ua,:) - pts_D(vb,:));
        if d_dom < 1e-12 || d_dom > r_max, continue; end
        nbrs(end+1) = vb;
    end
    radius_nbrs{ua} = nbrs;
end
fprintf('  Done.\n');

%% Build adjacency list for G_X
n_X   = size(pts_X,1);
E_X   = G_X.Edges.EndNodes;
adj_X = cell(n_X, 1);
for k = 1:size(E_X,1)
    u = E_X(k,1); v = E_X(k,2);
    adj_X{u} = [adj_X{u}, v];
    adj_X{v} = [adj_X{v}, u];
end

%% Distance cache
dist_cache = containers.Map('KeyType','int32','ValueType','any');
fprintf('Caching image-node distance rows...\n');
tic
for nd = unique(F)'
    gd_get_row(nd, G_X, dist_cache);
end
fprintf('  Done: %.3f sec\n', toc);

%% Compute initial Lip
fprintf('Computing initial Lip...\n');
tic
[lip_init, wu, wv] = global_lip(F, all_idx, pts_D, G_X, dist_cache);
fprintf('  Done: %.3f sec\n', toc);
fprintf('Initial Lip(c_bar_%d) = %.4f\n', N, lip_init);

%% Video setup
%vid_file = sprintf('descent_N%d.mp4', N);
vid_file = sprintf('descent_lane_N%d.mp4', N);
v = VideoWriter(vid_file, 'MPEG-4');
v.FrameRate = 3;
open(v);
fprintf('Recording video: %s\n', vid_file);

%% Open figure for video
figure('Position', [100 100 700 650]);

%% Capture frame 0 — initial filling
capture_frame(v, F, 0, lip_init, wu, wv, N, pts_D, pts_X, slits, is_bnd, G_X);

%% Coordinate descent
fprintf('\nRunning coordinate descent...\n');
lip_history = [lip_init];
lip_curr    = lip_init;

for iter = 1:max_iters

    % Find all optimizable nodes at worst ratio
    worst_nodes = [];
    for a = 1:n_all
        for b = a+1:n_all
            ua = all_idx(a); vb = all_idx(b);
            if pts_D(vb,1) - pts_D(ua,1) > r_max, break; end
            if abs(pts_D(ua,2) - pts_D(vb,2)) > r_max, continue; end
            d_dom = norm(pts_D(ua,:) - pts_D(vb,:));
            if d_dom < 1e-12 || d_dom > r_max, continue; end
            d_img = gd_dlookup(F(ua), F(vb), G_X, dist_cache);
            if isinf(d_img), continue; end
            if abs(d_img/d_dom - lip_curr) < 1e-6
                if ~is_bnd(ua), worst_nodes = [worst_nodes; ua]; end
                if ~is_bnd(vb), worst_nodes = [worst_nodes; vb]; end
            end
        end
    end
    worst_nodes = unique(worst_nodes);

    if isempty(worst_nodes)
        fprintf('  No optimizable nodes at worst ratio — converged.\n');
        break
    end

    fprintf('  Iter %3d: %d nodes at worst ratio %.4f\n', ...
        iter, length(worst_nodes), lip_curr);

    % Sweep
    n_swaps   = 0;
    idx_order = worst_nodes(randperm(length(worst_nodes)));

    for ii = 1:length(idx_order)
        u  = idx_order(ii);
        fu = F(u);

        current_worst = 0;
        nbrs_u = radius_nbrs{u};
        for bi = 1:length(nbrs_u)
            vb = nbrs_u(bi);
            d_dom = norm(pts_D(u,:) - pts_D(vb,:));
            if d_dom < 1e-12, continue; end
            d_img = gd_dlookup(F(u), F(vb), G_X, dist_cache);
            if isinf(d_img), continue; end
            current_worst = max(current_worst, d_img / d_dom);
            if current_worst >= lip_curr, break; end
        end

        candidates = [fu, adj_X{fu}];
        best_cost  = current_worst;
        best_node  = fu;

        for c = 1:length(candidates)
            cand = candidates(c);
            gd_get_row(cand, G_X, dist_cache);
            F(u) = cand;

            tentative_worst = 0;
            for bi = 1:length(nbrs_u)
                vb = nbrs_u(bi);
                d_dom = norm(pts_D(u,:) - pts_D(vb,:));
                if d_dom < 1e-12, continue; end
                d_img = gd_dlookup(cand, F(vb), G_X, dist_cache);
                if isinf(d_img), continue; end
                tentative_worst = max(tentative_worst, d_img / d_dom);
                if tentative_worst >= best_cost, break; end
            end

            if tentative_worst < best_cost
                best_cost = tentative_worst;
                best_node = cand;
            end
        end

        F(u) = best_node;
        if best_node ~= fu
            n_swaps = n_swaps + 1;
        end
    end

    % Recompute global Lip
    [lip_new, wu_new, wv_new] = global_lip(F, all_idx, pts_D, G_X, dist_cache);
    fprintf('           -> Lip = %.4f  (swaps: %d)\n', lip_new, n_swaps);

    if lip_new < lip_curr - 1e-9
        lip_curr = lip_new;
        wu = wu_new; wv = wv_new;
        lip_history(end+1) = lip_curr;
        % Capture frame for this iteration
        capture_frame(v, F, iter, lip_curr, wu, wv, N, pts_D, pts_X, slits, is_bnd, G_X);
    else
        fprintf('  Converged at iteration %d.\n', iter);
        break
    end
end

%% Close video
close(v);
fprintf('Video saved: %s\n', vid_file);

%% Results
fprintf('\n========== Results (N=%d) ==========\n', N);
fprintf('  Initial Lip: %.4f\n', lip_init);
fprintf('  Final   Lip: %.4f\n', lip_curr);
fprintf('  Improvement: %.4f\n', lip_init - lip_curr);

%% Convergence plot
figure;
plot(0:length(lip_history)-1, lip_history, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Iteration');
ylabel('Lip(\bar{c}_N)');
title(sprintf('\\bar{c}_%d convergence  (MaxN=%d)', N, MaxN), 'Interpreter', 'tex');
grid on;