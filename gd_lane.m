%% gd_lane.m
% Coordinate descent to minimize Lip(c_bar_N).
%
% Algorithm:
%   1. Compute global Lip and find all pairs achieving that ratio
%   2. Collect all nodes involved in those worst pairs
%   3. Sweep over those nodes, trying to move each image to a neighbor
%      that reduces its local worst ratio
%   4. Recompute global Lip
%   5. If improved: find new worst-ratio nodes, repeat
%   6. If not improved: converged
%
% Requires:
%   hl_draw_slits_iter.m
%   hl_overlay_grid_with_slit_holes.m
%   hl_build_domain_graph.m
%   hl_build_graph_from_points.m
%   hl_lane_filling.m
%   gd_get_row.m
%   gd_dlookup.m
%   global_lip.m

clear; clc;

%% ---- Parameters ----
N         = 2;
MaxN      = 4%N + 2;
max_iters = 100;
%% -------------------

r   = 1./(1:(N+2));
tol = 1e-10;

%% Build geometry
slits = hl_draw_slits_iter(r, N, 'DrawGrid', false);
close all;

fprintf('Building domain graph G_D...\n');
[G_D, pts_D] = hl_build_domain_graph(MaxN);

fprintf('Building slit domain graph G_X...\n');
pts_X = hl_overlay_grid_with_slit_holes(slits, N, 'MaxN', MaxN);
G_X   = hl_build_graph_from_points(pts_X, slits);
close all;

fprintf('Domain nodes: %d  |  X_N nodes: %d\n', size(pts_D,1), size(pts_X,1));

%% Build lane filling
fprintf('Building lane filling...\n');
%[F, ~] = hl_lane_filling(pts_D, slits, G_X, pts_X);
[F, ~] = hl_top_filling(pts_D, slits, G_X, pts_X);
F_init = F;

%% Slit-column restriction for Lip evaluation
ux_D = unique(pts_D(:,1));
h_D  = min(diff(ux_D));
is_bnd = (abs(pts_D(:,1)) < tol) | (abs(pts_D(:,1)-1) < tol) | ...
         (abs(pts_D(:,2)) < tol) | (abs(pts_D(:,2)-1) < tol);

keep = false(size(pts_D,1), 1);
for k = 1:numel(slits)
    in_x = abs(pts_D(:,1) - slits(k).x) <= h_D * 1.5;
    in_y = pts_D(:,2) >= slits(k).y0 - h_D & pts_D(:,2) <= slits(k).y1 + h_D;
    keep = keep | (in_x & in_y);
end
keep_idx = find(keep);


n_lip    = length(keep_idx);
fprintf('Lip evaluation nodes: %d\n', n_lip);

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
[lip_init, wu_init, wv_init] = global_lip(F, keep_idx, pts_D, G_X, dist_cache);
fprintf('  Done: %.3f sec\n', toc);
fprintf('Initial Lip(c_bar) = %.4f\n', lip_init);

%% Coordinate descent
fprintf('\nRunning coordinate descent...\n');
lip_history = [lip_init];
lip_curr    = lip_init;
wu_curr     = wu_init;
wv_curr     = wv_init;

for iter = 1:max_iters

    % Find ALL pairs achieving the current worst ratio
    worst_nodes = [];
    for a = 1:n_lip
        for b = a+1:n_lip
            ua = keep_idx(a); vb = keep_idx(b);
            d_dom = norm(pts_D(ua,:) - pts_D(vb,:));
            if d_dom < 1e-12, continue; end
            d_img = gd_dlookup(F(ua), F(vb), G_X, dist_cache);
            if isinf(d_img), continue; end
            if abs(d_img/d_dom - lip_curr) < 1e-6
                worst_nodes = [worst_nodes; ua; vb];
            end
        end
    end
    worst_nodes = unique(worst_nodes);
    fprintf('  Iter %3d: %d nodes at worst ratio %.4f\n', ...
        iter, length(worst_nodes), lip_curr);

    % Sweep over worst-ratio nodes in random order
    n_swaps   = 0;
    idx_order = worst_nodes(randperm(length(worst_nodes)));

    for ii = 1:length(idx_order)
        u  = idx_order(ii);
        fu = F(u);

        % Current local worst ratio for this node over all lip pairs
        current_worst = 0;
        for b = 1:n_lip
            vb = keep_idx(b);
            if vb == u, continue; end
            d_dom = norm(pts_D(u,:) - pts_D(vb,:));
            if d_dom < 1e-12, continue; end
            d_img = gd_dlookup(F(u), F(vb), G_X, dist_cache);
            if isinf(d_img), continue; end
            current_worst = max(current_worst, d_img / d_dom);
        end

        % Try neighbors of F(u) in G_X
        candidates = [fu, adj_X{fu}];
        best_cost  = current_worst;
        best_node  = fu;

        for c = 1:length(candidates)
            cand = candidates(c);
            gd_get_row(cand, G_X, dist_cache);
            F(u) = cand;

            tentative_worst = 0;
            for b = 1:n_lip
                vb = keep_idx(b);
                if vb == u, continue; end
                d_dom = norm(pts_D(u,:) - pts_D(vb,:));
                if d_dom < 1e-12, continue; end
                d_img = gd_dlookup(cand, F(vb), G_X, dist_cache);
                if isinf(d_img), continue; end
                tentative_worst = max(tentative_worst, d_img / d_dom);
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
    [lip_new, wu_new, wv_new] = global_lip(F, keep_idx, pts_D, G_X, dist_cache);
    fprintf('           -> Lip = %.4f  (swaps: %d)\n', lip_new, n_swaps);

    if lip_new < lip_curr - 1e-9
        lip_curr = lip_new;
        wu_curr  = wu_new;
        wv_curr  = wv_new;
        lip_history(end+1) = lip_curr;
    else
        fprintf('  Converged at iteration %d.\n', iter);
        break
    end
end

%% Results
fprintf('\n========== Results (N=%d) ==========\n', N);
fprintf('  Initial Lip: %.4f\n', lip_init);
fprintf('  Final   Lip: %.4f\n', lip_curr);
fprintf('  Improvement: %.4f\n', lip_init - lip_curr);

[~, wu_final, wv_final] = global_lip(F, keep_idx, pts_D, G_X, dist_cache);

%% Convergence plot
figure;
plot(0:length(lip_history)-1, lip_history, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Iteration');
ylabel('Lip(c-bar)');
title(sprintf('Coordinate descent convergence (N=%d)', N));
grid on;

%% Before/after plots
for pass = 1:2
    if pass == 1
        F_plot = F_init; wu_p = wu_init; wv_p = wv_init;
        lip_p  = lip_init; label = 'initial';
    else
        F_plot = F; wu_p = wu_final; wv_p = wv_final;
        lip_p  = lip_curr; label = 'optimized';
    end

    int_D = find(~is_bnd);
    img_p = pts_X(F_plot(int_D), :);
    dx_p  = img_p(:,1) - pts_D(int_D,1);
    dy_p  = img_p(:,2) - pts_D(int_D,2);
    mov_p = (abs(dx_p) + abs(dy_p)) > tol;

    figure('Position', [100 100 1200 500]);

    subplot(1,2,1); hold on; axis equal;
    xlim([0 1]); ylim([0 1]);
    title(sprintf('Domain Q — %s (Lip=%.4f)', label, lip_p));
    xlabel('x'); ylabel('y');
    plot(pts_D(int_D,1), pts_D(int_D,2), '.', 'Color', [0.88 0.88 0.88], 'MarkerSize', 3);
    plot(pts_D(keep_idx,1), pts_D(keep_idx,2), 'b.', 'MarkerSize', 5);
    plot(pts_D(wu_p,1), pts_D(wu_p,2), 'ro', 'MarkerSize', 14, 'LineWidth', 2.5);
    plot(pts_D(wv_p,1), pts_D(wv_p,2), 'rs', 'MarkerSize', 14, 'LineWidth', 2.5);
    hold off;

    subplot(1,2,2); hold on; axis equal;
    xlim([0 1]); ylim([0 1]);
    title(sprintf('Image c-bar(u) in X_N — %s', label));
    xlabel('x'); ylabel('y');
    for k = 1:numel(slits)
        plot([slits(k).x slits(k).x], [slits(k).y0 slits(k).y1], ...
            'k-', 'LineWidth', 2);
    end
    plot(pts_X(:,1), pts_X(:,2), '.', 'Color', [0.88 0.88 0.88], 'MarkerSize', 3);
    if any(mov_p)
        quiver(pts_D(int_D(mov_p),1), pts_D(int_D(mov_p),2), ...
               dx_p(mov_p), dy_p(mov_p), 0, ...
               'Color', [0.2 0.5 0.9], 'LineWidth', 1.2, 'MaxHeadSize', 0.5);
    end
    path_p = shortestpath(G_X, F_plot(wu_p), F_plot(wv_p));
    plot(pts_X(path_p,1), pts_X(path_p,2), 'g-', 'LineWidth', 2.5);
    plot(pts_X(F_plot(wu_p),1), pts_X(F_plot(wu_p),2), 'ro', 'MarkerSize', 14, 'LineWidth', 2.5);
    plot(pts_X(F_plot(wv_p),1), pts_X(F_plot(wv_p),2), 'rs', 'MarkerSize', 14, 'LineWidth', 2.5);
    hold off;
end