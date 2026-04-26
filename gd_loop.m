%% gd_loop.m
% Runs gradient descent for N = 2, 3, 4 at fixed MaxN = 6.
% Records initial and final Lip(c_bar) for each N.
% Plots L*_N vs N to show whether Lipschitz constant grows with N.
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
N_values  = 4:4;
MaxN      = 6;      % fixed for all N — keeps h constant for valid comparison
max_iters = 100;
%% -------------------

results = zeros(length(N_values), 3);  % [N, lip_init, lip_final]

for ni = 1:length(N_values)

    N = N_values(ni);
    fprintf('\n========== N = %d (MaxN = %d) ==========\n', N, MaxN);
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

    %% Build lane filling
    fprintf('Building lane filling...\n');
    [F, ~] = hl_lane_filling(pts_D, slits, G_X, pts_X);
    F_lane = F;  % save for before plot

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

    %% Initial Lip
    fprintf('Computing initial Lip...\n');
    tic
    [lip_init, wu, wv] = global_lip(F, keep_idx, pts_D, G_X, dist_cache);
    wu_init = wu; wv_init = wv;
    fprintf('  Done: %.3f sec\n', toc);
    fprintf('Initial Lip(c_bar) = %.4f\n', lip_init);

    %% Coordinate descent
    fprintf('Running coordinate descent...\n');
    lip_curr = lip_init;

    for iter = 1:max_iters

        % Find all nodes at worst ratio
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

        % Sweep over worst-ratio nodes
        n_swaps   = 0;
        idx_order = worst_nodes(randperm(length(worst_nodes)));

        for ii = 1:length(idx_order)
            u  = idx_order(ii);
            fu = F(u);

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
        [lip_new, wu, wv] = global_lip(F, keep_idx, pts_D, G_X, dist_cache);
        fprintf('  Iter %3d: %d nodes, Lip = %.4f (swaps: %d)\n', ...
            iter, length(worst_nodes), lip_new, n_swaps);

        if lip_new < lip_curr - 1e-9
            lip_curr = lip_new;
        else
            fprintf('  Converged at iteration %d.\n', iter);
            break
        end
    end

    fprintf('Final Lip(c_bar) = %.4f  (improvement: %.4f)\n', ...
        lip_curr, lip_init - lip_curr);

    results(ni,:) = [N, lip_init, lip_curr];

    %% Before/after plots for this N
    [~, wu_final, wv_final] = global_lip(F, keep_idx, pts_D, G_X, dist_cache);

    for pass = 1:2
        if pass == 1
            F_plot = F_lane; wu_p = wu_init; wv_p = wv_init;
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
        title(sprintf('Domain Q — %s  N=%d  Lip=%.4f', label, N, lip_p));
        xlabel('x'); ylabel('y');
        plot(pts_D(int_D,1), pts_D(int_D,2), '.', 'Color', [0.88 0.88 0.88], 'MarkerSize', 3);
        plot(pts_D(keep_idx,1), pts_D(keep_idx,2), 'b.', 'MarkerSize', 5);
        plot(pts_D(wu_p,1), pts_D(wu_p,2), 'ro', 'MarkerSize', 14, 'LineWidth', 2.5);
        plot(pts_D(wv_p,1), pts_D(wv_p,2), 'rs', 'MarkerSize', 14, 'LineWidth', 2.5);
        hold off;

        subplot(1,2,2); hold on; axis equal;
        xlim([0 1]); ylim([0 1]);
        title(sprintf('Image c-bar(u) in X_N — %s  N=%d', label, N));
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

end

%% Summary table
fprintf('\n========== Summary (MaxN=%d) ==========\n', MaxN);
fprintf('  N    Lip_init    Lip_final\n');
for ni = 1:size(results,1)
    fprintf('  %d    %.4f      %.4f\n', results(ni,1), results(ni,2), results(ni,3));
end

%% Summary plot
figure; hold on;
plot(results(:,1), results(:,2), 'rs--', 'LineWidth', 2, 'MarkerSize', 8);
plot(results(:,1), results(:,3), 'bo-',  'LineWidth', 2, 'MarkerSize', 8);
xlabel('N (dyadic depth)');
ylabel('Lip(c-bar)');
title(sprintf('L*_N vs N  (MaxN=%d, fixed resolution)', MaxN));
legend({'Lane filling (initial)', 'After gradient descent'}, 'Location', 'northwest');
grid on;
hold off;

% ========== N = 2 (MaxN = 6) ==========
% Building domain graph G_D...
% Domain graph: 16641 nodes, 65792 edges
% Building slit domain graph G_X...
% Grid nodes: 16348
% Graph edges: 63720
% Domain nodes: 16641  |  X_N nodes: 16348
% Building lane filling...
% Lane filling: boundary violations = 0
% Lane filling: nodes that moved    = 143 / 16641
% Lip evaluation nodes: 1005
% Caching image-node distance rows...
%   Done: 33.522 sec
% Computing initial Lip...
%   Done: 13.537 sec
% Initial Lip(c_bar) = 66.8284
% Running coordinate descent...
%   Iter   1: 2 nodes, Lip = 64.8284 (swaps: 2)
%   Iter   2: 6 nodes, Lip = 62.8284 (swaps: 6)
%   Iter   3: 10 nodes, Lip = 60.8284 (swaps: 10)
%   Iter   4: 14 nodes, Lip = 58.8284 (swaps: 14)
%   Iter   5: 18 nodes, Lip = 56.8284 (swaps: 18)
%   Iter   6: 22 nodes, Lip = 54.8284 (swaps: 22)
%   Iter   7: 26 nodes, Lip = 52.8284 (swaps: 26)
%   Iter   8: 30 nodes, Lip = 50.8284 (swaps: 30)
%   Iter   9: 34 nodes, Lip = 48.8284 (swaps: 34)
%   Iter  10: 38 nodes, Lip = 47.2548 (swaps: 38)
%   Iter  11: 4 nodes, Lip = 46.8284 (swaps: 2)
%   Iter  12: 38 nodes, Lip = 46.5477 (swaps: 38)
%   Iter  13: 4 nodes, Lip = 45.8406 (swaps: 3)
%   Iter  14: 2 nodes, Lip = 45.1335 (swaps: 1)
%   Iter  15: 4 nodes, Lip = 44.8284 (swaps: 2)
%   Iter  16: 42 nodes, Lip = 44.4264 (swaps: 42)
%   Iter  17: 4 nodes, Lip = 43.8284 (swaps: 3)
%   Iter  18: 2 nodes, Lip = 43.0122 (swaps: 1)
%   Iter  19: 4 nodes, Lip = 42.8284 (swaps: 2)
%   Iter  20: 46 nodes, Lip = 42.3051 (swaps: 46)
%   Iter  21: 4 nodes, Lip = 41.5980 (swaps: 2)
%   Iter  22: 4 nodes, Lip = 40.8909 (swaps: 2)
%   Iter  23: 4 nodes, Lip = 40.8284 (swaps: 4)
%   Iter  24: 50 nodes, Lip = 39.8284 (swaps: 50)
%   Iter  25: 4 nodes, Lip = 39.8284 (swaps: 1)
%   Converged at iteration 25.
% Final Lip(c_bar) = 39.8284  (improvement: 27.0000)
% 
% ========== Summary (MaxN=6) ==========
%   N    Lip_init    Lip_final
%   2    66.8284      39.8284


% ========== N = 3 (MaxN = 6) ==========
% Building domain graph G_D...
% Domain graph: 16641 nodes, 65792 edges
% Building slit domain graph G_X...
% Grid nodes: 16156
% Graph edges: 62312
% Domain nodes: 16641  |  X_N nodes: 16156
% Building lane filling...
% Lane filling: boundary violations = 0
% Lane filling: nodes that moved    = 167 / 16641
% Lip evaluation nodes: 1965
% Caching image-node distance rows...
%   Done: 33.607 sec
% Computing initial Lip...
%   Done: 52.478 sec
% Initial Lip(c_bar) = 66.8284
% Running coordinate descent...
%   Iter   1: 2 nodes, Lip = 64.8284 (swaps: 2)
%   Iter   2: 6 nodes, Lip = 62.8284 (swaps: 6)
%   Iter   3: 10 nodes, Lip = 60.8284 (swaps: 10)
%   Iter   4: 14 nodes, Lip = 58.8284 (swaps: 14)
%   Iter   5: 18 nodes, Lip = 56.8284 (swaps: 18)
%   Iter   6: 22 nodes, Lip = 54.8284 (swaps: 22)
%   Iter   7: 26 nodes, Lip = 52.8284 (swaps: 26)
%   Iter   8: 30 nodes, Lip = 50.8284 (swaps: 30)
%   Iter   9: 34 nodes, Lip = 48.8284 (swaps: 34)
%   Iter  10: 38 nodes, Lip = 47.2548 (swaps: 38)
%   Iter  11: 4 nodes, Lip = 46.8284 (swaps: 2)
%   Iter  12: 38 nodes, Lip = 46.5477 (swaps: 38)
%   Iter  13: 4 nodes, Lip = 45.8406 (swaps: 3)
%   Iter  14: 2 nodes, Lip = 45.1335 (swaps: 1)
%   Iter  15: 4 nodes, Lip = 44.8284 (swaps: 2)
%   Iter  16: 42 nodes, Lip = 44.4264 (swaps: 42)
%   Iter  17: 4 nodes, Lip = 43.7193 (swaps: 2)
%   Iter  18: 4 nodes, Lip = 43.0122 (swaps: 3)
%   Iter  19: 2 nodes, Lip = 42.8284 (swaps: 1)
%   Iter  20: 46 nodes, Lip = 42.3051 (swaps: 46)
%   Iter  21: 4 nodes, Lip = 41.5980 (swaps: 2)
%   Iter  22: 4 nodes, Lip = 40.8284 (swaps: 4)
%   Iter  23: 54 nodes, Lip = 39.8284 (swaps: 52)
%   Iter  24: 4 nodes, Lip = 39.4767 (swaps: 2)
%   Iter  25: 4 nodes, Lip = 39.4767 (swaps: 0)
%   Converged at iteration 25.
% Final Lip(c_bar) = 39.4767  (improvement: 27.3518)
% 
% ========== Summary (MaxN=6) ==========
%   N    Lip_init    Lip_final
%   3    66.8284      39.4767

% ========== N = 4 (MaxN = 6) ==========
% Building domain graph G_D...
% Domain graph: 16641 nodes, 65792 edges
% Building slit domain graph G_X...
% Grid nodes: 15900
% Graph edges: 60264
% Domain nodes: 16641  |  X_N nodes: 15900
% Building lane filling...
% Lane filling: boundary violations = 0
% Lane filling: nodes that moved    = 183 / 16641
% Lip evaluation nodes: 4269
% Caching image-node distance rows...
%   Done: 37.197 sec
% Computing initial Lip...
%   Done: 245.511 sec
% Initial Lip(c_bar) = 66.8284
% Running coordinate descent...
%   Iter   1: 2 nodes, Lip = 64.8284 (swaps: 2)
%   Iter   2: 6 nodes, Lip = 62.8284 (swaps: 6)
%   Iter   3: 10 nodes, Lip = 60.8284 (swaps: 10)
%   Iter   4: 14 nodes, Lip = 58.8284 (swaps: 14)
%   Iter   5: 18 nodes, Lip = 56.8284 (swaps: 18)
%   Iter   6: 22 nodes, Lip = 54.8284 (swaps: 22)
%   Iter   7: 26 nodes, Lip = 52.8284 (swaps: 26)
%   Iter   8: 30 nodes, Lip = 50.8284 (swaps: 30)
%   Iter   9: 34 nodes, Lip = 48.8284 (swaps: 34)
%   Iter  10: 38 nodes, Lip = 47.2548 (swaps: 38)
%   Iter  11: 4 nodes, Lip = 46.8284 (swaps: 2)
%   Iter  12: 38 nodes, Lip = 46.5477 (swaps: 38)
%   Iter  13: 4 nodes, Lip = 45.8406 (swaps: 3)
%   Iter  14: 2 nodes, Lip = 45.8284 (swaps: 2)
%   Iter  15: 2 nodes, Lip = 45.1335 (swaps: 2)
%   Iter  16: 2 nodes, Lip = 44.8284 (swaps: 1)
%   Iter  17: 42 nodes, Lip = 44.4264 (swaps: 42)
%   Iter  18: 4 nodes, Lip = 43.7193 (swaps: 2)
%   Iter  19: 4 nodes, Lip = 43.0122 (swaps: 2)
%   Iter  20: 4 nodes, Lip = 42.8284 (swaps: 4)
%   Iter  21: 46 nodes, Lip = 41.8284 (swaps: 46)
%   Iter  22: 4 nodes, Lip = 40.8909 (swaps: 2)
%   Iter  23: 4 nodes, Lip = 40.8284 (swaps: 2)
%   Iter  24: 50 nodes, Lip = 40.1838 (swaps: 50)
%   Iter  25: 4 nodes, Lip = 39.8284 (swaps: 3)
%   Iter  26: 2 nodes, Lip = 39.4767 (swaps: 1)
%   Iter  27: 4 nodes, Lip = 39.4767 (swaps: 0)
%   Converged at iteration 27.
% Final Lip(c_bar) = 39.4767  (improvement: 27.3518)
% 
% ========== Summary (MaxN=6) ==========
%   N    Lip_init    Lip_final
%   4    66.8284      39.4767