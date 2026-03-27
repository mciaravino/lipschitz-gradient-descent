%% drivingScript.m
% Coordinate descent optimization of Lip(F) for N = 2..4
% in a Hakobyan-Li slit domain.
%
% Key fix: Lipschitz computation now includes cross-slit pairs --
% nodes that straddle a slit horizontally. These are NOT graph edges
% (correctly excluded from the graph) but are close in Euclidean
% distance and far in intrinsic distance, making them the critical
% pairs for detecting high Lipschitz ratios.
%
% Requires in same directory:
%   hl_draw_slits_iter.m
%   hl_overlay_grid_with_slit_holes.m
%   hl_build_graph_from_points.m
%   gd_get_row.m
%   gd_dlookup.m
%   gd_compute_lip.m
%   gd_local_lip.m

clear; clc;

%% Parameters
N_values  = 2:4;  % Run a single N value by typing N:N
max_iters = 100;
MaxN      = 5;    % Grid resolution — increase for finer grid, decrease for speed

%% Storage
results = zeros(length(N_values), 3);  % [N, lip_init, lip_final]

%% Main loop
for ni = 1:length(N_values)

    N = N_values(ni);
    fprintf('\n========== N = %d ==========\n', N);
    r = 1./(1:(N+2));

    %% Build domain
    slits = hl_draw_slits_iter(r, N, 'DrawGrid', false);
    pts   = hl_overlay_grid_with_slit_holes(slits, N, 'MaxN', MaxN);
    G     = hl_build_graph_from_points(pts, slits);
    n     = size(pts, 1);

    %% Boundary and interior nodes
    tol = 1e-10;
    is_boundary = (abs(pts(:,1)) < tol) | (abs(pts(:,1)-1) < tol) | ...
                  (abs(pts(:,2)) < tol) | (abs(pts(:,2)-1) < tol);
    boundary_idx = find(is_boundary);
    interior_idx = find(~is_boundary);
    fprintf('  Boundary: %d  Interior: %d\n', ...
        length(boundary_idx), length(interior_idx));

    top_bnd    = boundary_idx(abs(pts(boundary_idx,2) - 1) < tol);
    bottom_bnd = boundary_idx(abs(pts(boundary_idx,2))     < tol);

    %% Grid spacing
    dx = min(diff(unique(pts(:,1))));
    dy = min(diff(unique(pts(:,2))));
    h  = min(dx, dy);

    %% Precompute edge list and adjacency list
    E   = G.Edges.EndNodes;
    adj = cell(n, 1);
    for k = 1:size(E, 1)
        u = E(k,1); v = E(k,2);
        adj{u} = [adj{u}, v];
        adj{v} = [adj{v}, u];
    end

    %% Build cross-slit pairs
    % Nodes that straddle a slit horizontally: close in Euclidean distance
    % but separated by the slit obstacle, so intrinsically far apart.
    % These pairs are deliberately excluded from the graph but must be
    % included in the Lipschitz computation.
    fprintf('Building cross-slit pairs...\n');
    CS      = zeros(0, 2);
    cs_tol  = h * 1.5;
    cs_adj  = cell(n, 1);  % cross-slit neighbors per node

    for k = 1:numel(slits)
        sx  = slits(k).x;
        sy0 = slits(k).y0;
        sy1 = slits(k).y1;

        left_idx  = find(abs(pts(:,1) - (sx - h)) < h*0.6 & ...
                         pts(:,2) >= sy0 - cs_tol & pts(:,2) <= sy1 + cs_tol);
        right_idx = find(abs(pts(:,1) - (sx + h)) < h*0.6 & ...
                         pts(:,2) >= sy0 - cs_tol & pts(:,2) <= sy1 + cs_tol);

        for li = 1:length(left_idx)
            for ri = 1:length(right_idx)
                if abs(pts(left_idx(li),2) - pts(right_idx(ri),2)) < h*0.6
                    u = left_idx(li);
                    v = right_idx(ri);
                    CS(end+1,:) = [u, v];
                    cs_adj{u}   = [cs_adj{u}, v];
                    cs_adj{v}   = [cs_adj{v}, u];
                end
            end
        end
    end

    fprintf('  Cross-slit pairs: %d\n', size(CS,1));

    %% Distance cache
    dist_cache = containers.Map('KeyType','int32','ValueType','any');

    %% Build vertical filling (initial guess)
    fprintf('Building vertical filling...\n');
    tic

    F = zeros(n, 1);
    for i = 1:length(boundary_idx)
        F(boundary_idx(i)) = boundary_idx(i);
    end

    bt_pairs = zeros(length(interior_idx), 2);
    for i = 1:length(interior_idx)
        u  = interior_idx(i);
        xu = pts(u,1);
        [~, ti] = min(abs(pts(top_bnd,1)    - xu));
        [~, bi] = min(abs(pts(bottom_bnd,1) - xu));
        bt_pairs(i,:) = [bottom_bnd(bi), top_bnd(ti)];
    end

    unique_bnodes = unique(bt_pairs(:,1));
    fprintf('  Caching %d bottom-node distance rows...\n', length(unique_bnodes));
    for k = 1:length(unique_bnodes)
        gd_get_row(unique_bnodes(k), G, dist_cache);
    end

    path_cache = containers.Map('KeyType','char','ValueType','any');

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
            cum_dist = cum_dist + ...
                norm(pts(path_nodes(k),:) - pts(path_nodes(k-1),:));
            if cum_dist >= target_dist
                chosen = path_nodes(k);
                break
            end
        end
        F(u) = chosen;
    end

    fprintf('  Unique image nodes: %d\n', length(unique(F)));
    fprintf('  Done: %.3f sec\n', toc);

    %% Cache distance rows for all initial image nodes
    fprintf('Caching image-node distance rows...\n');
    tic
    for nd = unique(F)'
        gd_get_row(nd, G, dist_cache);
    end
    fprintf('  Done: %.3f sec  (cache size: %d rows)\n', toc, dist_cache.Count);

    %% Compute initial Lip (includes cross-slit pairs)
    [lip_init, ~, ~] = gd_compute_lip(F, E, CS, pts, G, dist_cache);
    fprintf('Initial Lip(F) = %.4f\n', lip_init);

    %% Coordinate descent
    fprintf('Running coordinate descent...\n');
    lip_history    = zeros(max_iters+1, 1);
    lip_history(1) = lip_init;
    lip_curr       = lip_init;

    for iter = 1:max_iters

        n_swaps   = 0;
        idx_order = interior_idx(randperm(length(interior_idx)));

        for ii = 1:length(idx_order)
            u    = idx_order(ii);
            fu   = F(u);

            % Local cost includes both graph neighbors and cross-slit partners
            current_worst = gd_local_lip(u, F, adj, cs_adj{u}, pts, G, dist_cache);

            candidates = [fu, adj{fu}];
            best_cost  = current_worst;
            best_node  = fu;

            for c = 1:length(candidates)
                cand = candidates(c);
                if is_boundary(cand), continue; end

                gd_get_row(cand, G, dist_cache);
                F(u) = cand;
                tentative_worst = gd_local_lip(u, F, adj, cs_adj{u}, pts, G, dist_cache);

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

        % Recompute global Lip after full sweep
        [lip_curr, ~, ~] = gd_compute_lip(F, E, CS, pts, G, dist_cache);
        lip_history(iter+1) = lip_curr;
        fprintf('  Iter %3d:  Lip = %.4f  (swaps: %d)\n', iter, lip_curr, n_swaps);

        if n_swaps == 0
            fprintf('  Converged at iteration %d.\n', iter);
            lip_history = lip_history(1:iter+1);
            break
        end
    end

    fprintf('Final Lip(F) = %.4f  (improvement: %.4f)\n', ...
        lip_curr, lip_init - lip_curr);

    results(ni,:) = [N, lip_init, lip_curr];

    %% Per-N convergence plot
    figure;
    plot(0:length(lip_history)-1, lip_history, 'bo-', ...
        'LineWidth', 2, 'MarkerSize', 6);
    xlabel('Iteration');
    ylabel('Lip(F)');
    title(sprintf('Coordinate descent convergence (N=%d)', N));
    grid on;

    %% Displacement plot
    img_pts = pts(F(interior_idx), :);
    dx = img_pts(:,1) - pts(interior_idx,1);
    dy = img_pts(:,2) - pts(interior_idx,2);
    moved = (abs(dx) + abs(dy)) > 1e-10;

    fprintf('Nodes with arrows: %d / %d\n', sum(moved), length(interior_idx));

    figure; hold on; axis equal;
    xlim([0 1]); ylim([0 1]);
    title(sprintf('Displacement of F* (N=%d, Lip=%.4f)', N, lip_curr));
    xlabel('x'); ylabel('y');

    for k = 1:numel(slits)
        plot([slits(k).x slits(k).x], [slits(k).y0 slits(k).y1], ...
            'k-', 'LineWidth', 2);
    end
    plot(pts(interior_idx,1), pts(interior_idx,2), '.', ...
        'Color', [0.85 0.85 0.85], 'MarkerSize', 3);
    quiver(pts(interior_idx(moved),1), pts(interior_idx(moved),2), ...
           dx(moved), dy(moved), 0, ...
           'Color', [0.2 0.5 0.9], 'LineWidth', 1.2, 'MaxHeadSize', 0.5);
    hold off;

end

%% Summary table
fprintf('\n========== Summary ==========\n');
fprintf('  N    Lip_init    Lip_final\n');
for ni = 1:size(results,1)
    fprintf('  %d    %.4f      %.4f\n', ...
        results(ni,1), results(ni,2), results(ni,3));
end

%% Summary plot
figure; hold on;
plot(results(:,1), results(:,2), 'rs--', 'LineWidth', 2, 'MarkerSize', 8);
plot(results(:,1), results(:,3), 'bo-',  'LineWidth', 2, 'MarkerSize', 8);
xlabel('N (dyadic depth)');
ylabel('Lip(F)');
title('Lipschitz constant vs dyadic depth: initial and optimized');
legend({'Vertical filling (initial)', 'After coordinate descent'}, ...
    'Location', 'northwest');
grid on;
hold off;