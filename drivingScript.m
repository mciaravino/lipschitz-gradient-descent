%% march4Script_gd_loop.m
% Coordinate descent optimization of Lip(F) for N = 2..5
% in a Hakobyan-Li slit domain.
%
% Uses sparse distance computation to avoid memory overflow:
%   - gd_get_row(src, G, cache) computes and caches one distance row
%   - gd_dlookup(i, j, G, cache) retrieves distance using the cache
%   - No full n x n distance matrix is ever stored
%
% Requires in same directory:
%   hl_draw_slits_iter.m
%   hl_overlay_grid_with_slit_holes.m
%   hl_build_graph_from_points.m
%   gd_get_row.m
%   gd_dlookup.m

clear; clc;

%% Parameters
N_values  = 2:4; %Run a single N value by typing N:N
max_iters = 100;

%% Storage
results = zeros(length(N_values), 3);  % [N, lip_init, lip_final]

%% Main loop
for ni = 1:length(N_values)

    N = N_values(ni);
    fprintf('\n========== N = %d ==========\n', N);
    r = 1./(1:(N+2));

    %% Build domain
    slits = hl_draw_slits_iter(r, N, 'DrawGrid', false);
    pts   = hl_overlay_grid_with_slit_holes(slits, N, 'MaxN', 6);
    G     = hl_build_graph_from_points(pts);
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

    %% Precompute edge list and adjacency list
    E   = G.Edges.EndNodes;
    adj = cell(n, 1);
    for k = 1:size(E, 1)
        u = E(k,1); v = E(k,2);
        adj{u} = [adj{u}, v];
        adj{v} = [adj{v}, u];
    end

    %% Distance cache
    % Stores one distance row per source node, computed on demand.
    % Never stores the full n x n matrix.
    dist_cache = containers.Map('KeyType','int32','ValueType','any');

    %% Build vertical filling (initial guess)
    fprintf('Building vertical filling...\n');
    tic

    F = zeros(n, 1);
    for i = 1:length(boundary_idx)
        F(boundary_idx(i)) = boundary_idx(i);
    end

    % Find (b_node, t_node) pair for each interior node
    bt_pairs = zeros(length(interior_idx), 2);
    for i = 1:length(interior_idx)
        u  = interior_idx(i);
        xu = pts(u,1);
        [~, ti] = min(abs(pts(top_bnd,1)    - xu));
        [~, bi] = min(abs(pts(bottom_bnd,1) - xu));
        bt_pairs(i,:) = [bottom_bnd(bi), top_bnd(ti)];
    end

    % Cache distance rows for all unique bottom nodes
    unique_bnodes = unique(bt_pairs(:,1));
    fprintf('  Caching %d bottom-node distance rows...\n', length(unique_bnodes));
    for k = 1:length(unique_bnodes)
        gd_get_row(unique_bnodes(k), G, dist_cache);
    end

    % Compute shortest paths, caching per unique (b,t) pair
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

    %% Compute initial Lip
    lip_init = 0; wu_init = -1; wv_init = -1;
    for k = 1:size(E,1)
        u = E(k,1); v = E(k,2);
        d_dom = norm(pts(u,:) - pts(v,:));
        if d_dom < 1e-12, continue; end
        d_img = gd_dlookup(F(u), F(v), G, dist_cache);
        if isinf(d_img), continue; end
        ratio = d_img / d_dom;
        if ratio > lip_init
            lip_init = ratio; wu_init = u; wv_init = v;
        end
    end
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
            nbrs = adj{u};

            % Current local worst ratio for edges incident to u
            current_worst = 0;
            for j = 1:length(nbrs)
                v     = nbrs(j);
                d_dom = norm(pts(u,:) - pts(v,:));
                if d_dom < 1e-12, continue; end
                d_img = gd_dlookup(F(u), F(v), G, dist_cache);
                if isinf(d_img), continue; end
                current_worst = max(current_worst, d_img / d_dom);
            end

            % Try current image node and each of its graph neighbors
            candidates = [fu, adj{fu}];
            best_cost  = current_worst;
            best_node  = fu;

            for c = 1:length(candidates)
                cand = candidates(c);
                if is_boundary(cand), continue; end

                % Ensure candidate distance row is cached
                gd_get_row(cand, G, dist_cache);

                tentative_worst = 0;
                for j = 1:length(nbrs)
                    v     = nbrs(j);
                    d_dom = norm(pts(u,:) - pts(v,:));
                    if d_dom < 1e-12, continue; end
                    d_img = gd_dlookup(cand, F(v), G, dist_cache);
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

        % Recompute global Lip after full sweep
        lip_curr = 0;
        for k = 1:size(E,1)
            u = E(k,1); v = E(k,2);
            d_dom = norm(pts(u,:) - pts(v,:));
            if d_dom < 1e-12, continue; end
            d_img = gd_dlookup(F(u), F(v), G, dist_cache);
            if isinf(d_img), continue; end
            lip_curr = max(lip_curr, d_img / d_dom);
        end

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