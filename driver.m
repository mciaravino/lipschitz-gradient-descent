%% driver.m
% Runs gradient descent n_runs times for each N.
% Uses MaxN = N+2 for each run.
% Records Lip_final and Lip*h for each run.
% Prints summary table and plots results.
%
% Requires all files in gd_top.m's requirements list.

clear; clc;

%% ---- Parameters ----
N_values  = [2,3,4];  % add 5 if desired
n_runs    = 5;
max_iters = 100;
%% -------------------

tol     = 1e-10;
results = [];  % [N, MaxN, h, run, lip_final, lip_h]

for ni = 1:length(N_values)
    N    = N_values(ni);
    MaxN = 6;
    h    = 1 / 2^(MaxN+1);
    r    = 1./(1:(N+2));

    fprintf('\n========== N=%d, MaxN=%d, h=%.6f ==========\n', N, MaxN, h);

    %% Build geometry once per N
    slits = hl_draw_slits_iter(r, N, 'DrawGrid', false);
    close all;

    [~, pts_D] = hl_build_domain_graph(MaxN);
    pts_X = hl_overlay_grid_with_slit_holes(slits, N, 'MaxN', MaxN);
    G_X   = hl_build_graph_from_points(pts_X, slits);
    close all;

    fprintf('Domain nodes: %d  |  X_N nodes: %d\n', size(pts_D,1), size(pts_X,1));

    %% Node setup
    is_bnd  = (abs(pts_D(:,1)) < tol) | (abs(pts_D(:,1)-1) < tol) | ...
              (abs(pts_D(:,2)) < tol) | (abs(pts_D(:,2)-1) < tol);
    all_idx = (1:size(pts_D,1))';
    [~, ord] = sort(pts_D(all_idx,1));
    all_idx  = all_idx(ord);
    n_all    = length(all_idx);

    ux_D  = unique(pts_D(:,1));
    h_D   = min(diff(ux_D));
    r_max = 3 * h_D;

    %% Radius neighbor lists (built once per N)
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

    %% Adjacency list for G_X (built once per N)
    E_X   = G_X.Edges.EndNodes;
    adj_X = cell(size(pts_X,1), 1);
    for k = 1:size(E_X,1)
        u = E_X(k,1); v = E_X(k,2);
        adj_X{u} = [adj_X{u}, v];
        adj_X{v} = [adj_X{v}, u];
    end

    %% Multiple runs
    for run = 1:n_runs
        fprintf('\n  -- Run %d/%d --\n', run, n_runs);

        %% Fresh lane filling each run
        [F, ~] = hl_lane_filling(pts_D, slits, G_X, pts_X);

        %% Fresh distance cache each run
        dist_cache = containers.Map('KeyType','int32','ValueType','any');
        for nd = unique(F)'
            gd_get_row(nd, G_X, dist_cache);
        end

        %% Initial Lip
        [lip_curr, ~, ~] = global_lip(F, all_idx, pts_D, G_X, dist_cache, r_max);
        fprintf('  Initial Lip = %.4f\n', lip_curr);

        %% Descent
        for iter = 1:max_iters
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
            if isempty(worst_nodes), break; end

            idx_order = worst_nodes(randperm(length(worst_nodes)));
            for ii = 1:length(idx_order)
                u      = idx_order(ii);
                fu     = F(u);
                nbrs_u = radius_nbrs{u};

                current_worst = 0;
                for bi = 1:length(nbrs_u)
                    vb    = nbrs_u(bi);
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
                        vb    = nbrs_u(bi);
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
            end

            [lip_new, ~, ~] = global_lip(F, all_idx, pts_D, G_X, dist_cache, r_max);
            if lip_new < lip_curr - 1e-9
                lip_curr = lip_new;
            else
                break
            end
        end

        lip_h = lip_curr * h;
        fprintf('  Final Lip = %.4f,  Lip*h = %.6f\n', lip_curr, lip_h);
        results(end+1,:) = [N, MaxN, h, run, lip_curr, lip_h];
    end
end

%% Summary table
fprintf('\n========== Summary ==========\n');
fprintf('  N   MaxN      h         run   Lip_final   Lip*h\n');
for i = 1:size(results,1)
    fprintf('  %d    %d    %.6f    %d     %.4f    %.6f\n', ...
        results(i,1), results(i,2), results(i,3), ...
        results(i,4), results(i,5), results(i,6));
end

fprintf('\n  N   MaxN   mean(Lip*h)   std(Lip*h)   min(Lip*h)   max(Lip*h)\n');
for ni = 1:length(N_values)
    N    = N_values(ni);
    MaxN = N + 2;
    idx  = results(:,1) == N;
    vals = results(idx,6);
    fprintf('  %d    %d      %.6f      %.6f      %.6f      %.6f\n', ...
        N, MaxN, mean(vals), std(vals), min(vals), max(vals));
end

%% Save
save('driver_results.mat', 'results', 'N_values');
fprintf('\nResults saved to driver_results.mat\n');

%% Plot
figure; hold on;
colors = {'b','r','g','m'};
labels = {};
for ni = 1:length(N_values)
    N   = N_values(ni);
    idx = results(:,1) == N;
    vals = results(idx,6);
    scatter(N*ones(sum(idx),1), vals, 60, colors{ni}, 'filled', 'MarkerFaceAlpha', 0.5);
    plot(N, mean(vals), 's', 'Color', colors{ni}, 'MarkerSize', 14, ...
        'LineWidth', 2.5, 'MarkerFaceColor', colors{ni});
    labels{end+1} = sprintf('N=%d runs', N);
    labels{end+1} = sprintf('N=%d mean', N);
end
xlabel('N');
ylabel('$\mathrm{Lip}(\bar{c}_N) \times h$', 'Interpreter', 'latex');
title('$\mathrm{Lip}(\bar{c}_N) \times h$ vs $N$', 'Interpreter', 'latex');
legend(labels, 'Location', 'northeast');
grid on;
ylim([0 max(results(:,6))*1.3]);
hold off;