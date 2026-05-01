%% gd_top.m
% Coordinate descent to minimize Lip(c_bar_N) starting from lane filling.
% r_max = 3*h_D for all runs.
%
% Toggles:
%   record_video = true/false  — capture mp4 of descent
%   show_plots   = true/false  — show before/after plots and convergence
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
%   draw_colored_arrows.m
%   capture_frame.m

clear; clc;

%% ---- Parameters ----
N            = 2;
MaxN         = N + 2;
max_iters    = 100;
record_video = false;  % set true to capture mp4
show_plots   = false;  % set true to show mapping plots
%% -------------------

r   = 1./(1:(N+2));
tol = 1e-10;
h   = 1 / 2^(MaxN+1);  % grid spacing for Lip*h scaling

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
F_init = F;

%% Node sets
is_bnd  = (abs(pts_D(:,1)) < tol) | (abs(pts_D(:,1)-1) < tol) | ...
          (abs(pts_D(:,2)) < tol) | (abs(pts_D(:,2)-1) < tol);
all_idx = (1:size(pts_D,1))';
int_D   = find(~is_bnd);

ux_D  = unique(pts_D(:,1));
h_D   = min(diff(ux_D));
r_max = 3 * h_D;

[~, ord] = sort(pts_D(all_idx,1));
all_idx  = all_idx(ord);
n_all    = length(all_idx);

%% Build radius neighbor lists
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
E_X   = G_X.Edges.EndNodes;
adj_X = cell(size(pts_X,1), 1);
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
[lip_init, wu, wv] = global_lip(F, all_idx, pts_D, G_X, dist_cache, r_max);
fprintf('  Done: %.3f sec\n', toc);
fprintf('Initial Lip(c_bar_%d) = %.4f\n', N, lip_init);

%% Video setup
if record_video
    vid_file = sprintf('descent_lane_N%d.mp4', N);
    v = VideoWriter(vid_file, 'MPEG-4');
    v.FrameRate = 3;
    open(v);
    fprintf('Recording video: %s\n', vid_file);
    fig_vid = figure('Position', [100 100 700 650]);
    capture_frame(v, F, 0, lip_init, wu, wv, N, pts_D, pts_X, slits, is_bnd, G_X);
end

%% Coordinate descent
fprintf('\nRunning coordinate descent...\n');
lip_history = [lip_init];
lip_curr    = lip_init;

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

    if isempty(worst_nodes)
        fprintf('  No optimizable nodes at worst ratio — converged.\n');
        break
    end

    fprintf('  Iter %3d: %d nodes at worst ratio %.4f\n', ...
        iter, length(worst_nodes), lip_curr);

    n_swaps   = 0;
    idx_order = worst_nodes(randperm(length(worst_nodes)));

    for ii = 1:length(idx_order)
        u  = idx_order(ii);
        fu = F(u);
        nbrs_u = radius_nbrs{u};

        current_worst = 0;
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

    [lip_new, wu_new, wv_new] = global_lip(F, all_idx, pts_D, G_X, dist_cache, r_max);
    fprintf('           -> Lip = %.4f  (swaps: %d)\n', lip_new, n_swaps);

    if lip_new < lip_curr - 1e-9
        lip_curr = lip_new;
        wu = wu_new; wv = wv_new;
        lip_history(end+1) = lip_curr;
        if record_video
            capture_frame(v, F, iter, lip_curr, wu, wv, N, pts_D, pts_X, slits, is_bnd, G_X);
        end
    else
        fprintf('  Converged at iteration %d.\n', iter);
        break
    end
end

%% Close video
if record_video
    close(v);
    fprintf('Video saved: %s\n', vid_file);
end

%% Results
lip_h = lip_curr * h;
fprintf('\n========== Results (N=%d, MaxN=%d) ==========\n', N, MaxN);
fprintf('  Initial Lip:  %.4f\n', lip_init);
fprintf('  Final   Lip:  %.4f\n', lip_curr);
fprintf('  h:            %.6f\n', h);
fprintf('  Lip * h:      %.6f\n', lip_h);

%% Plots (optional)
if show_plots
    [~, wu_final, wv_final] = global_lip(F, all_idx, pts_D, G_X, dist_cache, r_max);

    for pass = 1:2
        if pass == 1
            F_plot = F_init; wu_p = wu; wv_p = wv;
            lip_p  = lip_init; label = 'initial';
        else
            F_plot = F; wu_p = wu_final; wv_p = wv_final;
            lip_p  = lip_curr; label = 'optimized';
        end

        int_D_ = find(~is_bnd);
        img_p  = pts_X(F_plot(int_D_),:);
        dx_p   = img_p(:,1) - pts_D(int_D_,1);
        dy_p   = img_p(:,2) - pts_D(int_D_,2);
        mov_p  = (abs(dx_p) + abs(dy_p)) > tol;

        figure('Position', [100 100 1200 500]);
        subplot(1,2,1); hold on; axis equal;
        xlim([0 1]); ylim([0 1]);
        title(sprintf('Domain Q — %s  N=%d  Lip=%.4f', label, N, lip_p));
        xlabel('x'); ylabel('y');
        plot(pts_D(int_D_,1), pts_D(int_D_,2), '.', 'Color', [0.88 0.88 0.88], 'MarkerSize', 3);
        plot(pts_D(wu_p,1), pts_D(wu_p,2), 'ro', 'MarkerSize', 14, 'LineWidth', 2.5);
        plot(pts_D(wv_p,1), pts_D(wv_p,2), 'rs', 'MarkerSize', 14, 'LineWidth', 2.5);
        hold off;

        subplot(1,2,2); hold on; axis equal;
        xlim([0 1]); ylim([0 1]);
        title(['$\bar{c}_' num2str(N) '$ — ' label], 'Interpreter', 'latex');
        xlabel('x'); ylabel('y');
        for k = 1:numel(slits)
            plot([slits(k).x slits(k).x], [slits(k).y0 slits(k).y1], 'k-', 'LineWidth', 2);
        end
        plot(pts_X(:,1), pts_X(:,2), '.', 'Color', [0.88 0.88 0.88], 'MarkerSize', 3);
        if any(mov_p)
            draw_colored_arrows(pts_D(int_D_,1), pts_D(int_D_,2), dx_p, dy_p, 'MaxMag', 1);
        end
        path_p = shortestpath(G_X, F_plot(wu_p), F_plot(wv_p));
        plot(pts_X(path_p,1), pts_X(path_p,2), 'g-', 'LineWidth', 2.5);
        plot(pts_X(F_plot(wu_p),1), pts_X(F_plot(wu_p),2), 'ro', 'MarkerSize', 14, 'LineWidth', 2.5);
        plot(pts_X(F_plot(wv_p),1), pts_X(F_plot(wv_p),2), 'rs', 'MarkerSize', 14, 'LineWidth', 2.5);
        hold off;
    end

    figure;
    plot(0:length(lip_history)-1, lip_history, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Iteration');
    ylabel('$\mathrm{Lip}(\bar{c}_N)$', 'Interpreter', 'latex');
    title(['$\bar{c}_' num2str(N) '$ convergence (MaxN=' num2str(MaxN) ')'], 'Interpreter', 'latex');
    grid on;
end