%% test_lane_filling.m
% Visualizes the lane filling map with optional Lipschitz computation.
%
% TOGGLES (edit the Parameters block below):
%   compute_lip  true  = compute Lip, circle worst pair, trace shortest path
%                false = visualization only, fast for any N
%   exhaustive   true  = check ALL n*(n-1)/2 domain pairs  [WARNING: MaxN<=4]
%                false = check graph edges + cross-slit pairs only
%   slit_cols    true  = restrict domain to 3 columns near each slit
%                false = use full domain grid
%
% Requires:
%   hl_draw_slits_iter.m
%   hl_overlay_grid_with_slit_holes.m
%   hl_build_domain_graph.m
%   hl_build_graph_from_points.m
%   hl_lane_filling.m
%   gd_get_row.m
%   gd_dlookup.m

clear; clc;

%% ---- Parameters (edit these) ----
N            = 2;
MaxN_min     = N + 2;   % minimum to resolve X_N (do not go below this)
MaxN_max     = 6;       % memory cap (beyond this, distance matrix too large)
MaxN         = 4 %min(max(N + 2, MaxN_min), MaxN_max);  % auto, override manually
compute_lip  = true;
exhaustive   = true;   % WARNING: only use when MaxN <= 4
slit_cols    = false;   % restrict domain to slit-adjacent columns
%% ---------------------------------

r = 1./(1:(N+2));
tol = 1e-10;

%% Build slit geometry
slits = hl_draw_slits_iter(r, N, 'DrawGrid', false);
close all;

%% Build domain graph G_D (full uniform grid)
fprintf('Building domain graph G_D...\n');
[G_D, pts_D_full] = hl_build_domain_graph(MaxN);

%% Build slit domain graph G_X
fprintf('Building slit domain graph G_X...\n');
pts_X = hl_overlay_grid_with_slit_holes(slits, N, 'MaxN', MaxN);
G_X   = hl_build_graph_from_points(pts_X, slits);
close all;
fprintf('Domain nodes: %d  |  X_N nodes: %d\n', size(pts_D_full,1), size(pts_X,1));

% %% Build lane filling on full domain
% fprintf('Building lane filling...\n');
% tic
% [F_full, ~] = hl_lane_filling(pts_D_full, slits, G_X, pts_X);
% fprintf('  Done: %.3f sec\n', toc);

%% Build top filling on full domain
fprintf('Building top filling...\n');
tic
[F_full, ~] = hl_top_filling(pts_D_full, slits, G_X, pts_X, slit_cols);
fprintf('  Done: %.3f sec\n', toc);


%% Optionally restrict domain to slit-adjacent columns
% Keeps points within x = slit_x +/- h and y within slit y-range +/- h.
% All peg holes in X_N are kept — we just use fewer pegs in D.
if slit_cols
    ux_D = unique(pts_D_full(:,1));
    h_D  = min(diff(ux_D));
    keep = false(size(pts_D_full,1), 1);
    for k = 1:numel(slits)
        sx  = slits(k).x;
        sy0 = slits(k).y0;
        sy1 = slits(k).y1;
        in_x = abs(pts_D_full(:,1) - sx) <= h_D * 1.5;
        in_y = pts_D_full(:,2) >= sy0 - h_D & pts_D_full(:,2) <= sy1 + h_D;
        keep = keep | (in_x & in_y);
    end
    keep_idx = find(keep);
    pts_D = pts_D_full(keep_idx, :);
    F     = F_full(keep_idx);
    fprintf('Slit-column restriction: %d / %d domain nodes kept\n', ...
        length(keep_idx), size(pts_D_full,1));
else
    pts_D    = pts_D_full;
    F        = F_full;
    keep_idx = (1:size(pts_D_full,1))';
end

%% Interior nodes of (possibly restricted) domain
is_bnd = (abs(pts_D(:,1)) < tol) | (abs(pts_D(:,1)-1) < tol) | ...
         (abs(pts_D(:,2)) < tol) | (abs(pts_D(:,2)-1) < tol);
interior_D = find(~is_bnd);

%% Displacement vectors for visualization
img_pts = pts_X(F(interior_D), :);
dx = img_pts(:,1) - pts_D(interior_D,1);
dy = img_pts(:,2) - pts_D(interior_D,2);
moved = (abs(dx) + abs(dy)) > tol;

%% Lipschitz computation
if compute_lip

    %% Cross-slit pairs (indices into pts_X)
    ux = unique(pts_X(:,1));
    h  = min(diff(ux));
    CS = [];
    for k = 1:numel(slits)
        sx  = slits(k).x;
        sy0 = slits(k).y0;
        sy1 = slits(k).y1;
        left_idx  = find(abs(pts_X(:,1) - (sx-h)) < h*0.6 & ...
                         pts_X(:,2) >= sy0-h & pts_X(:,2) <= sy1+h);
        right_idx = find(abs(pts_X(:,1) - (sx+h)) < h*0.6 & ...
                         pts_X(:,2) >= sy0-h & pts_X(:,2) <= sy1+h);
        for li = 1:length(left_idx)
            for ri = 1:length(right_idx)
                if abs(pts_X(left_idx(li),2) - pts_X(right_idx(ri),2)) < h*0.6
                    CS(end+1,:) = [left_idx(li), right_idx(ri)];
                end
            end
        end
    end
    fprintf('Cross-slit pairs: %d\n', size(CS,1));

    %% Cache image node distances
    dist_cache = containers.Map('KeyType','int32','ValueType','any');
    fprintf('Caching image-node distance rows...\n');
    tic
    for nd = unique(F)'
        gd_get_row(nd, G_X, dist_cache);
    end
    fprintf('  Done: %.3f sec\n', toc);

    %% Compute Lip
    n_pts = size(pts_D, 1);
    lip = 0; wu = -1; wv = -1;
    d_dom_worst = 0; d_img_worst = 0;

    if exhaustive
        fprintf('Exhaustive mode: checking %d pairs...\n', n_pts*(n_pts-1)/2);
        tic
        for ui = 1:n_pts
            for vi = ui+1:n_pts
                d_dom = norm(pts_D(ui,:) - pts_D(vi,:));
                if d_dom < 1e-12, continue; end
                d_img = gd_dlookup(F(ui), F(vi), G_X, dist_cache);
                if isinf(d_img), continue; end
                ratio = d_img / d_dom;
                if ratio > lip
                    lip = ratio; wu = ui; wv = vi;
                    d_dom_worst = d_dom; d_img_worst = d_img;
                end
            end
        end
        fprintf('  Done: %.3f sec\n', toc);
    else
        % Graph edges — remap E_D to restricted pts_D indices
        E_D = G_D.Edges.EndNodes;
        % Build reverse map: full index -> restricted index
        rev = zeros(size(pts_D_full,1), 1);
        rev(keep_idx) = 1:length(keep_idx);

        for k = 1:size(E_D,1)
            ui_full = E_D(k,1); vi_full = E_D(k,2);
            ui = rev(ui_full); vi = rev(vi_full);
            if ui == 0 || vi == 0, continue; end  % not in restricted set
            d_dom = norm(pts_D(ui,:) - pts_D(vi,:));
            if d_dom < 1e-12, continue; end
            d_img = gd_dlookup(F(ui), F(vi), G_X, dist_cache);
            if isinf(d_img), continue; end
            ratio = d_img / d_dom;
            if ratio > lip
                lip = ratio; wu = ui; wv = vi;
                d_dom_worst = d_dom; d_img_worst = d_img;
            end
        end

        % Cross-slit pairs: find restricted domain nodes mapping to CS image pairs
        for k = 1:size(CS,1)
            xi = CS(k,1); xj = CS(k,2);
            du = find(F == xi);
            dv = find(F == xj);
            for a = 1:length(du)
                for b = 1:length(dv)
                    ui = du(a); vi = dv(b);
                    d_dom = norm(pts_D(ui,:) - pts_D(vi,:));
                    if d_dom < 1e-12, continue; end
                    d_img = gd_dlookup(F(ui), F(vi), G_X, dist_cache);
                    if isinf(d_img), continue; end
                    ratio = d_img / d_dom;
                    if ratio > lip
                        lip = ratio; wu = ui; wv = vi;
                        d_dom_worst = d_dom; d_img_worst = d_img;
                    end
                end
            end
        end
    end

    %% Report
    fprintf('\n========== Results (N=%d) ==========\n', N);
    fprintf('  Worst domain pair:\n');
    fprintf('    u = (%.4f, %.4f)\n', pts_D(wu,1), pts_D(wu,2));
    fprintf('    v = (%.4f, %.4f)\n', pts_D(wv,1), pts_D(wv,2));
    fprintf('  d_D(u,v)             = %.6f  [Euclidean in domain]\n', d_dom_worst);
    fprintf('  d_X(cbar(u),cbar(v)) = %.6f  [intrinsic in X_N]\n',   d_img_worst);
    fprintf('  Lip(c_bar)           = %.4f\n', lip);

    %% Shortest path in X_N between image nodes of worst pair
    path_nodes  = shortestpath(G_X, F(wu), F(wv));
    path_coords = pts_X(path_nodes, :);

end

%% ---- Plots ----
figure('Position', [100 100 1200 500]);

% Left: domain — gray empty pegs, blue active pegs
subplot(1,2,1); hold on; axis equal;
xlim([0 1]); ylim([0 1]);
if compute_lip
    title(sprintf('Domain Q  (N=%d,  Lip=%.4f)', N, lip));
else
    title(sprintf('Domain Q  (N=%d)', N));
end
xlabel('x'); ylabel('y');

% All domain grid points in gray
is_bnd_full = (abs(pts_D_full(:,1)) < tol) | (abs(pts_D_full(:,1)-1) < tol) | ...
              (abs(pts_D_full(:,2)) < tol) | (abs(pts_D_full(:,2)-1) < tol);
interior_full = find(~is_bnd_full);
plot(pts_D_full(interior_full,1), pts_D_full(interior_full,2), '.', ...
    'Color', [0.88 0.88 0.88], 'MarkerSize', 3);

% Active (restricted) interior nodes in blue
plot(pts_D(interior_D,1), pts_D(interior_D,2), 'b.', 'MarkerSize', 5);

if compute_lip
    plot(pts_D(wu,1), pts_D(wu,2), 'ro', 'MarkerSize', 14, 'LineWidth', 2.5);
    plot(pts_D(wv,1), pts_D(wv,2), 'rs', 'MarkerSize', 14, 'LineWidth', 2.5);
    legend({'Empty pegs','Active pegs','Worst u','Worst v'}, 'Location','northeast');
else
    legend({'Empty pegs','Active pegs'}, 'Location','northeast');
end
hold off;

% Right: image in X_N
subplot(1,2,2); hold on; axis equal;
xlim([0 1]); ylim([0 1]);
title(sprintf('Image c-bar(u) in X_N  (N=%d)', N));
xlabel('x'); ylabel('y');
for k = 1:numel(slits)
    plot([slits(k).x slits(k).x], [slits(k).y0 slits(k).y1], ...
        'k-', 'LineWidth', 2);
end
plot(pts_X(:,1), pts_X(:,2), '.', 'Color', [0.88 0.88 0.88], 'MarkerSize', 3);
if any(moved)
    quiver(pts_D(interior_D(moved),1), pts_D(interior_D(moved),2), ...
           dx(moved), dy(moved), 0, ...
           'Color', [0.2 0.5 0.9], 'LineWidth', 1.2, 'MaxHeadSize', 0.5);
end
if compute_lip
    plot(path_coords(:,1), path_coords(:,2), 'g-', 'LineWidth', 2.5);
    plot(pts_X(F(wu),1), pts_X(F(wu),2), 'ro', 'MarkerSize', 14, 'LineWidth', 2.5);
    plot(pts_X(F(wv),1), pts_X(F(wv),2), 'rs', 'MarkerSize', 14, 'LineWidth', 2.5);
    legend({'Slits','Grid','Displacement','Shortest path','F(u)','F(v)'}, ...
        'Location','northeast');
end
hold off;