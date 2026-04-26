%% compute_lip_lane.m
% Computes the Lipschitz constant of the lane filling map.
% Checks both graph edges and cross-slit pairs.
%
% Requires in same directory:
%   hl_draw_slits_iter.m
%   hl_overlay_grid_with_slit_holes.m
%   hl_build_graph_from_points.m
%   hl_lane_filling.m
%   gd_get_row.m
%   gd_dlookup.m

clear; clc;

%% Parameters
N    = 2;
MaxN = 6;
r    = 1./(1:(N+2));

%% Build domain
slits = hl_draw_slits_iter(r, N, 'DrawGrid', false);
pts   = hl_overlay_grid_with_slit_holes(slits, N, 'MaxN', MaxN);
G     = hl_build_graph_from_points(pts, slits);
n     = size(pts, 1);
close all;

%% Build lane filling
fprintf('Building lane filling...\n');
tic
[F, ~] = hl_lane_filling(pts, slits, G);
fprintf('  Done: %.3f sec\n', toc);

%% Edge list
E = G.Edges.EndNodes;

%% Build cross-slit pairs
ux = unique(pts(:,1));
h  = min(diff(ux));

CS = [];
for k = 1:numel(slits)
    sx  = slits(k).x;
    sy0 = slits(k).y0;
    sy1 = slits(k).y1;

    left_idx  = find(abs(pts(:,1) - (sx - h)) < h*0.6 & ...
                     pts(:,2) >= sy0 - h & pts(:,2) <= sy1 + h);
    right_idx = find(abs(pts(:,1) - (sx + h)) < h*0.6 & ...
                     pts(:,2) >= sy0 - h & pts(:,2) <= sy1 + h);

    for li = 1:length(left_idx)
        for ri = 1:length(right_idx)
            if abs(pts(left_idx(li),2) - pts(right_idx(ri),2)) < h*0.6
                CS(end+1,:) = [left_idx(li), right_idx(ri)];
            end
        end
    end
end
fprintf('Cross-slit pairs: %d\n', size(CS,1));

%% Distance cache
dist_cache = containers.Map('KeyType','int32','ValueType','any');

% Cache rows for all image nodes
fprintf('Caching image-node distance rows...\n');
tic
for nd = unique(F)'
    gd_get_row(nd, G, dist_cache);
end
fprintf('  Done: %.3f sec\n', toc);

%% Compute Lipschitz constant over edges
lip = 0; wu = -1; wv = -1;

for k = 1:size(E,1)
    u = E(k,1); v = E(k,2);
    d_dom = norm(pts(u,:) - pts(v,:));
    if d_dom < 1e-12, continue; end
    d_img = gd_dlookup(F(u), F(v), G, dist_cache);
    if isinf(d_img), continue; end
    ratio = d_img / d_dom;
    if ratio > lip
        lip = ratio; wu = u; wv = v;
    end
end

% Check cross-slit pairs
for k = 1:size(CS,1)
    u = CS(k,1); v = CS(k,2);
    d_dom = norm(pts(u,:) - pts(v,:));
    if d_dom < 1e-12, continue; end
    d_img = gd_dlookup(F(u), F(v), G, dist_cache);
    if isinf(d_img), continue; end
    ratio = d_img / d_dom;
    if ratio > lip
        lip = ratio; wu = u; wv = v;
    end
end

%% Report
fprintf('\n========== Results (N=%d) ==========\n', N);
fprintf('  Lip(c_bar) = %.4f\n', lip);
fprintf('  Worst pair: node %d (%.4f,%.4f) and node %d (%.4f,%.4f)\n', ...
    wu, pts(wu,1), pts(wu,2), wv, pts(wv,1), pts(wv,2));
fprintf('  domain distance:  %.6f\n', norm(pts(wu,:) - pts(wv,:)));
fprintf('  intrinsic dist:   %.6f\n', gd_dlookup(F(wu), F(wv), G, dist_cache));