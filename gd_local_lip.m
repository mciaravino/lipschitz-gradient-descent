function worst = gd_local_lip(u, F_map, adj_cell, CS_u, pts_mat, G_ref, cache)
% GD_LOCAL_LIP  Worst Lipschitz ratio over edges and cross-slit pairs incident to u.
%
%   worst = gd_local_lip(u, F_map, adj_cell, CS_u, pts_mat, G_ref, cache)
%
%   Checks both graph neighbors of u (from adj_cell) and cross-slit
%   partners of u (from CS_u). This ensures the local cost correctly
%   captures the high-distortion pairs near slit boundaries.
%
%   Inputs:
%     u        - domain node index being evaluated
%     F_map    - n x 1 current filling map
%     adj_cell - cell array of graph neighbors, adj_cell{i} = neighbor indices
%     CS_u     - vector of cross-slit partners of node u (may be empty)
%     pts_mat  - n x 2 node coordinates
%     G_ref    - MATLAB graph object
%     cache    - containers.Map distance cache
%
%   Output:
%     worst - max ratio d_image / d_domain over all pairs incident to u

worst = 0;

% Graph neighbors
nbrs = adj_cell{u};
for j = 1:length(nbrs)
    v = nbrs(j);

    d_dom = norm(pts_mat(u,:) - pts_mat(v,:));
    if d_dom < 1e-12, continue; end

    d_img = gd_dlookup(F_map(u), F_map(v), G_ref, cache);
    if isinf(d_img), continue; end

    ratio = d_img / d_dom;
    if ratio > worst
        worst = ratio;
    end
end

% Cross-slit partners
for j = 1:length(CS_u)
    v = CS_u(j);

    d_dom = norm(pts_mat(u,:) - pts_mat(v,:));
    if d_dom < 1e-12, continue; end

    d_img = gd_dlookup(F_map(u), F_map(v), G_ref, cache);
    if isinf(d_img), continue; end

    ratio = d_img / d_dom;
    if ratio > worst
        worst = ratio;
    end
end

end