function worst = gd_local_lip(u, F_map, adj_cell, pts_mat, G_ref, cache)
% GD_LOCAL_LIP  Worst Lipschitz ratio over edges incident to node u only.
%
%   worst = gd_local_lip(u, F_map, adj_cell, pts_mat, G_ref, cache)
%
%   Only recomputes ratios for edges touching u — used inside the
%   coordinate descent loop for efficiency, since changing F_map(u)
%   only affects these edges.
%
%   Inputs:
%     u        - domain node index being evaluated
%     F_map    - n x 1 current filling map
%     adj_cell - cell array, adj_cell{i} = neighbor indices of node i
%     pts_mat  - n x 2 node coordinates
%     G_ref    - MATLAB graph object
%     cache    - containers.Map distance cache (see gd_get_row, gd_dlookup)
%
%   Output:
%     worst - max ratio d_image / d_domain over edges incident to u

worst = 0;
nbrs  = adj_cell{u};

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

end