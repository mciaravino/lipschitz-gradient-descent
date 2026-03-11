function [lip, wu, wv] = gd_compute_lip(F_map, E_mat, pts_mat, G_ref, cache)
% GD_COMPUTE_LIP  Compute global Lipschitz constant of filling map F.
%
%   [lip, wu, wv] = gd_compute_lip(F_map, E_mat, pts_mat, G_ref, cache)
%
%   Computes the maximum ratio d_image / d_domain over all graph edges.
%
%   Inputs:
%     F_map   - n x 1 vector, F_map(i) = image node index of node i
%     E_mat   - m x 2 edge list from G.Edges.EndNodes
%     pts_mat - n x 2 node coordinates
%     G_ref   - MATLAB graph object
%     cache   - containers.Map distance cache (see gd_get_row, gd_dlookup)
%
%   Outputs:
%     lip  - maximum ratio d_image / d_domain over all edges
%     wu   - domain node index of worst pair (first)
%     wv   - domain node index of worst pair (second)

lip = 0; wu = -1; wv = -1;

for k = 1:size(E_mat, 1)
    u = E_mat(k,1);
    v = E_mat(k,2);

    d_dom = norm(pts_mat(u,:) - pts_mat(v,:));
    if d_dom < 1e-12, continue; end

    d_img = gd_dlookup(F_map(u), F_map(v), G_ref, cache);
    if isinf(d_img), continue; end

    ratio = d_img / d_dom;
    if ratio > lip
        lip = ratio;
        wu  = u;
        wv  = v;
    end
end

end