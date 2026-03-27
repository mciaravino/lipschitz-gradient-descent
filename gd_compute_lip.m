function [lip, wu, wv] = gd_compute_lip(F_map, E_mat, CS_mat, pts_mat, G_ref, cache)
% GD_COMPUTE_LIP  Compute global Lipschitz constant of filling map F.
%
%   [lip, wu, wv] = gd_compute_lip(F_map, E_mat, CS_mat, pts_mat, G_ref, cache)
%
%   Checks both graph edges (E_mat) and cross-slit pairs (CS_mat).
%   Cross-slit pairs are nodes that straddle a slit — they are close
%   in Euclidean distance but far in intrinsic distance, and are the
%   critical pairs for detecting high Lipschitz ratios.
%
%   Inputs:
%     F_map   - n x 1 filling map (F_map(i) = image node of node i)
%     E_mat   - m x 2 graph edge list from G.Edges.EndNodes
%     CS_mat  - p x 2 cross-slit pair list
%     pts_mat - n x 2 node coordinates
%     G_ref   - MATLAB graph object
%     cache   - containers.Map distance cache
%
%   Outputs:
%     lip  - maximum ratio d_image / d_domain
%     wu   - worst pair node 1
%     wv   - worst pair node 2

lip = 0; wu = -1; wv = -1;

% Check graph edges
for k = 1:size(E_mat, 1)
    u = E_mat(k,1);
    v = E_mat(k,2);

    d_dom = norm(pts_mat(u,:) - pts_mat(v,:));
    if d_dom < 1e-12, continue; end

    d_img = gd_dlookup(F_map(u), F_map(v), G_ref, cache);
    if isinf(d_img), continue; end

    ratio = d_img / d_dom;
    if ratio > lip
        lip = ratio; wu = u; wv = v;
    end
end

% Check cross-slit pairs
for k = 1:size(CS_mat, 1)
    u = CS_mat(k,1);
    v = CS_mat(k,2);

    d_dom = norm(pts_mat(u,:) - pts_mat(v,:));
    if d_dom < 1e-12, continue; end

    d_img = gd_dlookup(F_map(u), F_map(v), G_ref, cache);
    if isinf(d_img), continue; end

    ratio = d_img / d_dom;
    if ratio > lip
        lip = ratio; wu = u; wv = v;
    end
end

end