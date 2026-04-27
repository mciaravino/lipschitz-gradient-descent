function [lip, wu, wv] = global_lip(F_map, keep_i, pts_D_ref, G_X_ref, cache, r_max)
% GLOBAL_LIP  Compute global Lipschitz constant over node pairs.
%
%   [lip, wu, wv] = global_lip(F_map, keep_i, pts_D_ref, G_X_ref, cache)
%   [lip, wu, wv] = global_lip(F_map, keep_i, pts_D_ref, G_X_ref, cache, r_max)
%
%   Exhaustive search over all pairs of nodes in keep_i within
%   Euclidean radius r_max in the domain. Pairs farther than r_max
%   are skipped — they have large denominators and rarely achieve
%   the worst ratio.
%
%   If r_max is not provided, defaults to 10 * grid spacing.
%
%   Denominator: Euclidean distance in domain (pts_D_ref)
%   Numerator:   intrinsic distance in X_N via G_X_ref

if nargin < 6 || isempty(r_max)
    ux = unique(pts_D_ref(:,1));
    h  = min(diff(ux));
    r_max = 3 * h;
end

lip = 0; wu = -1; wv = -1;
n   = length(keep_i);

for a = 1:n
    for b = a+1:n
        ua = keep_i(a); vb = keep_i(b);
        d_dom = norm(pts_D_ref(ua,:) - pts_D_ref(vb,:));
        if d_dom < 1e-12, continue; end
        if d_dom > r_max,  continue; end
        d_img = gd_dlookup(F_map(ua), F_map(vb), G_X_ref, cache);
        if isinf(d_img), continue; end
        ratio = d_img / d_dom;
        if ratio > lip
            lip = ratio; wu = ua; wv = vb;
        end
    end
end

end