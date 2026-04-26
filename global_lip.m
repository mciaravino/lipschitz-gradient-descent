function [lip, wu, wv] = global_lip(F_map, keep_i, pts_D_ref, G_X_ref, cache)
% GLOBAL_LIP  Compute global Lipschitz constant over slit-column pairs.
%
%   [lip, wu, wv] = global_lip(F_map, keep_i, pts_D_ref, G_X_ref, cache)
%
%   Exhaustive search over all pairs of nodes in keep_i.
%   Denominator: Euclidean distance in domain (pts_D_ref)
%   Numerator:   intrinsic distance in X_N via G_X_ref

lip = 0; wu = -1; wv = -1;
n = length(keep_i);

for a = 1:n
    for b = a+1:n
        ua = keep_i(a); vb = keep_i(b);
        d_dom = norm(pts_D_ref(ua,:) - pts_D_ref(vb,:));
        if d_dom < 1e-12, continue; end
        d_img = gd_dlookup(F_map(ua), F_map(vb), G_X_ref, cache);
        if isinf(d_img), continue; end
        ratio = d_img / d_dom;
        if ratio > lip
            lip = ratio; wu = ua; wv = vb;
        end
    end
end

end