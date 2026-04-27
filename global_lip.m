function [lip, wu, wv] = global_lip(F_map, keep_i, pts_D_ref, G_X_ref, cache, r_max)
% GLOBAL_LIP  Compute global Lipschitz constant over node pairs.
%
%   Nodes are sorted by x so the inner loop can break early once
%   x-distance exceeds r_max, giving near-linear average performance.
%
%   r_max defaults to 3 * grid spacing if not provided.

if nargin < 6 || isempty(r_max)
    ux = unique(pts_D_ref(:,1));
    h  = min(diff(ux));
    r_max = 3 * h;
end

lip = 0; wu = -1; wv = -1;

% Sort by x coordinate for early break
[~, ord] = sort(pts_D_ref(keep_i, 1));
sorted_i = keep_i(ord);
n = length(sorted_i);

for a = 1:n
    ua = sorted_i(a);
    for b = a+1:n
        vb = sorted_i(b);
        dx = pts_D_ref(vb,1) - pts_D_ref(ua,1);
        if dx > r_max, break; end  % sorted by x so all further b are too far
        d_dom = norm(pts_D_ref(ua,:) - pts_D_ref(vb,:));
        if d_dom < 1e-12 || d_dom > r_max, continue; end
        d_img = gd_dlookup(F_map(ua), F_map(vb), G_X_ref, cache);
        if isinf(d_img), continue; end
        ratio = d_img / d_dom;
        if ratio > lip
            lip = ratio; wu = ua; wv = vb;
        end
    end
end

end