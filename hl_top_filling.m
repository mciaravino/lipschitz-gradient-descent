function [F, c_bar] = hl_top_filling(pts_D, slits, G_X, pts_X, slit_cols_only)
% HL_TOP_FILLING  Top-lane filling map c_bar_N : D -> X_N.
%
%   Pushes interior points upward to just above the highest slit tip,
%   so everything stacks in the top lane. Gradient descent then
%   spreads points back out into valid positions.
%
%   slit_cols_only = false (default):
%     ALL interior points get pushed to the top lane
%
%   slit_cols_only = true:
%     Only points in the slit column +/- one column on either side
%     get pushed. Everything else maps to itself.
%
%   Boundary nodes always map to themselves.
%   y never exceeds 1 — points already in the top lane stay put.
%
%   Inputs:
%     pts_D         - n_D x 2 full domain grid coordinates
%     slits         - struct array with fields x, y0, y1
%     G_X           - slit domain graph (unused)
%     pts_X         - n_X x 2 slit domain grid coordinates
%     slit_cols_only - logical, default false
%
%   Outputs:
%     F     - n_D x 1 vector, F(i) = index into pts_X of image of pts_D(i)
%     c_bar - n_D x 2 continuous image coordinates before snapping

if nargin < 5, slit_cols_only = false; end

n_D = size(pts_D, 1);
tol = 1e-10;

% Grid spacing
ux_D = unique(pts_D(:,1));
h_D  = min(diff(ux_D));
uy_D = unique(pts_D(:,2));
h_y  = min(diff(uy_D));

% Global top target: just above the highest slit tip across all slits
y_top_global = max([slits.y1]) + h_y;

% Unique slit x-positions
slit_xs = unique([slits.x]);

c_bar = pts_D;  % default: identity

for i = 1:n_D
    x = pts_D(i,1);
    y = pts_D(i,2);

    % Boundary nodes never move
    if abs(x) < tol || abs(x-1) < tol || abs(y) < tol || abs(y-1) < tol
        continue;
    end

    % If slit_cols_only: only act on slit column +/- one column
    if slit_cols_only
        near_slit = any(abs(x - slit_xs) <= h_D * 1.5);
        if ~near_slit
            continue;
        end
    end

    % Target y: just above highest slit tip
    y_target = y_top_global;

    % Already in top lane — stay put
    if y >= y_target
        continue;
    end

    % Push up to top lane
    c_bar(i,2) = y_target;
end

%% Snap to nearest node in pts_X
F = zeros(n_D, 1);
for i = 1:n_D
    dists = (pts_X(:,1) - c_bar(i,1)).^2 + (pts_X(:,2) - c_bar(i,2)).^2;
    [~, F(i)] = min(dists);
end

%% Diagnostics
is_bnd = (abs(pts_D(:,1)) < tol) | (abs(pts_D(:,1)-1) < tol) | ...
         (abs(pts_D(:,2)) < tol) | (abs(pts_D(:,2)-1) < tol);
bnd_idx = find(is_bnd);

n_viol = 0;
for i = 1:length(bnd_idx)
    bi     = bnd_idx(i);
    mapped = pts_X(F(bi),:);
    if ~(abs(mapped(1)) < tol || abs(mapped(1)-1) < tol || ...
         abs(mapped(2)) < tol || abs(mapped(2)-1) < tol)
        n_viol = n_viol + 1;
    end
end

fprintf('Top filling: boundary violations = %d\n', n_viol);
fprintf('Top filling: nodes that moved    = %d / %d\n', ...
    sum(abs(c_bar(:,2) - pts_D(:,2)) > tol), n_D);

end