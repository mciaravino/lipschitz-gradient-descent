function [F, c_bar] = hl_lane_filling(pts_D, slits, G_X, pts_X)
% HL_LANE_FILLING  Simple filling map c_bar_N : D -> X_N.
%
%   For each point (x,y) in D:
%     - If x matches a slit x-position: push left one grid step
%     - Otherwise: map to itself
%
%   This is well-defined for any N since slit x-positions are
%   determined by N. Boundary maps to itself since slit x-positions
%   are never at x=0 or x=1.
%
%   Inputs:
%     pts_D  - n_D x 2 full domain grid coordinates
%     slits  - struct array with fields x, y0, y1
%     G_X    - slit domain graph (unused)
%     pts_X  - n_X x 2 slit domain grid coordinates
%
%   Outputs:
%     F     - n_D x 1 vector, F(i) = index into pts_X of image of pts_D(i)
%     c_bar - n_D x 2 continuous image coordinates before snapping

n_D = size(pts_D, 1);
tol = 1e-10;

% Grid spacing
ux_D = unique(pts_D(:,1));
h_D  = min(diff(ux_D));

% Unique slit x-positions
slit_xs = unique([slits.x]);

c_bar = pts_D;  % default: identity

for i = 1:n_D
    x = pts_D(i,1);
    y = pts_D(i,2);

    % Check if this point sits on a slit column
    for k = 1:numel(slits)
        if abs(x - slits(k).x) < tol
            % Only push if within this slit's y-range
            if y >= slits(k).y0 && y <= slits(k).y1
                c_bar(i,1) = x - h_D;  % push left one grid step
            end
            break;
        end
    end
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
bnd_idx  = find(is_bnd);
n_viol   = 0;
for i = 1:length(bnd_idx)
    bi     = bnd_idx(i);
    mapped = pts_X(F(bi),:);
    if ~(abs(mapped(1)) < tol || abs(mapped(1)-1) < tol || ...
         abs(mapped(2)) < tol || abs(mapped(2)-1) < tol)
        n_viol = n_viol + 1;
    end
end

fprintf('Lane filling: boundary violations = %d\n', n_viol);
fprintf('Lane filling: nodes that moved    = %d / %d\n', ...
    sum(abs(c_bar(:,1) - pts_D(:,1)) > tol), n_D);

end