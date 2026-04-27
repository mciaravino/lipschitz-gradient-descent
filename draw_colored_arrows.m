function draw_colored_arrows(x, y, dx, dy, varargin)
% DRAW_COLORED_ARROWS  Draw displacement arrows colored by magnitude.
%   Blue = small, red = large.
%   Groups arrows into color buckets for fast rendering.
%
%   draw_colored_arrows(x, y, dx, dy, 'MaxMag', 1)

p = inputParser;
p.addParameter('MaxMag', []);
p.addParameter('NBuckets', 20);
p.parse(varargin{:});
max_mag   = p.Results.MaxMag;
n_buckets = p.Results.NBuckets;

mag   = sqrt(dx.^2 + dy.^2);
moved = mag > 1e-10;
if ~any(moved), return; end

if isempty(max_mag), max_mag = max(mag(moved)); end
if max_mag < 1e-10, return; end

norm_mag = min(mag(moved) / max_mag, 1);
xm = x(moved); ym = y(moved);
dxm = dx(moved); dym = dy(moved);

bucket = min(floor(norm_mag * n_buckets), n_buckets-1) + 1;

for b = 1:n_buckets
    idx = bucket == b;
    if ~any(idx), continue; end
    t   = (b-1) / (n_buckets-1);
    col = [t, 0, 1-t];
    qh  = quiver(xm(idx), ym(idx), dxm(idx), dym(idx), 0, ...
        'LineWidth', 1.2, 'MaxHeadSize', 0.5);
    qh.Color = col;
end

end