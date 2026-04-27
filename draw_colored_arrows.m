function draw_colored_arrows(x, y, dx, dy, varargin)
% DRAW_COLORED_ARROWS  Draw displacement arrows colored by magnitude.
%   Blue = small displacement, red = large displacement.
%   Uses quiver grouped by color buckets for proper arrowheads.
%
%   draw_colored_arrows(x, y, dx, dy, 'MaxMag', 1)

p = inputParser;
p.addParameter('MaxMag', []);
p.parse(varargin{:});
max_mag = p.Results.MaxMag;

mag   = sqrt(dx.^2 + dy.^2);
moved = mag > 1e-10;
if ~any(moved), return; end

if isempty(max_mag), max_mag = max(mag(moved)); end
if max_mag < 1e-10, return; end

norm_mag = min(mag(moved) / max_mag, 1);

xm = x(moved); ym = y(moved);
dxm = dx(moved); dym = dy(moved);

n_buckets = 50;
bucket = min(floor(norm_mag * n_buckets), n_buckets-1) + 1;

for i = 1:length(xm)
    t   = norm_mag(i);
    col = [t, 0, 1-t];
    quiver(xm(i), ym(i), dxm(i), dym(i), 0, ...
        'Color', col, 'LineWidth', 1.2, 'MaxHeadSize', 0.5);
end

end