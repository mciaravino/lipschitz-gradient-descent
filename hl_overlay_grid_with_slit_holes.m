function pts = hl_overlay_grid_with_slit_holes(slits, N, varargin)

p = inputParser;
p.addParameter('PointSize', 14);
p.addParameter('Color', 'r');
p.addParameter('Tol', 1e-12);
p.addParameter('MaxN', 6);        % <-- NEW: grid fixed to N=6 depth
p.parse(varargin{:});
opt = p.Results;

% Grid resolution fixed to MaxN regardless of current N
refine = 0;
L = 2^(opt.MaxN + 1 + refine);   % 2^(6+1+0) = 2^7 = 128

coords = (0:L)'/L;

[X,Y] = meshgrid(coords, coords);
pts = [X(:), Y(:)];

x = pts(:,1);
y = pts(:,2);

kill = false(size(pts,1),1);

slit_x  = [slits.x]';
slit_y0 = [slits.y0]';
slit_y1 = [slits.y1]';

for k = 1:numel(slit_x)
    on_x = abs(x - slit_x(k)) <= opt.Tol;
    in_y = (y >= slit_y0(k)-opt.Tol) & (y <= slit_y1(k)+opt.Tol);
    kill = kill | (on_x & in_y);
end

pts = pts(~kill,:);

% Plot
hold on
plot(pts(:,1), pts(:,2), '.', 'Color', opt.Color, 'MarkerSize', opt.PointSize)
hold off

fprintf('Grid nodes: %d\n', size(pts,1))

end