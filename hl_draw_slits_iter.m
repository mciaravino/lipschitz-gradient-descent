function [slits] = hl_draw_slits_iter(r, N, varargin)
%HL_DRAW_SLITS_ITER Draw stage-N dyadic slit domain with levels 0..N.
%
% Convention (level ell):
%  - cells have side length 2^{-ell}
%  - in each cell insert one vertical slit at the cell center
%  - slit length = r(ell+2) * 2^{-ell}
%
% This indexing makes the level-0 slit use r2, so with r_n=1/n it has length 1/2.

p = inputParser;
p.addParameter('DrawGrid', false, @(x)islogical(x) && isscalar(x));
p.addParameter('LineWidth', 1.0, @(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter('GridAlpha', 0.15, @(x)isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.parse(varargin{:});
opt = p.Results;

% Need r(2..N+2)
if numel(r) < (N+2)
    error("r must have length at least N+2. Need length(r)>=%d.", N+2);
end

figure; hold on; axis equal;
xlim([0 1]); ylim([0 1]);
set(gca,'YDir','normal');
box on;
title(sprintf('Dyadic slit domain, stage N=%d (levels 0..%d)', N, N));
xlabel('x'); ylabel('y');

% Optional: draw dyadic grid (levels 0..N)
if opt.DrawGrid
    for ell = 0:N
        m = 2^ell;
        xs = (0:m)/m;
        ys = (0:m)/m;
        for k = 2:numel(xs)-1
            plot([xs(k) xs(k)],[0 1],'k-','LineWidth',0.5,'Color',[0 0 0 opt.GridAlpha]);
        end
        for k = 2:numel(ys)-1
            plot([0 1],[ys(k) ys(k)],'k-','LineWidth',0.5,'Color',[0 0 0 opt.GridAlpha]);
        end
    end
end

slits = struct('ell',{},'x',{},'y0',{},'y1',{});
idx = 0;

for ell = 0:N
    m = 2^ell;
    cellSize = 1/m;                 % 2^{-ell}
    slitLen  = r(ell+2) * cellSize; % <-- shifted indexing

    for ix = 0:(m-1)
        for iy = 0:(m-1)
            x0 = ix*cellSize; x1 = (ix+1)*cellSize;
            y0 = iy*cellSize; y1 = (iy+1)*cellSize;

            xc = (x0+x1)/2;
            yc = (y0+y1)/2;

            ya = max(y0, yc - slitLen/2);
            yb = min(y1, yc + slitLen/2);

            idx = idx + 1;
            slits(idx).ell = ell;
            slits(idx).x  = xc;
            slits(idx).y0 = ya;
            slits(idx).y1 = yb;

            plot([xc xc],[ya yb],'k-','LineWidth',opt.LineWidth);
        end
    end
end

hold off;
end