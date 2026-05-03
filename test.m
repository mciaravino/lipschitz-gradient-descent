N    = 3;
MaxN = 6;
r    = 1./(1:(N+2));
slits = hl_draw_slits_iter(r, N, 'DrawGrid', false);
close all;

figure('Position', [100 100 600 600]);
hold on; axis equal;
xlim([0 1]); ylim([0 1]);
title(sprintf('$X_%d$', N), 'Interpreter', 'latex');
xlabel('x'); ylabel('y');

for k = 1:numel(slits)
    plot([slits(k).x slits(k).x], [slits(k).y0 slits(k).y1], ...
        'k-', 'LineWidth', 2);
end
box on;
hold off;