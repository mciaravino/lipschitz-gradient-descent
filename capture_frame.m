function capture_frame(v, F_map, iter_num, lip_val, wu_p, wv_p, ...
                        N, pts_D, pts_X, slits, is_bnd, G_X)
% CAPTURE_FRAME  Draw current filling state and write frame to video.
%
%   capture_frame(v, F_map, iter_num, lip_val, wu_p, wv_p, ...
%                 N, pts_D, pts_X, slits, is_bnd, G_X)
%
%   Inputs:
%     v         - VideoWriter object (already open)
%     F_map     - current filling map (index into pts_X)
%     iter_num  - iteration number (0 = initial)
%     lip_val   - current Lip value
%     wu_p, wv_p - worst pair domain node indices (-1 if none)
%     N         - dyadic depth
%     pts_D     - domain grid coordinates
%     pts_X     - slit domain grid coordinates
%     slits     - slit struct array
%     is_bnd    - logical boundary mask for pts_D
%     G_X       - slit domain graph

tol   = 1e-10;
int_D = find(~is_bnd);

img_p = pts_X(F_map(int_D), :);
dx_p  = img_p(:,1) - pts_D(int_D,1);
dy_p  = img_p(:,2) - pts_D(int_D,2);

clf;
hold on; axis equal;
xlim([0 1]); ylim([0 1]);

if iter_num == 0
    if iter_num == 0
        ttl = ['$\bar{c}_' num2str(N) '$ --- initial (Lip = ' sprintf('%.4f', lip_val) ')'];
    else
        ttl = ['$\bar{c}_' num2str(N) '$ after ' num2str(iter_num) ' iterations (Lip = ' sprintf('%.4f', lip_val) ')'];
    end
    title(ttl, 'Interpreter', 'latex');
else
    title(sprintf('$\\bar{c}_%d$ after %d iterations  (Lip = %.4f)', ...
        N, iter_num, lip_val), 'Interpreter', 'latex');
end
xlabel('x'); ylabel('y');

% Slits
for k = 1:numel(slits)
    plot([slits(k).x slits(k).x], [slits(k).y0 slits(k).y1], ...
        'k-', 'LineWidth', 2);
end

% Grid
plot(pts_X(:,1), pts_X(:,2), '.', 'Color', [0.88 0.88 0.88], 'MarkerSize', 3);

% Colored displacement arrows
draw_colored_arrows(pts_D(int_D,1), pts_D(int_D,2), dx_p, dy_p, 'MaxMag', 1);

% Worst pair and shortest path
if wu_p > 0 && wv_p > 0
    path_p = shortestpath(G_X, F_map(wu_p), F_map(wv_p));
    plot(pts_X(path_p,1), pts_X(path_p,2), 'g-', 'LineWidth', 2.5);
    plot(pts_X(F_map(wu_p),1), pts_X(F_map(wu_p),2), 'ro', ...
        'MarkerSize', 14, 'LineWidth', 2.5);
    plot(pts_X(F_map(wv_p),1), pts_X(F_map(wv_p),2), 'rs', ...
        'MarkerSize', 14, 'LineWidth', 2.5);
end

hold off;
drawnow;
frame = getframe(gcf);
for rep = 1:5
    writeVideo(v, frame);
end

end