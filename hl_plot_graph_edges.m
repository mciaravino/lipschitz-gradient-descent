function hl_plot_graph_edges(G, pts, slits, varargin)
% Plot graph edges on top of the current slit drawing.
p = inputParser;
p.addParameter('MaxEdges', 5000, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.parse(varargin{:});
opt = p.Results;

hold on;

E = G.Edges.EndNodes;
m = size(E,1);
if m > opt.MaxEdges
    idx = randperm(m, opt.MaxEdges);
    E = E(idx,:);
end

for k = 1:size(E,1)
    a = pts(E(k,1),:);
    b = pts(E(k,2),:);
    plot([a(1) b(1)],[a(2) b(2)], '-', 'Color', [0 0.6 0], 'LineWidth', 0.5);
end

% redraw slits on top (black)
for k = 1:numel(slits)
    plot([slits(k).x slits(k).x],[slits(k).y0 slits(k).y1], 'k-', 'LineWidth', 1.5);
end

hold off;
end