function G = hl_build_graph_from_points(pts)

n = size(pts,1);

% Build lookup table
key = containers.Map;

for i = 1:n
    k = sprintf('%.10f_%.10f',pts(i,1),pts(i,2));
    key(k) = i;
end

% Grid spacing
dx = unique(diff(unique(pts(:,1))));
dy = unique(diff(unique(pts(:,2))));

h = min([dx;dy]);

% Neighbor offsets
nbr = [
-1 -1
-1  0
-1  1
 0 -1
 0  1
 1 -1
 1  0
 1  1
];

S = [];
T = [];
W = [];

for i = 1:n

    x = pts(i,1);
    y = pts(i,2);

    for k = 1:8

        xn = x + nbr(k,1)*h;
        yn = y + nbr(k,2)*h;

        keystr = sprintf('%.10f_%.10f',xn,yn);

        if isKey(key,keystr)

            j = key(keystr);

            if j > i

                S(end+1) = i;
                T(end+1) = j;
                W(end+1) = hypot(x-xn,y-yn);

            end

        end

    end

end

G = graph(S,T,W,n);

fprintf('Graph edges: %d\n', numedges(G));

end