function row = gd_get_row(src, G_ref, cache)
% GD_GET_ROW  Compute and cache the distance row from source node src.
%
%   row = gd_get_row(src, G_ref, cache)
%
%   If the row for src is already in cache, returns it immediately.
%   Otherwise computes distances(G_ref, src), stores it, and returns it.
%
%   Inputs:
%     src    - source node index (integer)
%     G_ref  - MATLAB graph object
%     cache  - containers.Map with KeyType int32, ValueType any
%
%   Output:
%     row    - 1 x n vector of distances from src to all nodes

key = int32(src);
if ~cache.isKey(key)
    cache(key) = distances(G_ref, src);
end
row = cache(key);
end