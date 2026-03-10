function d = gd_dlookup(i, j, G_ref, cache)
% GD_DLOOKUP  Return intrinsic distance between nodes i and j.
%
%   d = gd_dlookup(i, j, G_ref, cache)
%
%   Checks cache for a row from i or j before computing anything new.
%   If neither row is cached, computes and caches the row for i.
%
%   Inputs:
%     i, j   - node indices
%     G_ref  - MATLAB graph object
%     cache  - containers.Map with KeyType int32, ValueType any
%
%   Output:
%     d      - scalar intrinsic distance between nodes i and j

if cache.isKey(int32(i))
    row = cache(int32(i));
    d   = row(j);
elseif cache.isKey(int32(j))
    row = cache(int32(j));
    d   = row(i);
else
    row = distances(G_ref, i);
    cache(int32(i)) = row;
    d = row(j);
end
end