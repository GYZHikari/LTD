function edges = edges_between(inds)
% connect all the inds
% Copyright (C) Guangyu Zhong
% All rights reserved.
if isempty(inds)
    edges = [];
end
num = length(inds);
mat = tril(ones(num), -1);
[row, col] = find(mat);
edges = [inds(row), inds(col)];