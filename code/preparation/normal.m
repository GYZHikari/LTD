function normal_fea = normal(fea)
% fea    m*n    n is the feature number m is the feature dimension
[m,n] = size(fea);
min_fea = repmat(min(fea), m, 1);
max_fea = repmat(max(fea), m, 1);
normal_fea=(fea-min_fea)./(max_fea - min_fea);