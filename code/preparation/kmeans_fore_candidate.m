function eval = kmeans_fore_candidate(cluster_num,cluster_num2, sp_center, sp_fea, m,n)
%	Generate candidate foreground seeds by kmeans
%
% -------
% INPUTS:
%
%   cluster_num     : cluster number for kmeans
%   cluster_num2    : the number of the selected seeds in each cluster
%   sp_center       : center location of each superpixel ([sp_num x 2])
%   sp_fea          : feature of each superpixel ([sp_num x 3])
%   m, n            : height and width of the image
% --------
% OUTPUTS:
%
%  eval            : candidate foreground seed
% Copyright (C) Guangyu Zhong
% All rights reserved.
cluster_id = kmeans(sp_fea, cluster_num);

for ii = 1:cluster_num
    cluster_ii = find(cluster_id==ii);
    [~,candi_center] = kmeans(sp_fea(cluster_ii,:),cluster_num2);
    all_dist = slmetric_pw([candi_center'],[sp_fea(cluster_ii,:)'],'sqdist');
    [~, index] = sort(all_dist,2);
    eval(:,ii) = cluster_ii(index(:,1));
    indx = sp_center(eval(:,ii), 1);
    indy =sp_center(eval(:,ii), 2);
    tmp(ii) = sum((indx - repmat(m/2, cluster_num2, 1)).^2 + (indy - repmat(n/2, cluster_num2, 1)).^2);
end
[~, index] = min(tmp);
eval = eval(:,index);