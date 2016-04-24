function affmat = gene_affmat(edges, spnum, sp_center,sp_fea, opts)
% Generate affinity matrix
%
% -------
% INPUTS:
%
%   edges           : connected edges
%   spnum           : superpixel number
%   sp_center       : center location of each superpixel ([sp_num x 2]) 
%   sp_fea          : feature of each superpixel ([sp_num x 3]) 
%   opts            : graph options
% --------
% OUTPUTS:
%
%   affmat          : affinity matrix ([spnum x spnum])
% -------------
% Copyright (C) Guangyu Zhong
% All rights reserved.


gauss = getoptions(opts, 'gauss', 0);
affmat = zeros(spnum,spnum);
row = edges(:,1); col = edges(:,2);
ind = sub2ind(size(affmat), row, col);
similarity_type = getoptions(opts,'similarity_type','aff_similar');

tmp = sum( (sp_fea(row,:) - sp_fea(col,:)).^2, 2);
valDistances = sqrt(tmp)+eps;
valDistances = normal(valDistances);

switch  lower(similarity_type)
    case 'dissimilar'
        weights = exp(valDistances/max(valDistances));
    case 'similar'
        weights = exp(-valDistances/max(valDistances));
    case 'aff_similar'
        weights = exp(-opts.theta*valDistances);
    otherwise
        warning('not implement!!')
end

if  gauss
    max_dist_ratio_sp = getoptions(opts, 'max_dist_ratio_sp', 1.0);
    tmp = sqrt( sum( (sp_center(row,:) - sp_center(col,:)).^2, 2));
    w = exp( -tmp.^2/(max_dist_ratio_sp * max(tmp).^2) );
else
    w = ones(numel(row),1);
end

affmat(ind) = weights.*w;
ind = sub2ind(size(affmat), col, row);
affmat(ind) = weights.*w;




