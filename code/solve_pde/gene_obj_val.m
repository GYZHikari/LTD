function obj_val = gene_obj_val(ranks, label_inds, opt, prior_source, P,objectness)
% Generate objective function value
%
% -------
% INPUTS:
%
%   ranks           : diffusion result ([N x 1]) where N is superpixel
%                     number
%   label_inds      : selected heat source
%   prior_source    : guidence value ([N x 1]) 
%   P               : transfer matrix ([N x N]) 
%   objectness      : regressed value ([N x 1])
% --------
% OUTPUTS:
%
%   obj_val         : objective function value
% -------------
% Copyright (C) Guangyu Zhong
% All rights reserved.
switch lower(opt.objfun)
    case 'sum'
        obj_val = sum(ranks);
    case 'prior'
        obj_val = sum(ranks) - sum(1./(prior_source(label_inds).^2+eps));
    case'regre'
        [temp_seed, temp_id] = intersect([1:length(ranks)]', label_inds);
        candi_seeds = [1:length(ranks)]';
        candi_seeds(temp_id) = [];
        pa_b = P(label_inds, candi_seeds);
        pb = repmat(objectness(label_inds), 1, length(candi_seeds));
        p = pb .* pa_b .* log(pa_b);
        p(isnan(p)) = 0;
        ig = -sum(p(:));
        obj_val = sum(ranks) + opt.lambda*ig;
    case'regre+w'
        [temp_seed, temp_id] = intersect([1:length(ranks)]', label_inds);
        candi_seeds = [1:length(ranks)]';
        candi_seeds(temp_id) = [];
        pa_b = P(label_inds, candi_seeds);
        pb = repmat(objectness(label_inds), 1, length(candi_seeds));
        p = pb .* pa_b .* log(pa_b);
        p(isnan(p)) = 0;
        ig = -sum(p(:));
        obj_val = sum(ranks) + opt.lambda*ig-sum(1./(objectness(label_inds).^2+eps));       
    otherwise
         error('Unknown objective function.');
end