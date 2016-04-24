function [superpixels, sp_adj_mat,sp_num, sp_inds, sp_center, sp_npix] = gene_superpixel(imdir, imname, sp_num_max, spdir, m, n, sp_mode)
%	Generate superpixels of an image
%
% -------
% INPUTS:
%
%   imdir           : path of the input image
%   imname          : image name
%   sp_num_max      : max superpixel number 
%   spdir           : path of the superpixel data
%   m, n            : height and width of the image
%   sp_mode         : generate, regenerate superpixel data 
%                     load, load known superpixel data
% --------
% OUTPUTS:
%
%   sp_num          : actual superpixel number 
%   superpixels     : superpixel matrix ([m x n]) 
%   sp_adj_mat      : adjacent matrix ([sp_num x sp_num])
%   sp_inds         : superpixels in each label (cell{sp_num})
%   sp_center       : center location of each superpixel ([sp_num x 2])
%   sp_npix         : pixel number of each superpixel ([sp_num x 1])
% -------------
% Copyright (C) Guangyu Zhong and Junjie Cao
% All rights reserved.

if ~strcmp(imname(end-3:end),'.bmp')
    imname=[imname(1:end-4) '.bmp'];% the slic software support only the '.bmp' image
end
if isempty(sp_mode)
    sp_mode = 'generate';
end
if strcmp(sp_mode, 'generate')
    comm=['.\code\graph\SLICSuperpixelSegmentation' ' ' [imdir imname] ' ' ...
        int2str(10) ' ' int2str(sp_num_max) ' ' spdir];
    system(comm);
end
spname=[spdir imname(1:end-4)  '.dat'];
superpixels=read_dat([m,n],spname); % superpixel label matrix
sp_num=max(superpixels(:)); % the actual superpixel number

sp_adj_mat = build_sp_adjacent_matrix(superpixels,sp_num);

if nargout > 3
    sp_inds=cell(sp_num,1);
    sp_center = zeros(sp_num,2);
    sp_npix  = zeros(sp_num,1);
    for j=1:sp_num,
        sp_inds{j}=find(superpixels==j);
        [r c] = ind2sub([m,n], sp_inds{j});
        sp_npix(j) = length(r);
        sp_center(j,:) = round([mean(r) mean(c)]);
    end
end