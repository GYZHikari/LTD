function [ranks, pos_label_inds_choosen, obj_val] = solve_submodular_naive(...
    aff_mat,seed_num, pos_label_inds_all, neg_label_inds, prior_source,input_im,superpixels, opts,P,objectness, SHOW)
% solve submodular
% -------
% INPUTS:
%
%   aff_mat             : affinity matrix ([N x N]) where N is the patch number
%   seed_num            : max heat source number
%   pos_label_inds_all  : candidate heat sources
%   prior_source        : prior source function of all the nodes ([N x 1])
%   P                   : transfer matrix ([N x N]) 
%   objectness          : regressed value ([N x 1])
% --------
% OUTPUTS:
%
%   ranks                  : stable temperature ([N x 1])
%   pos_label_inds_choosen : heat sources
%   obj_val                : objective values
% -------------
% Copyright (C) Risheng Liu, Junjie Cao and Guangyu Zhong
% All rights reserved.
num_vtx = size(aff_mat, 1);
num_pos_label = length(pos_label_inds_all);
if seed_num > num_pos_label;
    seed_num = num_pos_label;
    disp(['new posotive seed number: ' num2str(seed_num)]);
end

obj_val = [];
current_pos_label_inds = [];
current_ranks = prior_source;

for i= 1:seed_num
    cand_pos_label_inds = setdiff(pos_label_inds_all, current_pos_label_inds);
    num_cand = length(cand_pos_label_inds);
    cand_ranks = zeros(num_vtx, num_cand);
    cand_obj_val = zeros(num_cand, 1);
    for j = 1:num_cand
        pos_label_inds = [current_pos_label_inds; cand_pos_label_inds(j)];
        switch lower(opts.seed)
            case 'initial'
                pos_label_ranks = prior_source(pos_label_inds);
            case 'ones'
                pos_label_ranks = ones(length(pos_label_inds), 1);
        end
        y_source = prior_source;
        cand_ranks(:, j) = solve_les_with_dirichlet(...
            aff_mat, pos_label_inds, pos_label_ranks,neg_label_inds, y_source, opts);
        if strcmp(opts.objfun, 'regre')||strcmp(opts.objfun, 'regre+w')          
                if i==1
                    cand_obj_val(j) =  sum(cand_ranks(:, j));
                else
                    cand_obj_val(j) = gene_obj_val(cand_ranks(:, j), pos_label_inds, opts, prior_source, P, objectness);
                end
        else
                cand_obj_val(j) = gene_obj_val(cand_ranks(:, j), pos_label_inds, opts, prior_source, P, objectness);
        end
    end
    
    [max_val, max_ind] = max(cand_obj_val);
    if i > 1 && (max_val<obj_val(i - 1))
        break;
    end
    current_pos_label_inds = [current_pos_label_inds; cand_pos_label_inds(max_ind)];
    obj_val = [obj_val; max_val];
    current_ranks = cand_ranks(:, max_ind);
    current_ranks=(current_ranks-min(current_ranks(:)))/(max(current_ranks(:))-min(current_ranks(:)));
    if SHOW
        [label_im] = disp_draw_intra_agglo_real(input_im, superpixels, current_pos_label_inds, 0,[],3) ;
        [rank_im] = sp2img(superpixels, current_ranks);
        if i==1
            h = figure('name', 'heat source & diffusion results'); clf
            h =draw_maps(label_im,rank_im,h);
        else
            h =draw_maps(label_im,rank_im,h);
        end
    end
end
ranks = current_ranks;
pos_label_inds_choosen = current_pos_label_inds;



