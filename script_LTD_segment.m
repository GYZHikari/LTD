%% A implementation for the work proposed by 
%% Risheng Liu, Guangyu Zhong, Junjie Cao, Zhouchen Lin, Shiguang Shan and Zhongxuan Luo.
%% Learn-To-Diffuse (LTD): A New Perspective to Design PDE System for Visual Analysis   
%% Copyright (C) Guangyu Zhong, Junjie Cao and Risheng Liu
%% All rights reserved.
%% Any problem, please contact Guangyu Zhong: guangyuzhonghikari@gmail.com
close all; clear all; clc;
addpath(genpath('code'));
imdir = './img/input/';
imname = '1_19_s.bmp';
spdir = './img/output/sp/';
outdir = './img/output/seg/';
if ~isdir(spdir)
    mkdir(spdir);
end
if ~isdir(outdir)
    mkdir(outdir);
end
SHOW = 1;
%% superpixel parameters
sp_num_max = 800; % max superpixel number
sp_mode = 'generate'; % generate: regenerate data load: load the known data
seed_mode = 'generate';  % generate: regenerate data load: load the known data
%% graph & candidate seeds parameters
% Intra-image Netwrok is a K-nn graph between superpixels. 
WGT_SCALE_PCT = 0.5; % scale of gaussian weights
CONST_Z = 0.1 ; % The contant for ground conductance
numNgh = 5; % K
graphOpts.max_dist_ratio_sp = 1.4;
graphOpts.gauss = 0;
graphOpts.similarity_type = 'aff_similar';
graphOpts.theta = 10;
cluster_num = 2;
cluster_num2 = 30;
%% segmentation parameters
opts.alpha = 0;
opts.seed = 'ones';
opts.seednum = 3;
opts.lambda = 25;
train_mode = 'generate';  % generate: regenerate data load: load the known data
sub_mode = 'sum_r_single_w'; % objective function modes sum : H = L  sum_r_all, sum_r_single: H = L + lambda*R  sum_r_single_w: H = L + lambda*R - W 


%% segment the input image into superpixels
input_im = imread([imdir, imname]);
[m, n, k] = size(input_im);
[superpixels, sp_adjmat, sp_num, sp_inds, sp_center, sp_npix] = gene_superpixel(imdir, imname, sp_num_max, spdir, m, n, sp_mode);
if SHOW
    [label_img] = disp_draw_intra_agglo_real(input_im, superpixels, [], 0, [spdir, imname(1:end-4) '_sp.png'],3) ;
    figure( 'name', 'SLIC');imshow(label_img,[]);
end
[~, sp_fea] = gene_sp_fea(input_im, superpixels);
sp_fea = normal(sp_fea);

%% generate affinity matrix
[aff_mat, Z_SP, edges1, edges2] = gen_intra_net(sp_fea, [],...
    sp_center, WGT_SCALE_PCT,CONST_Z,1,numNgh);
edges(:,1) = edges1(:);
edges(:,2) = edges2(:);
bd = unique([superpixels(1,:)'; superpixels(end,:)'; superpixels(:,1); superpixels(:,end)]);
bd_edges = edges_between(bd);
edges = [bd_edges;edges];
aff_mat = gene_affmat(edges, sp_num, [],sp_fea, graphOpts);
aff_mat = (aff_mat + aff_mat')./2;

%% generate candidate seeds
if strcmp(seed_mode, 'load')&&exist('fore_candi_seeds.mat', 'file')
    load('fore_candi_seeds.mat');
else
    eval = kmeans_fore_candidate(cluster_num,cluster_num2, sp_center, sp_fea, m, n);
end
if SHOW
    [label_img] = disp_draw_intra_agglo_real(input_im, superpixels, eval, 0, [outdir, imname(1:end-4) '_forecandi.png'],3) ;
    figure( 'name', 'foreground candidate seeds');imshow(label_img,[]);
end

%% segmentation samples
switch lower(sub_mode)
    case 'sum'
        opts.objfun = 'sum'; %
        P = [];
        objectness = [];
    case 'sum_r_all'
        opts.objfun = 'regre';
        mode = 'all';
        bin_num = 8;
        objectness = gene_objectness(train_mode, sp_num_max, bin_num, mode, input_im, superpixels);
        obj_im = sp2img(superpixels, objectness);
        if SHOW
            figure;imshow(obj_im);
            imwrite(obj_im, [outdir, imname(1:end-4), mode, '_g'  '.png']);
        end
        P = aff_mat ./ repmat(sum(aff_mat,2),1,size(aff_mat,2));
    case'sum_r_single'
        opts.objfun = 'regre';
        mode = 'single';
        bin_num = 8;
        objectness = gene_objectness(train_mode, sp_num_max, bin_num, mode, input_im, superpixels);
        obj_im = sp2img(superpixels, objectness);
        if SHOW
            figure;imshow(obj_im);
            imwrite(obj_im, [outdir, imname(1:end-4), '_', mode, '_g'  '.png']);
        end
        P = aff_mat ./ repmat(sum(aff_mat,2),1,size(aff_mat,2));
    case'sum_r_single_w'
        opts.objfun = 'regre+w';
        mode = 'single';
        bin_num = 8;
        objectness = gene_objectness(train_mode, sp_num_max, bin_num, mode, input_im, superpixels);
        obj_im = sp2img(superpixels, objectness);
        if SHOW
            figure('name', 'testing results');imshow(obj_im);
            imwrite(obj_im, [outdir, imname(1:end-4), '_', mode, '_g'  '.png']);
        end
        P = aff_mat ./ repmat(sum(aff_mat,2),1,size(aff_mat,2));
end
prior_source = zeros(sp_num,1);
prior_source(eval) = 1;
[ranks, pos_label_inds_choosen, obj_val] = solve_submodular_naive(...
    aff_mat,3, eval, bd, prior_source, input_im, superpixels,opts,P,objectness,SHOW);
if SHOW
    [seed_img] = disp_draw_intra_agglo_real(input_im, superpixels, pos_label_inds_choosen, 0, [outdir, imname(1:end-4), '_', sub_mode, '_seed.png'] , 3);
    figure('name', 'heat source'); imshow(seed_img);
end
seg_ranks = normal(ranks);
seg_ranks(ranks<0.5) = 0;
seg_ranks(ranks>=0.5) = 1;
seg_img = sp2img(superpixels, seg_ranks);
seg_name = [outdir, imname(1:end-4), '_', sub_mode '_seg' '.png'];
seg_map = draw_segment(seg_img, input_im, seg_name);
if SHOW
    figure('name', 'segmentation');imshow(seg_map,[]);
end
