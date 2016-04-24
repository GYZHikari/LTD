function objectness = gene_objectness(train_mode, sp_num_max, bin_num, mode, input_im, superpixels)
% train & test the regression
% -------
% INPUTS:
%
%   train_mode      : load, load the known SVR 
%                     generate, retrain SVR                   
%   sp_num_max      : max superpixel number 
%   bin_num         : bin number of color histogram
%   mode            : regression mode for training  
%                     all, regard all the objects as the ground truth
%                     single, regard a particle object as the ground truth
%   input_im        : input image
%   superpixels     : superpixel matrix  
% --------
% OUTPUTS:
%
%   objectness      : regression results ([sp_num x 1])
% -------------
% Copyright (C) Guangyu Zhong
% All rights reserved.

switch lower(train_mode)
    case 'load'
        load(['svm_model_' mode '.mat']);
    case'generate'
        train_dir = './train/input/';
        gt_dir = './train/gt/';
        train_sp = './train/sp/';
        svm_model = train_svm(train_dir, gt_dir, train_sp, sp_num_max, bin_num, mode);
end
sp_num = max(superpixels(:));
[~, fea_hist_temp] = gene_fea_hist(im2double(input_im),superpixels,bin_num,sp_num);
[objectness,accuracy,decision_values] = svmpredict(rand(sp_num,1),fea_hist_temp',svm_model);
objectness = normal(objectness);
objectness(objectness<0.3) = 0;