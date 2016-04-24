function svm_model = train_svm(train_dir, gt_dir, train_sp, sp_num_max, bin_num, mode)
train_names = dir([train_dir,'*.bmp']);
fea_hist = [];
gt = [];
for ii = 1:length(train_names)
    name = train_names(ii).name;
    train_im = imread([train_dir, name]);
    [m,n,d] = size(train_im);
    [labels, ~,sp_num, ~, center, ~]  = gene_superpixel(train_dir, name, sp_num_max, train_sp, m, n,[]);
    [~, fea_hist_temp] = gene_fea_hist(im2double(train_im),labels,bin_num,sp_num);
    fea_hist = [fea_hist, fea_hist_temp];
    
    gt_im = im2double(imread([gt_dir, [ name(1:end-4), ['_GT_' mode '.png'] ] ]));
    gt_im = gt_im(:,:,1);
    [xx, yy] = find(gt_im(:,:,1)==1);
    ind = sub2ind(size(labels), xx, yy);
    
    gt_label = zeros(m,n);
    gt_label(ind) = 1;
    gt_temp = zeros(sp_num, 1);
    ind = sub2ind(size(labels), center(:,1), center(:,2));
    gt_temp = gt_label(ind);
    gt = [gt; gt_temp];
end
%% svm training
svm_model = svmtrain(gt,fea_hist', '-s 4 -t 0');
save(['svm_model_' mode '.mat'],'svm_model' );
