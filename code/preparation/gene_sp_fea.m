function [rgb_fea, lab_fea] = gene_sp_fea(image, label)

[m,n,k] = size(image);
sp_num = max(label(:));
image_vals=reshape(image, m*n, k);
rgb_vals=zeros(sp_num,1,3);
for i=1:sp_num
    tmp = find(label==i);
    rgb_vals(i,1,:) = mean(image_vals(tmp,:),1);
end
lab_vals = colorspace('Lab<-', rgb_vals);
lab_fea=reshape(lab_vals,sp_num,3);% feature for each superpixel
rgb_fea=reshape(rgb_vals,sp_num,3);
