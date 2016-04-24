function seg_map = draw_segment(seg_img, input_im, seg_name)
input_im = im2double(input_im);
alpha_img(:,:,1) = (1-seg_img(:,:,1)) + seg_img(:,:,1).*input_im(:,:,1);
alpha_img(:,:,2) = (1-seg_img(:,:,1)) + seg_img(:,:,1).*input_im(:,:,2);
alpha_img(:,:,3) = (1-seg_img(:,:,1)) + seg_img(:,:,1).*input_im(:,:,3);
alpha = 0.5;
seg_map = (1-alpha).*input_im + alpha.*alpha_img;
temp_img = double(seg_img(:,:,1));
boundary = bwperim(temp_img,8);
[xx, yy] = find(boundary==1);
for ii = 1:length(xx)
    seg_map(xx(ii), yy(ii),1) = 1;
    seg_map(xx(ii), yy(ii),2) = 1;
    seg_map(xx(ii), yy(ii),3) = 1;
end
imwrite(seg_map, seg_name);

