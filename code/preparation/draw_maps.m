function h = draw_maps(label_im,rank_im,h)
% draw seed and rankimage
% GuangyuZhong 08/10/2013
figure(h);clf;
subplot(1,2,1);

imshow(label_im);
subplot(1,2,2);

imshow(rank_im);
   