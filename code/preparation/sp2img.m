function matrix = sp2img(labels, prior)
% Assign the prior of each patch to the whole image
% Copyright (C) Guangyu Zhong
% All rights reserved.
for i = 1:size(labels,1)
    for j = 1:size(labels,2)
        matrix(i,j) =  prior(labels(i,j));
    end
end