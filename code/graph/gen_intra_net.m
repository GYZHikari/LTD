function [locNet, Z_SP, r_ind_f, sort_ind] = gen_intra_net( descSpCLR, descSpTEX, cntSeg, ...
    WGT_SCALE_PCT, CONST_Z,param,numNgh) 
%
%	Generate Intra-network of an image
%
% -------
% INPUTS:
%
%   descSpCLR           : Color descriptors of SPs ([M x d]) 
%                         where M: # ofSPs, d: feat dimension (ex. rgb: 3)
%   descSpTEX           : Texture descriptors of SPs. (same with 'descSpCLR')
%   cntSeg              : Center coordinates of SP ([M x 2]) (r,c)-indices 
%   numNgh              : Intra-image Netwrok is a K-nn graph between superpixels. 
%                         where K = numNgh
%   WGT_SCALE_PCT       : Sigma for Gaussian similarity
%   CONST_Z             : constant for ground conductance
%   param  = 1          :dissimilar
%   param  = 2          :similar
% --------
% OUTPUTS:
%
%   locNet              : Intra-network of image (M x M)
%   Z_SP                : ground conductance of each SP ([M x 1])
%
% -------------
%
% Gunhee Kim 
% adjusted by Guangyu Zhong

% Similarity threshold. Set to 0 if the similarity is below this.
SIM_THRES = 0.01 ;

nSP = size(cntSeg,1) ;

% I implement two different way to build the intra-network of
% an image. Use spatial information or not. 
% (1) Ignore spatial information
% [locNet] = get_net_feat(descSpCLR, C_INTRA_NGH, WGT_SCALE_PCT, SIM_THRES) ;

% (2) First filter out by location, and then apply feature
% similarity.
% numNgh = round(C_INTRA_NGH*log(nSP)) ;
[locNet, Z_SP ,r_ind_f, sort_ind] = get_net_spatial(descSpCLR, descSpTEX, cntSeg, numNgh, WGT_SCALE_PCT,param) ;
Z_SP = CONST_Z * Z_SP ;

end

% (1) Ignore spatial information
% It turned out that this approach is not good.
function [locNet] = get_net_feat(descSpCLR, numNeighbor, ...
                       WGT_SCALE_PCT, SIM_THRES,param) 

    % Pairwise similarity
    pair_dist_f = fn_weights_gauss(descSpCLR, descSpCLR, [], [], ...
        WGT_SCALE_PCT,param) ;
    pair_dist_f(pair_dist_f<SIM_THRES) = 0 ;
    % make it k-nn graph. 
    locNet = fn_get_knn_matrix(pair_dist_f, numNeighbor) ;

end

% (2) First filter out by location, and then applyf feature similarity.
function [locNet, Z_SP, r_ind_f, sort_ind] = get_net_spatial(descSpCLR, descSpTEX, cntSeg, numNeighbor, ...
                                                WGT_SCALE_PCT,param) 
    nSP = size(cntSeg,1) ;

    Extension = max(round(numNeighbor*2), 1) ;
    locK = min(nSP, numNeighbor*Extension) ;

    % For each SP, find the nearest neighbors in spatial space.
    pair_dist_l = fn_dist_l2(double(cntSeg)') ;
    [sort_val_l sort_ind_l] = sort(pair_dist_l,2) ;
    sort_val_l(:,1) = [];
    sort_ind_l(:,1) = [];

    r_ind_l = repmat([1:nSP]',1,locK) ;
    sort_ind_l = sort_ind_l(:,1:locK) ;
    sort_val_l = sort_val_l(:,1:locK) ;
    % compute similarity between them using color
    % (1) normal scaling
%     pair_dist_f = fn_weights_gauss(descSpCLR, descSpCLR, ...
%             r_ind_l, sort_ind_l, WGT_SCALE_PCT,param) ;
    % (2) local scaling
    pair_dist_f = fn_weights_gauss_ls(descSpCLR, ...
           r_ind_l, sort_ind_l, max(round(numNeighbor),1)) ;

    % compute similarity between them using texture
    if ~isempty(descSpTEX),
        % (1) normal scaling
        pair_dist_t = fn_weights_gauss(descSpTEX, descSpTEX, ...
               r_ind_l, sort_ind_l, WGT_SCALE_PCT,param) ;
        % (2) local scaling
        %pair_dist_t = fn_weights_gauss_ls(descSpTEX, ...
        %    r_ind_l, sort_ind_l, numNeighbor) ;
        pair_dist_f = pair_dist_f + 0.5*pair_dist_t ;
        %pair_dist_f = pair_dist_f/2 ;
    end

    pair_dist_f = reshape(pair_dist_f, [nSP locK]) ;
    
    W = zeros(nSP,nSP) ;
   
   W(sub2ind([nSP nSP], r_ind_l(:), sort_ind_l(:))) = sort_val_l(:) ;

%     sort_ind = sort_ind_l(sub2ind([nSP locK], r_ind_l(:), ...
%             sort_ind_l(:))) ;
%      W = zeros(nSP,nSP) ;
%     W(sub2ind([nSP nSP], r_ind_l(:), sort_ind(:))) = sort_val_l(:) ;

    % For each SP, find the nearest neighbors in feature space.
   [sort_val_f sort_ind_f] = sort(pair_dist_f,2, 'descend') ;
    r_ind_f = repmat([1:nSP]',1,numNeighbor) ;
   sort_ind_f = sort_ind_f(:,1:numNeighbor) ;
  sort_val_f = sort_val_f(:,1:numNeighbor) ;

    % std of each SP
    std_SP = std(sort_val_f,[],2) ;
    % mean of each SP
    mean_SP = mean(sort_val_f,2) ;
    % Z of each SP is defined by std_SP .* mean_SP
    Z_SP = std_SP .* mean_SP ;

   sort_ind = sort_ind_l(sub2ind([nSP locK], r_ind_f(:), ...
             sort_ind_f(:))) ;
    locNet = zeros(nSP,nSP) ;
   
   locNet(sub2ind([nSP nSP], r_ind_f(:), sort_ind(:))) = sort_val_f(:) ;
    locNet = sparse(locNet) ;

end

