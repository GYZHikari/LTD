function weights = fn_weights_gauss(featA, featB, indA, indB, pct_scale,param)
%
%  Compute Gaussian similarity 
%
% -------
% INPUTS:
%
%   featA       : The first Feature vectors ([# of SP x feat dim])
%   featB       : The second Feature vectors ([# of SP x feat dim])
%   edges       : An Mx2 list of M edges indexing into points
%                 of MxM pairwise l2 norm
%   pct_scale   : Scale for gaussian weights
%   param = 1   :  dissimilar
%   param = 2   :similar
% --------
% OUTPUTS:
%
%   weights:    An Mx2 weight list of M edges.
%
% ---------------
% Gunhee Kim
%

%Constants
EPS = 1e-5;

% normal scaling (fixed value at 'scale')
if ~isempty(indA) || ~isempty(indB), 
    [rA cA] = size(indA) ;
    [rB cB] = size(indB) ;
    % error check
    if rA ~= rB || cA ~= cB, 
        error('XX: indA and indB should have the same dimensions.') ;
    end
    indA = indA(:) ; 
    indB = indB(:) ;

    % compute the distances of pairs only.
    featDist = sum((featA(indA,:) - ...
        featB(indB,:)).^2,2) ;

    indA = reshape(indA, [rA cA]) ;
    indB = reshape(indB, [rB cB]) ;

    featDist = reshape(featDist, [rA cA]) ;
% if edge is empty, compute pairwise diatance
else
    featDist = fn_dist_l2_sqrt(featA, featB) ;
end

 % normalize to [0,1]
featDist = normalize(featDist) ;

[rf cf] = size(featDist) ;
nth = round(pct_scale * rf * cf ) ;

nth = min(max(1,nth), (rf*cf)) ;
scale = sort(featDist(:)) ;
scale = scale(nth) ;


% compute Gaussian weights
switch param,
    case 1,
weights = exp(-(featDist/scale)) + EPS;
    case 2,
 weights = exp((featDist/scale)) ;     
end

end


% normalize to [0,1]
function data = normalize(data)

    sr = size(data,1) ;
    sc = size(data,2) ;

    % normalize to [0,1]
    minDist = min(data(:)) ;
    maxDist = max(data(:)) ;
    data  = data - minDist ;
    diff = maxDist - minDist ; 

    if diff==0, 
        weights = ones(sr, sc) ;
    else
        data = data / diff ;
    end

end


