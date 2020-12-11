function Tbs = toeplitz2dAdj(Q, H, W, findIdx)
%
% Adjoint of twofold Toeplitz operator
%
% INPUTS:
%   Q       : input matrix
%   H       : number of reflecting elements along the vertical direction
%   W       : number of reflecting elements along the horizontal direction
%   findIdx : mapping indexes
%
% OUTPUTS:
%   Tbs     : output vector

%%
TbsLen = (H-1)*(2*W-1)+W;

Tbs = zeros(TbsLen,1);

%% extract & sum
for index=1:TbsLen
    
    indexTemp = cell2mat(findIdx(index));
    Tbs(index) = sum(Q(indexTemp));
    
end

end