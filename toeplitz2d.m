function Tb = toeplitz2d(u, H, W)
%
% Twofold Toeplitz operator
%
% INPUTS:
%   u       : input vector
%   H       : number of reflecting elements along the vertical direction
%   W       : number of reflecting elements along the horizontal direction
%
% OUTPUTS:
%   Tb      : twofold Toeplitz matrix

%%
TaBlks = cell(H,1);

TaBlks(1) = mat2cell(toeplitz([u(1);conj(u(2:W))]),W,W);

for a=2:H
    
    startIdx = (a-2)*(2*W-1)+W+1;
    rowIdx = startIdx+W-1 : -1 : startIdx;
    colIdx = startIdx+W-1 : startIdx+2*W-2;
    
    TaBlks(a) = mat2cell(toeplitz(u(colIdx),u(rowIdx)),W,W);
    
end

n = length(TaBlks);
Tb = cell2mat(TaBlks(toeplitz(1:n)));
Tb = (Tb.*tril(ones(H*W),-1))' + Tb.*tril(ones(H*W),0);

end