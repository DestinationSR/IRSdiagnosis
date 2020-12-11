function [h_rec, m_rec] = cSLRMD(H, W, F, y, lambda, delta)
%
% Compressed sparse and low-rank matrix recovery
%
% INPUTS:
%   H       : number of reflecting elements along the vertical direction
%   W       : number of reflecting elements along the horizontal direction
%   F       : sensing matrix
%   y       : measurement vector
%   lambda  : regularization parameter
%   delta   : noise level
%
% OUTPUTS:
%   h_rec   : recovered channel vector
%   m_rec   : recovered failure mask
%
% NOTICE:   CVX toolbox required. Available on: cvxr.com/cvx

%%
cvx_begin quiet

variable Hchan(H,W) complex
variable mask(H,W) complex

minimize(lambda*norm(mask(:),1) + norm_nuc(Hchan))
norm(y - F*(Hchan(:) + mask(:))) <= delta

cvx_end

h_rec = Hchan(:);
m_rec = vec(mask./Hchan + ones(H,W));

end