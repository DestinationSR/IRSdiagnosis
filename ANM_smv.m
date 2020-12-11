function [h_rec, m_rec] = ANM_smv(H, W, rho, lambda, tau, maxIter, F, y)
%
% SMV atomic norm minimization
%
% INPUTS:
%   H       : number of reflecting elements along the vertical direction
%   W       : number of reflecting elements along the horizontal direction
%   rho     : penalty parameter
%   lambda  : regularization parameter of failure mask
%   tau     : regularization parameter of channel
%   maxIter : max number of iterations
%   F       : sensing matrix
%   y       : measurements
%
% OUTPUTS:
%   h_rec   : recovered channel vector
%   m_rec   : recovered failure mask

%%
W_iter = diag([H*(W:-1:1) kron(H-1:-1:1,[1:W W-1:-1:1])]);
e1 = zeros((H-1)*(2*W-1)+W,1);
e1(1) = 1;

Z = zeros(H*W+1);
Lambda = zeros(H*W+1);
d = zeros(H*W,1);

%% cache indexes
TbsLen = (H-1)*(2*W-1)+W;
uTmp = ((1:TbsLen) + (1j*(1:TbsLen))).';
Tb = toeplitz2d(uTmp,H,W);
findIndex = cell(TbsLen,1);

for index=1:TbsLen
    
    findIndex(index,:) = mat2cell(Tb == index + 1j*index,H*W,H*W);
    
end

%% ADMM iterations
dOld = zeros(H*W,1);

for nIter=1:maxIter
    
    v = Z(H*W+1,H*W+1) + (Lambda(H*W+1,H*W+1) - tau/2)/rho;
    
    u = W_iter\(toeplitz2dAdj(Z(1:H*W,1:H*W) + Lambda(1:H*W,1:H*W)/rho,H,W,findIndex)) - tau*e1/2/rho;
    
    h_chan = (F'*F + 2*rho*eye(H*W))\(F'*(y-F*d) + 2*Lambda(1:H*W,H*W+1) + 2*rho*Z(1:H*W,H*W+1));
    
    d = lasso_ADMM(F,y-F*h_chan,lambda);
    
    Wtmp = [toeplitz2d(u,H,W),h_chan;h_chan',v];
    Q = Wtmp - Lambda/rho;
    [V_eig,E_eig] = eig((Q+Q')/2);
    e = diag(E_eig);
    idx = (e>0);
    Z = V_eig(:,idx)*diag(e(idx))*V_eig(:,idx)';
    Z = (Z+Z')/2;
    
    Lambda = Lambda + rho*(Z - Wtmp);
    
    if((norm(d - dOld) < 1e-3) && (nIter >= 200))
        break;
    end
    
    dOld = d;
    
end

h_rec = h_chan;
m_rec = d./h_chan + ones(H*W,1);

end