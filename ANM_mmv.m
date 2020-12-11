function [H_rec, m_rec] = ANM_mmv(H, W, NantRX, rho, lambda, tau, maxIter, F, Y)
%
% MMV atomic norm minimization
%
% INPUTS:
%   H       : number of reflecting elements along the vertical direction
%   W       : number of reflecting elements along the horizontal direction
%   NantRX  : number of antennas at the RX
%   rho     : penalty parameter
%   lambda  : regularization parameter of failure mask
%   tau     : regularization parameter of channel
%   maxIter : max number of iterations
%   F       : sensing matrix
%   Y       : measurements
%
% OUTPUTS:
%   H_rec   : recovered channel vectors
%   m_rec   : recovered failure mask

%%
W_iter = diag([H*(W:-1:1) kron(H-1:-1:1,[1:W W-1:-1:1])]);
e1 = zeros((H-1)*(2*W-1)+W,1);
e1(1) = 1;

Z = zeros(H*W + NantRX);
Lambda = zeros(H*W + NantRX);
D = zeros(H*W,NantRX);

%% cache indexes
TbsLen = (H-1)*(2*W-1)+W;
uTmp = ((1:TbsLen) + (1j*(1:TbsLen))).';
Tb = toeplitz2d(uTmp,H,W);
findIndex = cell(TbsLen,1);

for index=1:TbsLen
    
    findIndex(index,:) = mat2cell(Tb == index + 1j*index,H*W,H*W);
    
end

%% ADMM iterations
Dold = zeros(H*W,NantRX);

for nIter=1:maxIter
    
    V = Z(H*W+1:end,H*W+1:end) + (Lambda(H*W+1:end,H*W+1:end) - tau/2*eye(NantRX))/rho;
    
    u = W_iter\(toeplitz2dAdj(Z(1:H*W,1:H*W) + Lambda(1:H*W,1:H*W)/rho,H,W,findIndex)) - tau*e1/2/rho;
    
    Hchan = (F'*F + 2*rho*eye(H*W))\(F'*(Y-F*D) + 2*Lambda(1:H*W,H*W+1:end) + 2*rho*Z(1:H*W,H*W+1:end));
    
    D = MFOCUSS(F,Y-F*Hchan,lambda);
    
    Wtmp = [toeplitz2d(u,H,W),Hchan;Hchan',V];
    Z = Wtmp - Lambda/rho;
    [V_eig,E_eig] = eig((Z+Z')/2);
    e = diag(E_eig);
    idx = (e>0);
    Z = V_eig(:,idx)*diag(e(idx))*V_eig(:,idx)';
    Z = (Z+Z')/2;
    
    Lambda = Lambda + rho*(Z - Wtmp);
    
    if((norm(D - Dold,'fro')/NantRX < 1e-3) && (nIter>=200))
        break;
    end
    
    Dold = D;
    
end

H_rec = Hchan;
m_rec = mean(D./Hchan + ones(H*W,NantRX),2);

end