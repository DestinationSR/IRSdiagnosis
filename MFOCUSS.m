function X = MFOCUSS(F, Y, lambda)
%
% MFOCUSS algorithm for MMV sparse recovery problem: Y = FX + W
%
% INPUTS:
%   F       : sensing matrix
%   Y       : measurement matrix
%   lambda	: regularization parameter.
%             Generally it is close to the noise variance.
%             In the noiseless cases, simply setting lambda = 1e-10 leads to good performance.
%             In noisy cases, need the modified L-curve method to find optimal lambda.
%             See Ref [1] for details.
%
% OUTPUTS:
%   X       : estimated solution matrix
%
% REFERENCE:
%   [1] Cotter, Shane F., et al. "Sparse solutions to linear inverse problems
%       with multiple measurement vectors." IEEE Transactions on Signal Processing
%       53.7 (2005): 2477-2488.
%
%   Mainbody was written by Zhilin Zhang: https://sites.google.com/site/researchbyzhang/software

%%
% Dimension of the Problem
[~,N] = size(F);
[K,N_RX] = size(Y);

% Default Control Parameters
PRUNE_GAMMA = 1e-4;        % threshold for prunning small gamma_i
p           = 1;           % p-norm
EPSILON     = 1e-5;        % threshold for stopping iteration
MAX_ITERS   = 500;         % maximum iterations

% Initializations
gamma = ones(N,1);         % initialization of gamma_i
keep_list = (1:N)';        % record the index of nonzero gamma_i
m = length(keep_list);     % number of nonzero gamma_i
mu = zeros(N,N_RX);        % initialization of the solution matrix
count = 0;                 % record iterations

% Learning loop
while (1)
    
    % =========== Prune weights as their hyperparameters go to zero ===========
    if (min(gamma) < PRUNE_GAMMA)
        index = find(gamma > PRUNE_GAMMA);
        gamma = gamma(index);
        F = F(:,index); % corresponding columns in F
        keep_list = keep_list(index);
        m = length(gamma);
        
        if(m == 0)
            break;
        end
        
    end
    
    
    % ====== Compute new weights ======
    G = repmat(sqrt(gamma)',K,1);
    PhiG = F.*G;
    [U,S,V] = svd(PhiG,'econ');
    
    [d1,~] = size(S);
    if(d1 > 1)
        diag_S = diag(S);
    else
        diag_S = S(1);
    end
    
    U_scaled = U(:,1:min(K,m)).*repmat((diag_S./(diag_S.^2 + sqrt(lambda) + 1e-16))',K,1);
    Xi = G'.*(V*U_scaled');
    
    mu_old = mu;
    mu = Xi*Y;
    
    
    % *** Update hyperparameters ***
    mu2_bar = sum(abs(mu).^2,2);
    gamma = (mu2_bar/N_RX).^(1-p/2);
    
    
    % ========= Check stopping conditions, etc. =========
    count = count + 1;
    if(count >= MAX_ITERS)
        break
    end
    
    if (size(mu) == size(mu_old))
        dmu = max(max(abs(mu_old - mu)));
        if(dmu < EPSILON)
            break
        end
    end
    
end

% expand final solution
X = zeros(N,N_RX);
if(isempty(keep_list))
    return
end

X(keep_list,:) = mu;

end