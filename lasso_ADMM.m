function x = lasso_ADMM(A, y, lambda)
%
% ADMM-based algorithm for solving the LASSO problem: y = Ax + w
%
% INPUTS:
%   A       : sensing matrix
%   y       : measurements vector
%   lambda	: regularization parameter
%
% OUTPUTS:
%   x       : estimated solution vector
%
% REFERENCE:
%   [1] Boyd, Stephen, Neal Parikh, and Eric Chu.
%       Distributed optimization and statistical learning via the alternating
%       direction method of multipliers. Now Publishers Inc, 2011.

%%
[~,M] = size(A);
u = zeros(M,1);
z = zeros(M,1);

rho = 1;            % penalty parameter
iterMax = 500;      % max number of iterations
haltThres = 1e-5;   % halting threshold

invA = inv(A'*A + rho*eye(M)); % cache
Ay = A'*y; % cache

xold = 0;
for iter=1:iterMax
    
    x = invA*(Ay + rho*(z - u));
    
    z = soft(x+u,lambda/rho);
    
    u = u + x - z;
    
    if(norm(x - xold) < haltThres)
        break;
    end
    
    xold = x;
    
end

end

%% soft thresholding
function y = soft(x,lambda)

y = sign(x).*max(abs(x) - lambda,0);

end