% Incomplete Cholesky factorization with RBF kernel
%
% Description:
% Finds low rank approximation of RBF kernel Gram matrix K = PGG'P for the
% n x p data matrix X. Here, K is an n x n Gram matrix, G is n x m with m << n, 
% and P is a permutation matrix.
%
% Arguments:
% X       n x p data matrix 
% sigma   bandwidth for RBF kernel
% tol     threshold for remaining eigenvalues to consider
%
% Output:
% G       n x m matrix (m << n)
% P       n vector of permutation indices    
%
% 
% Adapted from Francis Bach's Cholesky with side information implementation
%
% Copyright (c) 2010  Robert Tillman  [rtillman@cmu.edu]
%               2005  Francis Bach    [francis.bach@ens.fr]
% All rights reserved.  See the file COPYING for license terms.
function [G,P] = inchol(X,sigma,tol)

if (sigma<=0)
   error('sigma must be > 0');
end
if (tol<=0)
   error('tol must be > 0');
end

n = size(X,1);
% begin with full matrix
G = zeros(n,n);
% using RBF kernel so diagonal entries are ones
diagK = ones(n,1);
% permutation indices;
P = 1:n;
% updated diagonal elements
D = diagK;

% construct columns of K until threshold is met
for k=1:n

    % select next most informative pivot
    best = D(k);
    bestInd = k;
    for j=k:n
        if (D(j) > best/.99)
            best = D(j);
            bestInd = j;
        end
    end

    % threshold met so remove columns to the right and break
    if best<tol
        G = G(:,1:(k-1));
        break;
    end

    % move pivot to the front
    pk = P(k);
    P(k)=P(bestInd);
    P(bestInd) = pk;
    dk = D(k);
    D(k)=D(bestInd);
    D(bestInd)=dk;

    % update existing columns	
    for j=1:(k-1)
       gk = G(k,j);
       G(k,j)=G(bestInd,j);
       G(bestInd,j)=gk;
    end  

    % compute new Cholesky column
    G(k,k)=sqrt(D(k));
    G(k+1:n,k)=1/G(k,k)*(rbf(X(P(k+1:n),:),repmat(X(P(k),:),(n-k),1),sigma) - G(k+1:n,1:k-1)*(G(k,1:k-1))');
    
    % update diagonal
    D(k+1:n) =  D(k+1:n) - G(k+1:n,k).^2;
    D(k) = 0;

end
