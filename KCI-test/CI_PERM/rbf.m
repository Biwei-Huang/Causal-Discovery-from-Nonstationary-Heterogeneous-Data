% RBF kernel evaluation
%
% Description:
% Evaluates RBF kernel for n points
%
% Input:
% x1      n x p matrix (n points with dimensionality p)
% x2      n x p matrix (n points with dimensionality p)
% sigma   kernel bandwidth
%
% Output:
% k       n x 1 matrix of k(x1,x2) evaluations
%
% Copyright (c) 2010  Robert Tillman  [rtillman@cmu.edu]
% All rights reserved.  See the file COPYING for license terms.

function k = rbf(x1,x2,sigma)

if (size(x1,1)~=size(x2,1))
  error('x1 and x2 must contain the same number of data points');
end
if (size(x1,2)~=size(x2,2)) 
  error('x1 and x2 must be of the same dimensionality');
end
if (sigma<=0) 
  error('sig must be > 0');
end

k = exp(-.5*sum((x1-x2).^2,2)/(sigma^2));
