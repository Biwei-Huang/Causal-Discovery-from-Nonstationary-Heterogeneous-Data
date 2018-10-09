function y = logdet(A)
% log(det(A)) where A is positive-definite.
% This is faster and more stable than using log(det(A)).
% Copyright (c) 2010-2011  ...
% All rights reserved.  See the file COPYING for license terms.
try
  U = chol(A);
  y = 2*sum(log(diag(U)));
catch
  [void, errid] = lasterr;
  if strcmp(errid, 'MATLAB:posdef')
    warning(['Matrix is not positive definite in logdet, using log(det())'])
    y = log(det(A));
    return
  else
    error(lasterr)
  end
end
