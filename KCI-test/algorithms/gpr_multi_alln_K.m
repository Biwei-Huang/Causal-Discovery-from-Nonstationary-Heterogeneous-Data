function [out1, out2] = gpr_multi(logtheta, covfunc, x, Ky, xstar);
% Here we change the function gpr to gpr_multi, in which y contains a set
% of vectors on which we do repression from x

% gpr - Gaussian process regression, with a named covariance function. Two
% modes are possible: training and prediction: if no test data are given, the
% function returns minus the log likelihood and its partial derivatives with
% respect to the hyperparameters; this mode is used to fit the hyperparameters.
% If test data are given, then (marginal) Gaussian predictions are computed,
% whose mean and variance are returned. Note that in cases where the covariance
% function has noise contributions, the variance returned in S2 is for noisy
% test targets; if you want the variance of the noise-free latent function, you
% must substract the noise variance.
%
% usage: [nlml dnlml] = gpr(logtheta, covfunc, x, y)
%    or: [mu S2]  = gpr(logtheta, covfunc, x, y, xstar)
%
% where:
%
%   logtheta is a (column) vector of log hyperparameters
%   covfunc  is the covariance function
%   x        is a n by D matrix of training inputs
%   y        is a (column) vector (of size n) of targets
%   xstar    is a nn by D matrix of test inputs
%   nlml     is the returned value of the negative log marginal likelihood
%   dnlml    is a (column) vector of partial derivatives of the negative
%                 log marginal likelihood wrt each log hyperparameter
%   mu       is a (column) vector (of size nn) of prediced means
%   S2       is a (column) vector (of size nn) of predicted variances
%
% For more help on covariance functions, see "help covFunctions".
%
% (C) copyright 2006 by Carl Edward Rasmussen (2006-03-20).

if ischar(covfunc), covfunc = cellstr(covfunc); end % convert to cell if needed
[n, D] = size(x);
[n, m] = size(Ky);
if eval(feval(covfunc{:})) ~= size(logtheta, 1)
  error('Error: Number of parameters do not agree with covariance function')
end

K = feval(covfunc{:}, logtheta, x);    % compute training set covariance matrix

L = chol(K)';                        % cholesky factorization of the covariance
% for i = 1:m
%     alpha(:,i) = solve_chol(L',y(:,i));
% end
% alpha = solve_chol(L',y);
K_inv = solve_chol(L',eye(n));

if nargin == 4 % if no test cases, compute the negative log marginal likelihood

%   out1 = 0.5* trace(y'*alpha) + n*sum(log(diag(L))) + 0.5*n*n*log(2*pi);
  out1 = 0.5* trace(K_inv * Ky) + n*sum(log(diag(L))) + 0.5*n*n*log(2*pi);

  if nargout == 2               % ... and if requested, its partial derivatives
    out2 = zeros(size(logtheta));       % set the size of the derivative vector
    W = K_inv * (n *eye(n) - Ky * K_inv);                % precompute for convenience
    for i = 1:length(out2)
      out2(i) = sum(sum(W.*feval(covfunc{:}, logtheta, x, i)))/2;
    end
  end

else                    % ... otherwise compute (marginal) test predictions ...

  [Kss, Kstar] = feval(covfunc{:}, logtheta, x, xstar);     %  test covariances

  out1 = Kstar' * alpha;                                      % predicted means

  if nargout == 2
    v = L\Kstar;
    out2 = Kss - sum(v.*v)';
  end  

end
