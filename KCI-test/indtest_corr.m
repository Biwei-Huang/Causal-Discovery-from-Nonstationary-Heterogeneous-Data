function [pval, stat] = indtest_corr(X, Y, Z, pars)
% function [pval, stat] = indtest_corr(X, Y, Z, pars)
%
% Uses: statistics toolbox matlab
%
% Performs either a correlation test or a partial correlation test
%
% INPUT:
%   X           Nxd1 matrix of samples (N data points, d1 dimensions)
%   Y           Nxd2 matrix of samples (N data points, d2 dimensions)
%   Z           Nxd3 matrix of samples (N data points, d3 dimensions)
%   pars        structure containing parameters for the independence test
%     .bonferroni   if true, bonferroni correction is performed (standard: false)
%
% OUTPUT:
%   pval      p value of the test
%   stat      test statistic
%
%
% Copyright (c) 2011-2011  Kun Zhang
%               2011-2011  Jonas Peters
% All rights reserved.  See the file COPYING for license terms.


if ~isfield(pars,'bonferroni')
    pars.bonferroni = false;
end;

if isempty(Z)
    [sta,pp] = corr(X,Y);
    stat=max(sta);
    pval = min(min(pp));
    if pars.bonferroni
	pval=size(X,2)*size(Y,2)*pval;
    end
else
    [sta, pp]=partialcorr(X,Y,Z);
    pval = min(min(pp));
    stat=max(sta);
    if pars.bonferroni
	pval=size(X,2)*size(Y,2)*pval;
    end
end

return
