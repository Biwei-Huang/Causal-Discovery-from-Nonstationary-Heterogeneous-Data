function [pval, stat] = indtest_hsic(X, Y, Z, pars)
% function [pval, stat] = indtest_hsic(X, Y, Z, pars)
%
% This function is a WRAPPER
% Performs either an HSIC test (Gretton et al.) or a conditional HSIC test (Fukumizu et. al)
%
% INPUT:
%   X          Nxd1 matrix of samples (N data points, d1 dimensions)
%   Y          Nxd2 matrix of samples (N data points, d2 dimensions)
%   Z          Nxd3 matrix of samples (N data points, d3 dimensions)
%   pars       structure containing parameters for the independence test
%     .pairwise	    if true, the test is performed pairwise if d1>1 (standard: false)
%     .bonferroni   if true, bonferroni correction is performed (standard: false)
%     .perm         # of bootstrap samples for cond. hsic test (standard: 500)
%
% OUTPUT:
%   pval      p value of the test
%   stat      test statistic
%
%
% Copyright (c) 2011-2011  Kun Zhang
%               2011-2011  Jonas Peters
% All rights reserved.  See the file COPYING for license terms.


if ~isfield(pars,'pairwise')
    pars.pairwise = false;
end;

if ~isfield(pars,'bonferroni')
    pars.bonferroni = false;
end;

if ~isfield(pars,'perm')
    pars.perm= 500;
end;

if isempty(Z) %unconditional HSIC
    if pars.pairwise
        p = zeros(size(X,2),size(Y,2));
        for i = 1:size(X,2);
            for j = 1:size(Y,2);
                [p(i,j) sta(i,j)] = fasthsic(X(:,i),Y(:,j));
            end
        end
        [pp iii] = min(p);
        [pval jj] = min(pp);
        stat = sta(iii(jj),jj);
        if pars.bonferroni
	        pval=size(X,2)*size(Y,2)*pval;
        end
    else
        [pval stat]= fasthsic(X, Y);
    end
else %conditional HSIC
    [aa, pval, stat]=hsiccondTestIC(X,Y,Z,0.8,pars.perm);
end

return



