% Median distance heuristic for setting RBF kernel bandwidth
%
% Description:
% Uses the median distance between points for setting the bandwidth for RBF kernels.
%
% Arguments:
% X             n x p matrix of n datapoints with dimensionality p
% maxpoints     maximum number of points to use when setting the bandwidth
%
% Output:
% sigma         value for bandwidth
%
%
% Copyright (c) 2010  Robert Tillman  [rtillman@cmu.edu]
%               2007  Arthur Gretton  [arthur.gretton@tuebingen.mpg.de]
% All rights reserved.  See the file COPYING for license terms.

function sigma = medbw(X, maxpoints)

if (maxpoints < 1 || maxpoints ~= int32(maxpoints))
  error('maxpoints must be a positive integer')
end

% truncates data if more points than maxpoints
n =  size(X,1);
if (n>maxpoints)
   med = X(1:maxpoints,:);
   n = maxpoints;
else
   med = X;
end

% finds median distance between points
G = sum((med.*med),2);
Q = repmat(G,1,n);
R = repmat(G',n,1);
dists = Q + R - 2*med*med';
dists = dists-tril(dists);
dists=reshape(dists,n^2,1);
sigma = sqrt(0.5*median(dists(dists>0)));
