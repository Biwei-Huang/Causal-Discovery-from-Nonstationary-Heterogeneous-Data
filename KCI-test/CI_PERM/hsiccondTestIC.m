% Statistical test for kernel conditional independence of X and Y given Z with 
% incomplete Cholesky factorization for low rank approximation of Gram matrices
%
% Arguments:
% X          n x p vector of data points
% Y          n x m vector of data points
% Z          n x r vector of data points
% alpha      significance level
% shuffles   number of shuffles for the permutation test
%
% Output:
% sig        boolean indicator of whether the test was significant for the given alpha
% p          resulting p-value
%
% Copyright (c) 2010  Robert Tillman  [rtillman@cmu.edu]
% All rights reserved.  See the file COPYING for license terms.

function [sig,p,testStat] = hsiccondTestIC(X,Y,Z,alpha,shuffles)

n = size(X,1);
if (n~=size(Y,1) || n~=size(Z,1))
   error('X, Y, and Z must have the same number of data points');
end
if (alpha<0 || alpha>1)
   error('alpha must be between 0 and 1');
end
if (shuffles<=0 || shuffles~=int32(shuffles))
   error('number of shuffles must be a positive integer');
end

% smoothing constant for conditional cross covariance operator
epsilon=1e-4;
% threshold for eigenvalues to consider in low rank Gram matrix approximations
tol = 1e-4;

% augment X and Y for conditional test
X = [X,Z];
Y = [Y,Z];

% set kernel size to median distance between points
maxpoints = 1000;
sigx = medbw(X, maxpoints);
sigy = medbw(Y, maxpoints);
sigz = medbw(Z, maxpoints);

% low rank approximation of Gram matrices using incomplete Cholesky factorization
[K, Pk] = inchol(X,sigx,tol);
[L, Pl] = inchol(Y,sigy,tol);
[M, Pm] = inchol(Z,sigz,tol);

% center Gram matrices factoring in permutations made during low rank approximation
Kc = K(Pk,:) - repmat((sum(K)/n),n,1);
Lc = L(Pl,:) - repmat((sum(L)/n),n,1);
Mc = M(Pm,:) - repmat((sum(M)/n),n,1);

% compute the U-statistic
%pairs = nchoosek(1:n,2);
%bz = n*(n-1)/sum(rbf(Z(pairs(:,1)),Z(pairs(:,2)),sigz).^2);

% compute HSIC dependence value
testStat = hsiccondIC(Kc,Lc,Mc,epsilon);

% first cluster Z;
nc = pickK(Z);
clusters = kmeans(Z,nc,'EmptyAction','drop','MaxIter',1000,'Display','off');
%[centers,clusters,datapoints] = MeanShiftCluster(Z,sigz,false);
%nc = length(centers);

% simulate null distribution and permutation test
nullapprox = zeros(shuffles,1);
for i=1:shuffles
    % permute within clusters
    Plnew = 1:n;
    for j=1:nc
        indj = find(clusters==j);
        pj = indj(randperm(length(indj)));
        Plnew(indj) = Plnew(pj);
    end
    % centered Gram matrix for new sample
    newLc = Lc(Plnew,:);
    % compute HSIC dependence value for new sample
    nullapprox(i)=hsiccondIC(Kc,newLc,Mc,epsilon);
end

% get p-value from empirical cdf
p = length(find(nullapprox>=testStat))/shuffles;

% determine significance
sig=(p<=alpha);
