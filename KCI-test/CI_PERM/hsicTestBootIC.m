%This function implements the HSIC independence test using a bootstrap approximation
%to the test threshold

%Inputs: 
%        X contains dx columns, m rows. Each row is an i.i.d sample
%        Y contains dy columns, m rows. Each row is an i.i.d sample
%        alpha is the level of the test
%        shuffles is number of shuffles to approximate null distribution

%Outputs: 
%        sig: boolean indicator of whether the test was significant
%        p: p-value

%Set kernel size to median distance between points, if no kernel specified

% Copyright (c) 2010  Robert Tillman  [rtillman@cmu.edu]
%               2007  Arthur Gretton  [arthur.gretton@tuebingen.mpg.de]
% All rights reserved.  See the file COPYING for license terms.

function [sig,p] = hsicTestBootIC(X,Y,alpha,shuffles);

% tolerance for incomplete Cholesky
tol = 1e-12;

m=size(X,1);

% set kernel size to median distance between points
maxpoints = 1000;
sigx = medbw(X, maxpoints);
sigy = medbw(Y, maxpoints);

%Compute the approximations of Gram matrices
[K, Pk] = inchol(X,sigx,tol);
[L, Pl] = inchol(Y,sigy,tol);

%bone = ones(m,1);
%H = eye(m)-1/m*ones(m,m);

% center Gram matrices and permute indices
Kc = K(Pk,:) - repmat((sum(K)/m),m,1);
Lc = L(Pl,:) - repmat((sum(L)/m),m,1);
    
testStat = (1/m^2)*sum(sum(Kc.*((Kc'*Lc)*Lc')'));

  HSICarr = zeros(shuffles,1);
  for whichSh=1:shuffles
   
    [notUsed,indL] = sort(rand(m,1));
    
    newLc = Lc(indL,:);
    HSICarr(whichSh) = (1/m^2)*sum(sum(Kc.*((Kc'*newLc)*newLc')'));
    
  end 

% get p-value from empirical cdf
p = length(find(HSICarr>=testStat))/shuffles;

% determine significance
sig=(p<=alpha);
