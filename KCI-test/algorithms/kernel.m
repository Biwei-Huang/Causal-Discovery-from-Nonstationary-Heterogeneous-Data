function [kx, bw_new] = kernel(x, xKern, theta)

% KERNEL Compute the rbf kernel
% Copyright (c) 2010-2011  ...
% All rights reserved.  See the file COPYING for license terms.
n2 = dist2(x, xKern);
if theta(1)==0
    theta(1)=2/median(n2(tril(n2)>0));
    theta_new=theta(1);
end
wi2 = theta(1)/2;
kx = theta(2)*exp(-n2*wi2);
bw_new=1/theta(1);
   
