% Conditional dependence operator empirical estimator with incomplete Choleksy 
% factorization for low rank approximation of Gram matrices
%
% Arguments:
% Gx        low rank approximation for centered Gram matrix for X
% Gy        low rank approximation for centered Gram matrix for Y
% Gz        low rank approximation for centered Gram matrix for Z
% epsilon   the smoothing constant
%
% Output: 
% emphsic   the test statistic
%
% Copyright (c) 2010  Robert Tillman  [rtillman@cmu.edu]
% All rights reserved.  See the file COPYING for license terms. 

function emphsic = hsiccondIC(Gx,Gy,Gz,epsilon)

n = size(Gx,1);
if (n~=size(Gy,1) || n~=size(Gz,1))
   error('Gx, Gy, and Gz must have the same number of rows');
end
if (epsilon<=0)
   error('epsilon must > 0');
end

mx = size(Gx,2);
my = size(Gy,2);
mz = size(Gz,2);

[Ux, Sx, Vx] = svd(Gx,'econ');
[Uy, Sy, Vy] = svd(Gy,'econ');
[Uz, Sz, Vz] = svd(Gz,'econ');

Sxsq = diag(Sx).^2;
Sysq = diag(Sy).^2;
Szsq = diag(Sz).^2;
Szsqe = Szsq + epsilon;
Szsqt = Szsq./Szsqe;

% first term  - GxGx'GyGy'
first = sum(sum((Ux*(diag(Sxsq)*(Ux'*Uy)*diag(Sysq))).*Uy));

% second term - 2GyGy'GzGz'(GzGz' - epsilonI)^(-2)GzGz'GxGx'
second1 = Ux*(diag(Sxsq)*(Ux'*Uz)*diag(Szsqt)*(Uz'*Uy)*diag(Sysq));
second = -2*sum(sum(second1.*Uy));

% third term  - 2GyGy'GzGz'(GzGz' - epsilonI)^(-2)GzGz'GxGx'GzGz'(GzGz' - epsilonI)^(-2)GzGz'
third = sum(sum((second1*(Uy'*Uz)*diag(Szsqt)).*Uz));

% compute test statistic using first, second, and third terms above with
% the U-statistic
emphsic = (first+second+third)/((n-1)^2);
