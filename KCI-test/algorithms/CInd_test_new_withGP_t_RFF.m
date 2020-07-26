function [p_val, Sta, Cri] = CInd_test_new_withGP_t_RFF(x, y, z, alpha, pars)
% To test if x and y are independent.
% INPUT:
%   The number of rows of x and y is the sample size.
%   alpha is the significance level (we suggest 1%).
%   pars contains the kernel width and whether to use GP to optimize the kernel width.
% Output:
%   Cri: the critical point at the p-value equal to alpha obtained by bootstrapping.
%   Sta: the statistic Tr(K_{\ddot{X}|Z} * K_{Y|Z}).
%   p_val: the p value obtained by bootstrapping.
% If Sta > Cri, the null hypothesis (x is independent from y) is rejected.
% Copyright (c) 2010-2011  ...
% All rights reserved.  See the file COPYING for license terms.

% Controlling parameters
width = pars.width;
if(pars.widthT==0) % kernel width on time index when IF_GP=0, need tunning!!!!!
    widthT = 0.1;
else
    widthT = pars.widthT;
end
IF_unbiased = 0;
Bootstrap = 1;

T = length(y); % the sample size
% Num_eig = floor(T/4); % how many eigenvalues are to be calculated?
Num_eig = T;
T_BS = 10000; % 5000
lambda = 1E-3; % the regularization paramter  %%%%Problem
Thresh = 1E-5;
% normalize the data
x = x - repmat(mean(x), T, 1);
x = x * diag(1./std(x));
y = y - repmat(mean(y), T, 1);
y = y * diag(1./std(y));
z = z - repmat(mean(z), T, 1);
z = z * diag(1./std(z));

D = size(z, 2);
logtheta_x = []; logtheta_y = [];  df_x = []; df_y = [];
Cri = []; Sta = []; p_val = []; Cri_appr = []; p_appr = [];

if width ==0
    if T <= 200
        width = 1.2; % 0.8
    elseif T < 1200
        width = 0.6;
    else
        width = 0.4; % 0.3
    end
end

Zx = transformFeatures([x z/2]/(width*sqrt(D)/sqrt(2))); % calculate random fourier featuress
Zx = Zx - repmat(mean(Zx,2), 1, T);
Zy = transformFeatures(y/(width*sqrt(D)/sqrt(2)));
Zy = Zy - repmat(mean(Zy,2), 1, T);
% check whether the last dimension of z is the time index
tmp = [1:T]';
tmp = tmp - repmat(mean(tmp), T, 1);
tmp = tmp * diag(1./std(tmp));
if(norm(z(:,end)-tmp)<1e-5)
    if(D>1)
        Zz = transformFeatures([z(:,1:end-1)/(width*sqrt(D-1)/sqrt(2)),z(:,end)/(widthT/sqrt(2))]);
    else
        Zz = transformFeatures(z/(widthT/sqrt(2)));
    end
else
    Zz = transformFeatures(z/(width*sqrt(D)/sqrt(2)));
end
Zz = Zz - repmat(mean(Zz,2), 1, T);


Cxz = Zx*Zz';
Cyz = Zy*Zz';
Czz = Zz*Zz';

Ex = Zx - (Cxz/(Czz+1e-5*eye(size(Czz,1))))*Zz;
Ex = Ex - repmat(mean(Ex,2), 1, T);
Ey = Zy - (Cyz/(Czz+1e-5*eye(size(Czz,1))))*Zz;
Ey = Ey - repmat(mean(Ey,2), 1, T);
C = Ex * Ey';
% calculate the statistic
Sta = norm(C,'fro')^2;

% calculate the eigenvalues that will be used later
[eig_Kxz, eivx] = eigdec(Ex'*Ex,Num_eig);
[eig_Kyz, eivy] = eigdec(Ey'*Ey,Num_eig);

% calculate the product of the square root of the eigvector and the eigenvector
IIx = find(eig_Kxz > max(eig_Kxz) * Thresh);
IIy = find(eig_Kyz > max(eig_Kyz) * Thresh);
eig_Kxz = eig_Kxz(IIx);
eivx = eivx(:,IIx);
eig_Kyz = eig_Kyz(IIy);
eivy = eivy(:,IIy);

eiv_prodx = eivx * diag(sqrt(eig_Kxz));
eiv_prody = eivy * diag(sqrt(eig_Kyz));
clear eivx eig_Kxz eivy eig_Kyz
% calculate their product
Num_eigx = size(eiv_prodx, 2);
Num_eigy = size(eiv_prody, 2);
Size_u = Num_eigx * Num_eigy;
uu = zeros(T, Size_u);
for i=1:Num_eigx
    for j=1:Num_eigy
        uu(:,(i-1)*Num_eigy + j) = eiv_prodx(:,i) .* eiv_prody(:,j);
    end
end
if Size_u > T
    uu_prod = uu * uu';
else
    uu_prod = uu' * uu;
end
if Bootstrap
    eig_uu = eigdec(uu_prod,min(T,Size_u));
    II_f = find(eig_uu > max(eig_uu) * Thresh);
    eig_uu = eig_uu(II_f);
end

Cri=-1;
p_val=-1;


if Bootstrap
    % use mixture of F distributions to generate the Null dstr
    if length(eig_uu) * T < 1E6
        f_rand1 = chi2rnd(1,length(eig_uu),T_BS);
        if IF_unbiased
            Null_dstr = T^2/(T-1-df_x)/(T-1-df_y) * eig_uu' * f_rand1; %%%%Problem
        else
            Null_dstr = eig_uu' * f_rand1;
        end
    else
        % iteratively calcuate the null dstr to save memory
        Null_dstr = zeros(1,T_BS);
        Length = max(floor(1E6/T),100);
        Itmax = floor(length(eig_uu)/Length);
        for iter = 1:Itmax
            f_rand1 = chi2rnd(1,Length,T_BS);
            if IF_unbiased
                Null_dstr = Null_dstr + T^2/(T-1-df_x)/(T-1-df_y) *... %%%%Problem
                    eig_uu((iter-1)*Length+1:iter*Length)' * f_rand1;
            else
                Null_dstr = Null_dstr + ... %%%%Problem
                    eig_uu((iter-1)*Length+1:iter*Length)' * f_rand1;
            end
            
        end
    end
    sort_Null_dstr = sort(Null_dstr);
    Cri = sort_Null_dstr(ceil((1-alpha)*T_BS));
    p_val = sum(Null_dstr>Sta)/T_BS;
end



