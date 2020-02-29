function [p_val Sta] = UInd_KCItest(x, y, width)
% function [p_val] = UInd_test(x, y)
% To test if x and y are unconditionally independent with bootstrap (which is
%       the same as in HSIC test) or with the finite-sample Gamma approximation.
% INPUT:
%   X and Y: data matrices of size number_of_samples * dimensionality.
%   width (optional): the kernel width for x and y.
% Output:
%   p_val: the p value obtained by bootstrapping (if the sample size is
%       smaller than 1000) or by Gamma approximation (if the sample size is
%       large).
% Copyright (c) 2010-2011  Kun Zhang, Jonas Peters.
% All rights reserved.  See the file COPYING for license terms.
%
% For details of the method, see K. Zhang, J. Peters, D. Janzing, and B. Schoelkopf,
%       "A kernel-based conditional independence test and application in causal discovery,"
%       In UAI 2011,
%         and
%       A. Gretton, K. Fukumizu, C.-H. Teo, L. Song, B. Schoelkopf and A. Smola, "A kernel
%       Statistical test of independence." In NIPS 21, 2007.

T = length(y); % the sample size

% Controlling parameters
Bootstrap = 1;
Approximate = 0;
Method_kernel_width = 1; % 1: empirical value; 2: median

% Num_eig = floor(T/4); % how many eigenvalues are to be calculated?
if T>1000
    Num_eig = floor(T/2);
else
    Num_eig = T;
end
T_BS = 2000;
lambda = 1E-3; % the regularization paramter
Thresh = 1E-6;
% normalize the data
x = x - repmat(mean(x), T, 1);
x = x * diag(1./std(x));
y = y - repmat(mean(y), T, 1);
y = y * diag(1./std(y));
Cri = []; Sta = []; p_val = []; Cri_appr = []; p_appr = [];

% use empirical kernel width instead of the median
if ~exist('width', 'var')|isempty(width)|width==0
    if T < 200
        width = 0.8;
    elseif T < 1200
        width =0.5;
    else
        width = 0.3;
    end
end
if Method_kernel_width == 1
    theta = 1/(width^2); % I use this parameter to construct kernel matices. Watch out!! width = sqrt(2) sigma  AND theta= 1/(2*sigma^2)
else
    theta = 0;
end
%    width = sqrt(2)*medbw(x, 1000); %use median heuristic for the band width.
%theta = 1/(width^2); % I use this parameter to construct kernel matices. Watch out!! width = sqrt(2) sigma  AND theta= 1/(2*sigma^2)

H =  eye(T) - ones(T,T)/T; % for centering of the data in feature space
% Kx = kernel([x], [x], [theta/size(x,2),1]); Kx = H * Kx * H; %%%%Problem
% Ky = kernel([y], [y], [theta/size(y,2),1]); Ky = H * Ky * H;  %%%%Problem
Kx = kernel([x], [x], [theta * size(x,2),1]); Kx = H * Kx * H; %%%%Problem
Ky = kernel([y], [y], [theta * size(y,2),1]); Ky = H * Ky * H;  %%%%Problem

Sta = trace(Kx * Ky);



Cri = -1;
p_val = -1;
if Bootstrap
    % calculate the eigenvalues that will be used later
    % Due to numerical issues, Kx and Ky may not be symmetric:
    [eig_Kx, eivx] = eigdec((Kx+Kx')/2,Num_eig);
    [eig_Ky, eivy] = eigdec((Ky+Ky')/2,Num_eig);
    % calculate Cri...
    % first calculate the product of the eigenvalues
    eig_prod = stack( (eig_Kx * ones(1,Num_eig)) .* (ones(Num_eig,1) * eig_Ky'));
    II = find(eig_prod > max(eig_prod) * Thresh);
    eig_prod = eig_prod(II); %%% new method
    
    % use mixture of F distributions to generate the Null dstr
    if length(eig_prod) * T < 1E6
        %     f_rand1 = frnd(1,T-2-df, length(eig_prod),T_BS);
        %     Null_dstr = eig_prod'/(T-1) * f_rand1;
        f_rand1 = chi2rnd(1,length(eig_prod),T_BS);
        Null_dstr = eig_prod'/T * f_rand1; %%%%Problem
    else
        % iteratively calcuate the null dstr to save memory
        Null_dstr = zeros(1,T_BS);
        Length = max(floor(1E6/T),100);
        Itmax = floor(length(eig_prod)/Length);
        for iter = 1:Itmax
            %         f_rand1 = frnd(1,T-2-df, Length,T_BS);
            %         Null_dstr = Null_dstr + eig_prod((iter-1)*Length+1:iter*Length)'/(T-1) * f_rand1;
            f_rand1 = chi2rnd(1,Length,T_BS);
            Null_dstr = Null_dstr + eig_prod((iter-1)*Length+1:iter*Length)'/T * f_rand1;
            
        end
        Null_dstr = Null_dstr + eig_prod(Itmax*Length+1:length(eig_prod))'/T *... %%%%Problem
            chi2rnd(1, length(eig_prod) - Itmax*Length,T_BS);
        %         frnd(1,T-2-df, length(eig_prod) - Itmax*Length,T_BS);
    end
    %         % use chi2 to generate the Null dstr:
    %         f_rand2 = chi2rnd(1, length(eig_prod),T_BS);
    %         Null_dstr = eig_prod'/(TT(epoch)-1) * f_rand2;
    sort_Null_dstr = sort(Null_dstr);
    p_val = sum(Null_dstr>Sta)/T_BS;
end

if Approximate
    mean_appr = trace(Kx) * trace(Ky) /T;
    var_appr = 2* trace(Kx*Kx) * trace(Ky*Ky)/T^2;
    k_appr = mean_appr^2/var_appr;
    theta_appr = var_appr/mean_appr;
    p_appr = 1-gamcdf(Sta, k_appr, theta_appr);
    p_val = p_appr;
end
