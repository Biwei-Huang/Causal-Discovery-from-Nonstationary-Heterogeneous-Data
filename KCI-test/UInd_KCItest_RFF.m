function [p_val Sta] = UInd_KCItest_RFF(x, y, pars)
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

T = length(y); % the sample size

% Controlling parameters
width = pars.width;
Bootstrap = 1;

Method_kernel_width = 1; % 1: empirical value; 2: median

% Num_eig = floor(T/4); % how many eigenvalues are to be calculated?
if T>1000
    Num_eig = floor(T/2);
else
    Num_eig = T;
end
T_BS = 5000;
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


Zx = transformFeatures(x/width); % calculate random fourier featuress
Zy = transformFeatures(y/width); 
C = (Zx - repmat(mean(Zx,2), 1, T)) * (Zy - repmat(mean(Zy,2), 1, T))';
Sta = norm(C,'fro')^2;


Cri = -1;
p_val = -1;
if Bootstrap
    % calculate the eigenvalues that will be used later
    [eig_Kx, eivx] = eigdec(Zx'*Zx,Num_eig);
    [eig_Ky, eivy] = eigdec(Zy'*Zy,Num_eig);
    % calculate Cri...
    % first calculate the product of the eigenvalues
    eig_prod = stack( (eig_Kx * ones(1,length(eig_Kx))) .* (ones(length(eig_Kx),1) * eig_Ky'));
    II = find(eig_prod > max(eig_prod) * Thresh);
    eig_prod = eig_prod(II); %%% new method
    
    % use mixture of F distributions to generate the Null dstr
    if length(eig_prod) * T < 1E6
        f_rand1 = chi2rnd(1,length(eig_prod),T_BS);
        Null_dstr = eig_prod'/T * f_rand1; %%%%Problem
    else
        % iteratively calcuate the null dstr to save memory
        Null_dstr = zeros(1,T_BS);
        Length = max(floor(1E6/T),100);
        Itmax = floor(length(eig_prod)/Length);
        for iter = 1:Itmax
            f_rand1 = chi2rnd(1,Length,T_BS);
            Null_dstr = Null_dstr + eig_prod((iter-1)*Length+1:iter*Length)'/T * f_rand1;
            
        end
        Null_dstr = Null_dstr + eig_prod(Itmax*Length+1:length(eig_prod))'/T *... %%%%Problem
            chi2rnd(1, length(eig_prod) - Itmax*Length,T_BS);
    end
    sort_Null_dstr = sort(Null_dstr);
    p_val = sum(Null_dstr>Sta)/T_BS;
end

