function [p_val, Sta, Cri] = CInd_test_new_withGP_t(x, y, z, alpha, pars)
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
IF_GP = pars.if_GP1; % Set IF_GP=0 to speed up the process
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
T_BS = 5000; 
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
    %      width = sqrt(2)*medbw(x, 1000); %use median heuristic for the band width.
end

theta = 1/(width^2 * D); % I use this parameter to construct kernel matices. Watch out!! width = sqrt(2) sigma  AND theta= 1/(2*sigma^2)

H =  eye(T) - ones(T,T)/T; % for centering of the data in feature space
% Kx = kernel([x z], [x z], [theta,1]); Kx = H * Kx * H;
Kx = kernel([x z/2], [x z/2], [theta,1]); Kx = H * Kx * H;
% Ky = kernel([y z], [y z], [theta,1]); %Ky = Ky * H;
% Kx = kernel([x], [x], [theta,1]); %Kx = Kx * H; %%%%Problem
Ky = kernel([y], [y], [theta,1]); Ky = H * Ky * H;  %%%%Problem
% later with optimized hyperparameters

if IF_GP
    % learning the hyperparameters
    [eig_Kx, eix] = eigdec((Kx+Kx')/2, min(400, floor(T/4))); % /2
    [eig_Ky, eiy] = eigdec((Ky+Ky')/2, min(200, floor(T/5))); % /3
    % disp('  covfunc = {''covSum'', {''covSEard'',''covNoise''}};')
    covfunc = {'covSum', {'covSEard','covNoise'}};
    logtheta0 = [log(width * sqrt(D))*ones(D,1); 0; log(sqrt(0.1))];
    fprintf('Optimizing hyperparameters in GP regression...\n');
    
    %old gpml-toolbox
    %
    IIx = find(eig_Kx > max(eig_Kx) * Thresh); eig_Kx = eig_Kx(IIx); eix = eix(:,IIx);
    IIy = find(eig_Ky > max(eig_Ky) * Thresh); eig_Ky = eig_Ky(IIy); eiy = eiy(:,IIy);
    [logtheta_x, fvals_x, iter_x] = minimize(logtheta0, 'gpr_multi', -350, covfunc, z, 2*sqrt(T) *eix * diag(sqrt(eig_Kx))/sqrt(eig_Kx(1)));
    [logtheta_y, fvals_y, iter_y] = minimize(logtheta0, 'gpr_multi', -350, covfunc, z, 2*sqrt(T) *eiy * diag(sqrt(eig_Ky))/sqrt(eig_Ky(1)));
    
    logtheta_x
    logtheta_y
    
    
    covfunc_z = {'covSEard'};
    Kz_x = feval(covfunc_z{:}, logtheta_x, z);
    Kz_y = feval(covfunc_z{:}, logtheta_y, z);
    
    % Note: in the conditional case, no need to do centering, as the regression
    % will automatically enforce that.
    
    % Kernel matrices of the errors
    P1_x = (eye(T) - Kz_x*pdinv(Kz_x + exp(2*logtheta_x(end))*eye(T)));
    Kxz = P1_x* Kx * P1_x';
    P1_y = (eye(T) - Kz_y*pdinv(Kz_y + exp(2*logtheta_y(end))*eye(T)));
    Kyz = P1_y* Ky * P1_y';
    % calculate the statistic
    Sta = trace(Kxz * Kyz);
    
    % degrees of freedom
    df_x = trace(eye(T)-P1_x);
    df_y = trace(eye(T)-P1_y);
else   
    % check whether the last dimension of z is the time index
    tmp = [1:T]';
    tmp = tmp - repmat(mean(tmp), T, 1);
    tmp = tmp * diag(1./std(tmp));
    if(norm(z(:,end)-tmp)<1e-5)
        covfunc_z = {'covSEard'};
        logtheta = [log(width * sqrt(D-1))*ones(D-1,1); log(widthT) ; 0];
        Kz = feval(covfunc_z{:}, logtheta, z);
    else
        Kz = kernel(z, z, [theta,1]);
    end
    Kz = H * Kz * H;
    % Kernel matrices of the errors
    P1 = (eye(T) - Kz*pdinv(Kz + lambda*eye(T)));
    Kxz = P1* Kx * P1';
    Kyz = P1* Ky * P1';
    % calculate the statistic
    Sta = trace(Kxz * Kyz);
    
    % degrees of freedom
    df = trace(eye(T)-P1);
end

% calculate the eigenvalues
% Due to numerical issues, Kxz and Kyz may not be symmetric:
[eig_Kxz, eivx] = eigdec((Kxz+Kxz')/2,Num_eig);
[eig_Kyz, eivy] = eigdec((Kyz+Kyz')/2,Num_eig);

% calculate the product of the square root of the eigvector and the eigen
% vector
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

