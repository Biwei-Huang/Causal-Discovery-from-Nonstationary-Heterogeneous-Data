function [testStat] = infer_nonsta_dir(X,Y,width,Wt)
% Copyright (C) 2017 Kun Zhang

% infer the causal direction between X and Y when their causal modules are
% both nonstationary but independent
% % width: the kernel width for Y; Wt is the initial kernel width on C (or
% T).
% X: parents; Y; effect; and you can send width to 1, and set Wt to 50
% Don't forget to normalize the data

% size of X:
GP_hyp = 1;
[T,d] = size(X);
X = X * diag(1./std(X));
Y = Y * diag(1./std(Y));
% Wt = 50; % 80
% width = 1; % 0.7
theta = 1/width^2; % 0.2
lambda = 2; % 0.05 0.3  10
Ml = [];

% size of Y should be T*1.
Kyy = kernel(Y, Y, [theta,1]);

%% P(Y|X)
if GP_hyp
    Thresh = 1E-4;
    [eig_Ky, eiy] = eigdec((Kyy+Kyy')/2, min(400, floor(T/4))); % /2
    covfunc = {'covSum', {'covSEard','covNoise'}}; 
    logtheta0 = [log(width)*ones(d,1);  log(Wt); 0; log(sqrt(0.1))];
    fprintf('Optimization hyperparameters in GP regression:\n');
    
    IIy = find(eig_Ky > max(eig_Ky) * Thresh); eig_Ky = eig_Ky(IIy); eiy = eiy(:,IIy);
    [logtheta_y, fvals_y, iter_y] = minimize(logtheta0, 'gpr_multi', -350, covfunc, [X (1:T)'], 2*sqrt(T) *eiy * diag(sqrt(eig_Ky))/sqrt(eig_Ky(1)));
%     exp(logtheta_y),
    if(logtheta_y(d+1)>log(1e4))
        logtheta_y(d+1)=log(1e4);
    end
    
    covfunc_z = {'covSEard'};
    Kxt = feval(covfunc_z{:}, logtheta_y, [X (1:T)']);
    % Note: in the conditional case, no need to do centering, as the regression
    % will automatically enforce that.
    
    % Kernel matrices of the errors
    invK = pdinv(Kxt + exp(2*logtheta_y(end))*eye(T));
    
%     Kxx = kernel(X, X, [1/exp(2*logtheta_y(1:d)),1]);
%     Ktt = kernel((1:T)', (1:T)', [1/exp(2*logtheta_y(d+1)),1]);
    Kxx = feval(covfunc_z{:}, logtheta_y([1:d,d+2]), X);
    Ktt = feval(covfunc_z{:}, logtheta_y([d+1,d+2]), (1:T)');
else
    Kxx = kernel(X, X, [theta,1]);
    Kyy = kernel(Y, Y, [theta,1]);
    Ktt = kernel((1:T)', (1:T)', [1/Wt^2,1]);
    invK = pdinv( Kxx.* Ktt +  lambda * eye(T));
end
Kxx3 = Kxx^2; %^3
prod_invK =  invK * Kyy * invK;
% now finding Ml
Ml = 1/T^2 * Ktt*( Kxx3 .* prod_invK) * Ktt;
% the square distance
D = diag(diag(Ml)) * ones(size(Ml)) + ones(size(Ml)) * diag(diag(Ml)) - 2*Ml;
% Gaussian kernel
sigma2_square = median( D(find(tril(ones(size(D)),-1))) );
Mg = exp(-D/sigma2_square/2);


%% P(X)
Kxx = kernel(X,X, [theta,1]);
[eig_Kx, eix] = eigdec((Kxx+Kxx')/2, min(400, floor(T/4))); % /2
covfunc = {'covSum', {'covSEard','covNoise'}};
logtheta0 = [log(Wt); 0; log(sqrt(0.1))];
fprintf('Optimization hyperparameters in GP regression:\n');
IIx = find(eig_Kx > max(eig_Kx) * Thresh); eig_Kx = eig_Kx(IIx); eix = eix(:,IIx);
[logtheta_x, fvals_x, iter_x] = minimize(logtheta0, 'gpr_multi', -350, covfunc, [(1:T)'], 2*sqrt(T) *eix * diag(sqrt(eig_Kx))/sqrt(eig_Kx(1)));
% exp(logtheta_x),
if(logtheta_x(1)>log(1e4))
    logtheta_x(1)=log(1e4);
end
    
% Ktt = kernel((1:T)', (1:T)', [1/exp(2*logtheta_x(1)),1]);
Ktt = feval(covfunc_z{:}, logtheta_x(1:2), (1:T)');
invK2 = pdinv(Ktt + exp(2*logtheta_x(end))*eye(T));
Ml2 = Ktt*invK2*Kxx*invK2*Ktt;
% the square distance
D2 = diag(diag(Ml2)) * ones(size(Ml2)) + ones(size(Ml2)) * diag(diag(Ml2)) - 2*Ml2;
% Gaussian kernel
sigma2_square2 = median( D2(find(tril(ones(size(D2)),-1))) );
Mg2 = exp(-D2/sigma2_square2/2);

%% 
H = eye(T)-1/T*ones(T,T);
Mg = H*Mg*H;
Mg2 = H*Mg2*H;
testStat = 1/T^2*sum(sum(Mg'.*Mg2));
% eta = 1e-6;
% Rg = Mg*pdinv(Mg+T*eta*eye(T));
% Rg2 = Mg2*pdinv(Mg2+T*eta*eye(T));
% % testStat = sum(sum(Rg'.*Rg2));
% testStat = trace(Rg*Rg2);




