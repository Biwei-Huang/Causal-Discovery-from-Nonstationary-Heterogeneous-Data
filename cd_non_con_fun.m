function [Yg,Yl,Mg,Ml,D,eigValueg,eigValuel] = cd_non_con_fun(X,Y,c_indx,width,IF_GP)
% learn the nonstationary driving force of the causal mechanism
% X: parents; Y; effect
% width: the kernel width for X and Y
% c_indx: surrogate variable to capture the distribution shift; 
% If If_GP = 1, learning the kernel width for P(Y|X). Set it to 0 can speed up the process!!!

if(width==0)
    width = 0.1;
end
Wt = 1;  % the initial kernel width on C (or T). May need tunning for different data sets!!!
[T,d] = size(X);
X = X * diag(1./std(X));
Y = Y/std(Y);
theta = 1/width^2; % 0.2
lambda = 1; % 0.05 0.3  10
Ml = [];

% size of Y should be T*1.
Kyy = kernel(Y, Y, [theta,1]);

if IF_GP
    Thresh = 1E-4;
    [eig_Ky, eiy] = eigdec((Kyy+Kyy')/2, min(400, floor(T/4))); % /2
    covfunc = {'covSum', {'covSEard','covNoise'}};
    %     covfunc = {'covSum', {'covMatern3ard','covNoise'}};
    
    logtheta0 = [log(width)*ones(d,1);  log(Wt); 0; log(sqrt(0.1))];
    fprintf('Optimization hyperparameters in GP regression:\n');
    
    IIy = find(eig_Ky > max(eig_Ky) * Thresh); eig_Ky = eig_Ky(IIy); eiy = eiy(:,IIy);
    [logtheta_y, fvals_y, iter_y] = minimize(logtheta0, 'gpr_multi', -350, covfunc, [X c_indx], 2*sqrt(T) *eiy * diag(sqrt(eig_Ky))/sqrt(eig_Ky(1)));
    exp(logtheta_y),
    
    covfunc_z = {'covSEard'};
    Kxt = feval(covfunc_z{:}, logtheta_y, [X c_indx]);
    
    % Note: in the conditional case, no need to do centering, as the regression
    % will automatically enforce that.
    
    % Kernel matrices of the errors
    invK = pdinv(Kxt + exp(2*logtheta_y(end))*eye(T));
    
%         Kxx = kernel(X, X, [1/exp(2*logtheta_y(1)),1]);
%         Ktt = kernel((1:T)', (1:T)', [1/exp(2*logtheta_y(d+1)),1]);
    Kxx = feval(covfunc_z{:}, logtheta_y([1:d,d+2]), X);
    Ktt = feval(covfunc_z{:}, logtheta_y([d+1,d+2]), c_indx);
else
    Kxx = kernel(X, X, [theta,1]);
    Kyy = kernel(Y, Y, [theta,1]);
    Ktt = kernel(c_indx, c_indx, [1/Wt^2,1]);
    invK = pdinv( Kxx.* Ktt +  lambda * eye(T));
end
Kxx3 = Kxx^3; %^3

prod_invK =  invK * Kyy * invK;
% now finding Ml

Ml = 1/T^2 * Ktt*( Kxx3 .* prod_invK) * Ktt;

% Len = floor(T/50);
% for c = 1:50
%     cc = Len*(c-1)+1;
%     for c1 = c:50
%         fprintf('.');
%         cc1 = Len*(c1-1)+1;
%         % Ml(c,c1) = trace(diag(Ktt(:,cc)) *  Kxx3 * diag(Ktt(:,cc1)) * prod_invK);
%         Ml(c,c1) = trace( ((Ktt(:,cc) * Ktt(:,cc1)' ) .* Kxx3) * prod_invK );
%         if c1>c
%             Ml(c1,c) = Ml(c,c1);
%         end
%     end
% end
% Ml = 1/T^2 * Ml;

% the square distance
D = diag(diag(Ml)) * ones(size(Ml)) + ones(size(Ml)) * diag(diag(Ml)) - 2*Ml;

% Gaussian kernel
sigma2_square = median( D(find(tril(ones(size(D)),-1))) );
Mg = exp(-D/sigma2_square/2);


[Yg, eigVectorg, eigValueg]=kPCA_kernel_orig(Mg,3);
[Yl, eigVectorl, eigValuel]=kPCA_kernel_orig(Ml,3);
% figure, plot(Yg(:,1),'b'); hold on; plot(Yg(:,2),'k--'); title('VIsualization of change in PA^i \rightarrow V^i (with Gaussian kernel)')
% legend(['First component of \lambda_i (eigenvalue: ' num2str(eigValueg(1)) ')'],['Second component of \lambda_i (eigenvalue: ' num2str(eigValueg(2)) ')']);

% figure, plot(Yl(:,1),'b'); hold on; plot(Yl(:,2),'k--'); title('VIsualization of change in PA^i \rightarrow V^i (with linear kernel)')
% legend(['First component of \lambda_i (eigenvalue: ' num2str(eigValuel(1)) ')'],['Second component of \lambda_i (eigenvalue: ' num2str(eigValuel(2)) ')']);
