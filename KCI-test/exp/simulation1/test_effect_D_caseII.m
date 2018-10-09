% script to fnd the prob. of type I and II errors in case II of Simulation 1.

% script to test how the p values change with the
% dimensionality of the conditional set

T = 200; % T = 200 or 400
It_max = 1000; % you may change the total number of replications... 5000...
Width = 0; % automatically determines the kernel width depending on the sample size
TT = zeros(5,It_max); % Time used by our method
TT1 = zeros(5,It_max); % Time used by CI_PERM

% Independent case: type I error
for iter = 1:It_max
    if ~mod(iter,50)
        fprintf('%d replications in finding the type I errors in Case II...\n', iter),
    end
    x = normrnd(0,1,T,1);
    y = normrnd(0,1,T,1);
    z = normrnd(0,1,T,1);
    zz1 = .7 * ( (z.^3)/5 + z/2);
    x = zz1 + tanh(x);
    x = x + (x.^3)/3 + tanh(x/3)/2;
    zz2 = (z.^3/4 + z)/3;
    y = y + zz2;
    y = y + tanh(y/3);
    x = x - mean(x); x = x/std(x);
    y = y - mean(y); y = y/std(y);
    z = z - mean(z); z = z/std(z);
    
    % testing...
    % conditioning on one variable
    t_t = cputime;
    % our method
    [tmp1, tmp2, tmp3, tmp4, tmp5] =...
        CInd_test_new_withGP(x, y, [z], 0.01, Width);
    TT(1,iter) = cputime - t_t;
    Sta_1(iter) = tmp1; Cri_1(iter) = tmp2; p_val_1(iter) = tmp3; Cri_appr_1(iter) = tmp4;
    p_appr_1(iter) = tmp5;
    t_t = cputime;
    % CI_PERM
    [pval, stat] = indtest_hsic(x, y, z, 'perm');
    TT1(1,iter) = cputime - t_t;
    p_hsic_1(iter) = pval;
    % with two conditioning variables
    %% this is for the non-spurious Z case
    z2 = normrnd(0,1,T,1);
    x = normrnd(0,1,T,1);
    y = normrnd(0,1,T,1);
    zz1_2 = zz1/2 + z2;
    zz1_2 = zz1_2/2 + .7 * tanh(zz1_2);
    x = zz1_2 + tanh(x);
    x = x + (x.^3)/3 + tanh(x/3)/2;
    zz2_2 = zz2/2 + z2;
    zz2_2 = zz2_2/2 + .7 * tanh(zz2_2);
    y = y + zz2_2;
    y = y + tanh(y/3);
    x = x - mean(x); x = x/std(x);
    y = y - mean(y); y = y/std(y);
    z2 = z2 - mean(z2); z2 = z2/std(z2); z = [z z2];
    
    t_t = cputime;
    [tmp1, tmp2, tmp3, tmp4, tmp5] =...
        CInd_test_new_withGP(x, y, [z], 0.01, Width);
    TT(2,iter) = cputime - t_t;
    Sta_2(iter) = tmp1; Cri_2(iter) = tmp2; p_val_2(iter) = tmp3; Cri_appr_2(iter) = tmp4;
    p_appr_2(iter) = tmp5;
    t_t = cputime;
    [pval, stat] = indtest_hsic(x, y, z, 'perm');
    TT1(2,iter) = cputime - t_t;
    p_hsic_2(iter) = pval;
    % with three conditioning variables
    %% this is for the non-spurious Z case
    z3 = normrnd(0,1,T,1);
    x = normrnd(0,1,T,1);
    y = normrnd(0,1,T,1);
    zz1_3 = zz1_2*2/3 + z3*5/6;
    zz1_3 = zz1_3/2 + .7 * tanh(zz1_3);
    x = zz1_3 + tanh(x);
    x = x + (x.^3)/3 + tanh(x/3)/2;
    zz2_3 = zz2_2*2/3 + z3*5/6;
    zz2_3 = zz2_3/2 + .7 * tanh(zz2_3);
    y = y + zz2_3;
    y = y + tanh(y/3);
    x = x - mean(x); x = x/std(x);
    y = y - mean(y); y = y/std(y);
    z3 = z3 - mean(z3); z3 = z3/std(z3); z = [z z3];
    
    t_t = cputime;
    [tmp1, tmp2, tmp3, tmp4, tmp5] =...
        CInd_test_new_withGP(x, y, [z], 0.01, Width);
    TT(3,iter) = cputime - t_t;
    Sta_3(iter) = tmp1; Cri_3(iter) = tmp2; p_val_3(iter) = tmp3; Cri_appr_3(iter) = tmp4;
    p_appr_3(iter) = tmp5;
    t_t = cputime;
    [pval, stat] = indtest_hsic(x, y, z, 'perm');
    TT1(3,iter) = cputime - t_t;
    p_hsic_3(iter) = pval;
    % with four conditioning variables
    %% this is for the non-spurious Z case
    z4 = normrnd(0,1,T,1);
    x = normrnd(0,1,T,1);
    y = normrnd(0,1,T,1);
    zz1_4 = zz1_3*2/3 + z4*5/6;
    zz1_4 = zz1_4/2 + .7 * tanh(zz1_4);
    x = zz1_4 + tanh(x);
    x = x + (x.^3)/3 + tanh(x/3)/2;
    zz2_4 = zz2_3*2/3 + z4*5/6;
    zz2_4 = zz2_4/2 + .7 * tanh(zz2_4);
    y = y + zz2_4;
    y = y + tanh(y/3);
    x = x - mean(x); x = x/std(x);
    y = y - mean(y); y = y/std(y);
    z4 = z4 - mean(z4); z4 = z4/std(z4); z = [z z4];
    
    t_t = cputime;
    [tmp1, tmp2, tmp3, tmp4, tmp5] =...
        CInd_test_new_withGP(x, y, [z], 0.01, Width);
    TT(4,iter) = cputime - t_t;
    Sta_4(iter) = tmp1; Cri_4(iter) = tmp2; p_val_4(iter) = tmp3; Cri_appr_4(iter) = tmp4;
    p_appr_4(iter) = tmp5;
    t_t = cputime;
    [pval, stat] = indtest_hsic(x, y, z, 'perm');
    TT1(4,iter) = cputime - t_t;
    p_hsic_4(iter) = pval;
    
    % with five conditioning variables
    %% this is for the non-spurious Z case
    z5 = normrnd(0,1,T,1);
    x = normrnd(0,1,T,1);
    y = normrnd(0,1,T,1);
    zz1_5 = zz1_4*2/3 + z5*5/6;
    zz1_5 = zz1_5/2 + .7 * tanh(zz1_5);
    x = zz1_5 + tanh(x);
    x = x + (x.^3)/3 + tanh(x/3)/2;
    zz2_5 = zz2_4*2/3 + z5*5/6;
    zz2_5 = zz2_5/2 + .7 * tanh(zz2_5);
    y = y + zz2_5;
    y = y + tanh(y/3);
    x = x - mean(x); x = x/std(x);
    y = y - mean(y); y = y/std(y);
    z5 = z5 - mean(z5); z5 = z5/std(z5); z = [z z5];
    
    t_t = cputime;
    [tmp1, tmp2, tmp3, tmp4, tmp5] =...
        CInd_test_new_withGP(x, y, [z], 0.01, Width);
    TT(5,iter) = cputime - t_t;
    Sta_5(iter) = tmp1; Cri_5(iter) = tmp2; p_val_5(iter) = tmp3; Cri_appr_5(iter) = tmp4;
    p_appr_5(iter) = tmp5;
    t_t = cputime;
    [pval, stat] = indtest_hsic(x, y, z, 'perm');
    TT1(5,iter) = cputime - t_t;
    p_hsic_5(iter) = pval;
    save new_compare_Ind_CaseII.mat
    save TT_ourmethod.txt TT -ASCII
    save TT1_CI_PERM.txt TT1 -ASCII
end

Error_bs(1) = sum(p_val_1(1:It_max)<0.01)/It_max; Error_a(1) = sum(p_appr_1(1:It_max)<0.01)/It_max;
Error_bs(2) = sum(p_val_2(1:It_max)<0.01)/It_max; Error_a(2) = sum(p_appr_2(1:It_max)<0.01)/It_max;
Error_bs(3) = sum(p_val_3(1:It_max)<0.01)/It_max; Error_a(3) = sum(p_appr_3(1:It_max)<0.01)/It_max;
Error_bs(4) = sum(p_val_4(1:It_max)<0.01)/It_max; Error_a(4) = sum(p_appr_4(1:It_max)<0.01)/It_max;
Error_bs(5) = sum(p_val_5(1:It_max)<0.01)/It_max; Error_a(5) = sum(p_appr_5(1:It_max)<0.01)/It_max;

Error_bs5(1) = sum(p_val_1(1:It_max)<0.05)/It_max; Error_a5(1) = sum(p_appr_1(1:It_max)<0.05)/It_max;
Error_bs5(2) = sum(p_val_2(1:It_max)<0.05)/It_max; Error_a5(2) = sum(p_appr_2(1:It_max)<0.05)/It_max;
Error_bs5(3) = sum(p_val_3(1:It_max)<0.05)/It_max; Error_a5(3) = sum(p_appr_3(1:It_max)<0.05)/It_max;
Error_bs5(4) = sum(p_val_4(1:It_max)<0.05)/It_max; Error_a5(4) = sum(p_appr_4(1:It_max)<0.05)/It_max;
Error_bs5(5) = sum(p_val_5(1:It_max)<0.05)/It_max; Error_a5(5) = sum(p_appr_5(1:It_max)<0.05)/It_max;
figure, plot(Error_bs, 'ro--'), hold on, plot(Error_a, 'kx--'),
plot(Error_bs5, 'r^--'), hold on, plot(Error_a5, 'kh--'), 

legend('Bootstrap (\alpha=0.01)','Gamma appr. (\alpha=0.01)', 'Bootstrap (\alpha=0.05)','Gamma appr. (\alpha=0.05)');
title(['Type I error (T=' int2str(T) ')']);
save new_compare_Ind_CaseII.mat

Dependent case: type II error
for iter = 1:It_max 
    if ~mod(iter,50)
        fprintf('%d replications in finding the type II errors in Case II...\n', iter),
    end
    x = randn(T,1); % normrnd(0,1,T,1);
    y = randn(T,1); % normrnd(0,1,T,1);
    z = randn(T,1); % normrnd(0,1,T,1);
    zz1 = .7 * ( (z.^3)/5 + z/2);
    x = zz1 + tanh(x);
    x = x + (x.^3)/3 + tanh(x/3)/2;
    zz2 = (z.^3/4 + z)/3;
    y = y + zz2;
    y = y + tanh(y/3);
    x = x - mean(x); x = x/std(x);
    y = y - mean(y); y = y/std(y);
    z = z - mean(z); z = z/std(z);
    ff = randn(T,1)*0.5;
    x = x + ff; y = y + ff; % to make them conditionally DEPENDENT
    
    % testing...
    % conditioning on one variable
    t_t = cputime;
    [tmp1, tmp2, tmp3, tmp4, tmp5] =...
        CInd_test_new_withGP(x, y, [z], 0.01, Width);
    TT(1,iter) = cputime - t_t;
    Sta_1(iter) = tmp1; Cri_1(iter) = tmp2; p_val_1(iter) = tmp3; Cri_appr_1(iter) = tmp4;
    p_appr_1(iter) = tmp5;
    t_t = cputime;
    [pval, stat] = indtest_hsic(x, y, z, 'perm');
    TT1(1,iter) = cputime - t_t;
    p_hsic_1(iter) = pval;
    % on two variables
    %% this is for the non-spurious Z case
    z2 = normrnd(0,1,T,1);
    x =  normrnd(0,1,T,1);
    y = normrnd(0,1,T,1);
    zz1_2 = zz1/2 + z2;
    zz1_2 = zz1_2/2 + .7 * tanh(zz1_2);
    x = zz1_2 + tanh(x);
    x = x + (x.^3)/3 + tanh(x/3)/2;
    zz2_2 = zz2/2 + z2;
    zz2_2 = zz2_2/2 + .7 * tanh(zz2_2);
    y = y + zz2_2;
    y = y + tanh(y/3);
    x = x - mean(x); x = x/std(x);
    y = y - mean(y); y = y/std(y);
    z2 = z2 - mean(z2); z2 = z2/std(z2); z = [z z2];
    x = x + ff; y = y + ff;
    
    t_t = cputime;
    [tmp1, tmp2, tmp3, tmp4, tmp5] =...
        CInd_test_new_withGP(x, y, [z], 0.01, Width);
    TT(2,iter) = cputime - t_t;
    Sta_2(iter) = tmp1; Cri_2(iter) = tmp2; p_val_2(iter) = tmp3; Cri_appr_2(iter) = tmp4;
    p_appr_2(iter) = tmp5;
    t_t = cputime;
    [pval, stat] = indtest_hsic(x, y, z, 'perm');
    TT1(2,iter) = cputime - t_t;
    p_hsic_2(iter) = pval;
    % on three variables
    %% this is for the non-spurious Z case
    z3 = normrnd(0,1,T,1);
    x = normrnd(0,1,T,1);
    y = normrnd(0,1,T,1);
    zz1_3 = zz1_2*2/3 + z3*5/6;
    zz1_3 = zz1_3/2 + .7 * tanh(zz1_3);
    x = zz1_3 + tanh(x);
    x = x + (x.^3)/3 + tanh(x/3)/2;
    zz2_3 = zz2_2*2/3 + z3*5/6;
    zz2_3 = zz2_3/2 + .7 * tanh(zz2_3);
    y = y + zz2_3;
    y = y + tanh(y/3);
    x = x - mean(x); x = x/std(x);
    y = y - mean(y); y = y/std(y);
    z3 = z3 - mean(z3); z3 = z3/std(z3); z = [z z3];
    x = x + ff; y = y + ff;
    
    t_t = cputime;
    [tmp1, tmp2, tmp3, tmp4, tmp5] =...
        CInd_test_new_withGP(x, y, [z], 0.01, Width);
    TT(3,iter) = cputime - t_t;
    Sta_3(iter) = tmp1; Cri_3(iter) = tmp2; p_val_3(iter) = tmp3; Cri_appr_3(iter) = tmp4;
    p_appr_3(iter) = tmp5;
    t_t = cputime;
    [pval, stat] = indtest_hsic(x, y, z, 'perm');
    TT1(3,iter) = cputime - t_t;
    p_hsic_3(iter) = pval;
    % on four variables
    %% this is for the non-spurious Z case
    z4 = normrnd(0,1,T,1);
    x = normrnd(0,1,T,1);
    y = normrnd(0,1,T,1);
    zz1_4 = zz1_3*2/3 + z4*5/6;
    zz1_4 = zz1_4/2 + .7 * tanh(zz1_4);
    x = zz1_4 + tanh(x);
    x = x + (x.^3)/3 + tanh(x/3)/2;
    zz2_4 = zz2_3*2/3 + z4*5/6;
    zz2_4 = zz2_4/2 + .7 * tanh(zz2_4);
    y = y + zz2_4;
    y = y + tanh(y/3);
    x = x - mean(x); x = x/std(x);
    y = y - mean(y); y = y/std(y);
    z4 = z4 - mean(z4); z4 = z4/std(z4); z = [z z4];
    x = x + ff; y = y + ff;
    
    t_t = cputime;
    [tmp1, tmp2, tmp3, tmp4, tmp5] =...
        CInd_test_new_withGP(x, y, [z], 0.01, Width);
    TT(4,iter) = cputime - t_t;
    Sta_4(iter) = tmp1; Cri_4(iter) = tmp2; p_val_4(iter) = tmp3; Cri_appr_4(iter) = tmp4;
    p_appr_4(iter) = tmp5;
    t_t = cputime;
    [pval, stat] = indtest_hsic(x, y, z, 'perm');
    TT1(4,iter) = cputime - t_t;
    p_hsic_4(iter) = pval;
    % conditional on 5 variables...
    %% this is for the non-spurious Z case
    z5 = normrnd(0,1,T,1);
    x = normrnd(0,1,T,1);
    y = normrnd(0,1,T,1);
    zz1_5 = zz1_4*2/3 + z5*5/6;
    zz1_5 = zz1_5/2 + .7 * tanh(zz1_5);
    x = zz1_5 + tanh(x);
    x = x + (x.^3)/3 + tanh(x/3)/2;
    zz2_5 = zz2_4*2/3 + z5*5/6;
    zz2_5 = zz2_5/2 + .7 * tanh(zz2_5);
    y = y + zz2_5;
    y = y + tanh(y/3);
    x = x - mean(x); x = x/std(x);
    y = y - mean(y); y = y/std(y);
    z5 = z5 - mean(z5); z5 = z5/std(z5); z = [z z5];
    x = x + ff; y = y + ff;
    
    t_t = cputime;
    [tmp1, tmp2, tmp3, tmp4, tmp5] =...
        CInd_test_new_withGP(x, y, [z], 0.01, Width);
    TT(5,iter) = cputime - t_t;
    Sta_5(iter) = tmp1; Cri_5(iter) = tmp2; p_val_5(iter) = tmp3; Cri_appr_5(iter) = tmp4;
    p_appr_5(iter) = tmp5;
    t_t = cputime;
    [pval, stat] = indtest_hsic(x, y, z, 'perm');
    TT1(5,iter) = cputime - t_t;
    p_hsic_5(iter) = pval;
    save new_compare_Dep_CaseII.mat
end

Error_bs(1) = sum(p_val_1(1:It_max)>0.01)/It_max; Error_a(1) = sum(p_appr_1(1:It_max)>0.01)/It_max;
Error_bs(2) = sum(p_val_2(1:It_max)>0.01)/It_max; Error_a(2) = sum(p_appr_2(1:It_max)>0.01)/It_max;
Error_bs(3) = sum(p_val_3(1:It_max)>0.01)/It_max; Error_a(3) = sum(p_appr_3(1:It_max)>0.01)/It_max;
Error_bs(4) = sum(p_val_4(1:It_max)>0.01)/It_max; Error_a(4) = sum(p_appr_4(1:It_max)>0.01)/It_max;
Error_bs(5) = sum(p_val_5(1:It_max)>0.01)/It_max; Error_a(5) = sum(p_appr_5(1:It_max)>0.01)/It_max;

Error_bs5(1) = sum(p_val_1(1:It_max)>0.05)/It_max; Error_a5(1) = sum(p_appr_1(1:It_max)>0.05)/It_max;
Error_bs5(2) = sum(p_val_2(1:It_max)>0.05)/It_max; Error_a5(2) = sum(p_appr_2(1:It_max)>0.05)/It_max;
Error_bs5(3) = sum(p_val_3(1:It_max)>0.05)/It_max; Error_a5(3) = sum(p_appr_3(1:It_max)>0.05)/It_max;
Error_bs5(4) = sum(p_val_4(1:It_max)>0.05)/It_max; Error_a5(4) = sum(p_appr_4(1:It_max)>0.05)/It_max;
Error_bs5(5) = sum(p_val_5(1:It_max)>0.05)/It_max; Error_a5(5) = sum(p_appr_5(1:It_max)>0.05)/It_max;

figure, plot(Error_bs, 'ro--'), hold on, plot(Error_a, 'kx--'),
plot(Error_bs5, 'r^--'), hold on, plot(Error_a5, 'kh--'), 
legend('Bootstrap (\alpha=0.01)','Gamma appr. (\alpha=0.01)', 'Bootstrap (\alpha=0.05)','Gamma appr. (\alpha=0.05)');title(['Type II error (T=' int2str(T) ')']);
title(['Type II error (T=' int2str(T) ')']);
save new_compare_Dep_CaseII.mat
