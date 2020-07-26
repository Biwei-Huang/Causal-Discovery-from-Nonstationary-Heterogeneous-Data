function [g_skeleton, g_inv, gns, SP] = nonsta_cd_new_multi(X,dlabel,cond_ind_test,c_indx,maxFanIn,alpha, Type, pars)
% Constraint-based causal Discovery from Nonstationary/heterogeneous Data
% (CD-NOD)
% INPUT:
%       X: - T*n matrix. T is number of data points and n is the number
%               of observed variables
%       dlabel: - In the case with multi-dimensional variables, we use dlable to indicat the index of each variable
%       cond_ind_test: - function handle that computes p-values for X ind. Y given Z:
%                 (p_val = cond_ind_test(X, Y, Z, pars))
%       c_indx: surrogate variable to capture the distribution shift. If
%               data is nonstationary, then it is the time index. If data
%               is from multiple domains, then it is the domain index
%       maxFanIn:  - maximum number of variables in the conditioning set
%       alpha: - significance level of the independence test
%       Type: - run corresponding phases of CD-NOD
%          If Type=0, run all phases of CD-NOD (including
%             phase 1: learning causal skeleton,
%             phase 2: identifying causal directions with generalization of invariance,
%             phase 3: identifying directions with independent change principle, and
%             phase 4: recovering the nonstationarity driving force) )
%          If Type = 1, perform phase 1 + phase 2 + phase 3
%          If Type = 2, perform phase 1 + phase 2
%          If Type = 3, only perform phase 1
%       pars: - including pars.pairwise, pars.bonferroni, pars.if_GP1,
%               pars.if_GP2, pars.width, and pars.widthT
%          If pars.if_GP1 = 1, optimize the kernel width with GP in conditional independence tests;
%             otherwise, use a fixed kernel width
%          If pars.if_GP2 = 1, optimize the kernel width with GP in direction determination with
%             independent change principle & nonstationary driving force visualization
%          pars.width: kernel width on observational variables (except the time index).
%             If it is 0, then use the default kernel width when IF_GP1 = 0
%          pars.widthT: kernel width on the time index


% OUTPUT:
%       g_skeleton: - (n+1)*(n+1) matrix to represent recovered causal
%       skeleton over augmented set of variables
%            i-j: gns(i,j)=-1 & gns(j,i)=-1; i j: gns(i,j)=0 & gns(j,i)=0
%          - the last row of gns indicates the connection of nonstationarity
%            indicator (C) with other observed variables
%       g_inv: - (n+1)*(n+1) matrix to represent recovered graph structure
%       up to the Markov equivalence class learning on augmented
%       causal graph, with directions inferred by generalization of
%       invariance
%            i->j: g_inv(i,j)=1; i-j: g_inv(i,j)=-1; i j: g(i,j)=0
%          - the last row of g indicates the connection of nonstationarity
%            indicator (C) with other observed variables
%       gns: - (n+1)*(n+1) matrix to represent recovered graph structure,
%       with directions inferred by generalization of invariance &
%       independent change principle
%            i->j: gns(i,j)=1; i-j: gns(i,j)=-1; i j: gns(i,j)=0
%          - the last row of gns indicates the connection of nonstationarity
%            indicator (C) with other observed variables
%       ("gns" should have more oriented edges than "g_inv")
%       SP: - details of each independence test

% Copyright (c) 2017  Biwei Huang & Kun Zhang
% All rights reserved.

if ~isfield(pars,'pairwise')
    pars.pairwise = false;
end
if ~isfield(pars,'bonferroni')
    pars.bonferroni = false;
end
if ~isfield(pars,'width')
    pars.width = 0;
end
if ~isfield(pars,'widthT')
    pars.widthT = 0;
end
if ~isfield(pars,'if_GP1')
    pars.if_GP1 = 1;
end
if ~isfield(pars,'if_GP2')
    pars.if_GP2 = 1;
end

X=[X,c_indx]; % concatenate the surrogate variable C with others
X=X-repmat(mean(X),size(X,1),1);
X=X*diag(1./std(X));

n = length(dlabel)+1;% number of variables
dlabel{n} = size(X,2);

g_skeleton = []; g_inv = []; gns = []; SP = [];

%% Phase 1: learn the causal skeleton
% construct complete (fully connected) graph
g = eye(n) - ones(n,n);

% witness set
witness = zeros(n,n,n);

count =0;
% find graph skeleton by CPC
for s=0:maxFanIn % iteratively increase size of conditioning set
    for i=1:n
        % nodes adjacent to i
        adjSet = find(g(i,:)~=0);
        if (length(adjSet)<=s)
            continue;
        end
        % test whether i ind j | s
        for j=adjSet
            % unconditional test
            if (s==0)
                p_val = feval(cond_ind_test, X(:,dlabel{i}), X(:,dlabel{j}), [], pars);
                count=count+1;
                SP{count}=strcat('i=', int2str(i),', j=',int2str(j),', pval=',num2str(p_val));
                if  p_val > alpha
                    fprintf('%d ind %d with p-value %d\n', i,j,p_val);
                    g(i,j)=0;
                    g(j,i)=0;
                else
                    fprintf('%d notind %d with p-value %d\n', i,j,p_val);
                end
                continue;
            end
            
            % conditional independence test
            combs = nchoosek(adjSet(adjSet~=j),s);
            for k=1:size(combs,1)
                condSet = combs(k,:);
                % if independent
                condSet_string='{';
                for i2=1:(length(condSet)-1)
                    condSet_string=[condSet_string,int2str(condSet(i2)),', '];
                end
                condSet_label = [];
                for i2 = 1:length(condSet)
                    condSet_label = [condSet_label, dlabel{condSet}];
                end
                condSet_string=[condSet_string,int2str(condSet(length(condSet))),'}'];
                p_val = feval(cond_ind_test, X(:,dlabel{i}), X(:,dlabel{j}), X(:,condSet_label), pars);
                count=count+1;
                SP{count}=strcat('i=',int2str(i),', j=',int2str(j),', condSet_string=',condSet_string,', pval=',num2str(p_val));
                if p_val > alpha
                    fprintf('%d ind %d | %s with p-value %d\n', i,j,condSet_string, p_val);
                    witness(i,j,condSet) = ones(1,s);
                    witness(j,i,condSet) = ones(1,s);
                    g(i,j)=0;
                    g(j,i)=0;
                else
                    fprintf('%d notind %d | %s with p-value %d\n', i,j,condSet_string, p_val);
                end
            end
        end
    end
end
g_skeleton = g; % the causal skeleton over observed variables and C


%% phase 2: infer causal directions by generalizaion of invariance
% infer causal direction: C - X => C -> X
if(Type<=2)
    g(n,find(g(n,:)~=0)) = 1;
    g(find(g(:,n)~=0),n) = 0;
    % infer V-structures: X - Y - Z  => X -> Y <- Z
    adjMatrix = g;
    for i=1:n-1
        adj = find(adjMatrix(:,i)~=0);
        c = length(adj);
        for j=1:(c-1)
            for k=(j+1):c
                % check if moral
                if (adjMatrix(adj(j), adj(k))~=0 | adjMatrix(adj(k), adj(j))~=0)
                    continue;
                end
                % check to see if in witness set
                if (witness(adj(j), adj(k), i)==1 | witness(adj(k), adj(j), i)==1)
                    continue;
                end
                % orient immorality
                g(adj(j),i) = 1;
                g(i,adj(j)) = 0;
                g(adj(k),i) = 1;
                g(i,adj(k)) = 0;
            end
        end
    end
    
    % meeks rules to learn the Markov equivalence class
    g = meeks(-g,-adjMatrix);
    g = -g;
    g_inv = g; % skeleton + identified direction with generalization of invariance
end


%% Phase 3: infer causal directions by independent change principle
% infer the causal directions between two connected variables whose causal
% modules are both nonstationary
% see paper : "Behind Distribution Shift: Mining Driving Forces of Changes
% and Causal Arrows"
% You may comment out the following code, if you are only interested in the
% Markov equivalence class
if(Type<=1)
    Vns = find(g(n,:)==1); % find nodes with nonstationary causal modules
    Vns_un = []; % nodes with nonstationary causal modules and undirected edges
    for i = 1:length(Vns)
        if(~isempty(find(g(Vns(i),:)==-1)))
            Vns_un = [Vns_un,Vns(i)];
        end
    end
    Vns = Vns_un;
    gns = g;
    
    while(length(Vns)>1)
        score = [];
        hypo_eff = [];
        hypo_cau = [];
        for i = 1:length(Vns)
            hypo_eff{i} = Vns(i);
            hypo_cau{i} = union(find(gns(Vns(i),1:end-1)==-1),find(gns(1:end-1,Vns(i))==1)');
            hypo_eff_label = [];
            for j = 1:length(hypo_eff{i})
                hypo_eff_label = [hypo_eff_label,dlabel{hypo_eff{i}(j)}];
            end
            score(i) = infer_nonsta_dir(X(:,dlabel{hypo_cau{i}}),X(:,hypo_eff_label),c_indx,pars.width,pars.if_GP2);
        end
        [~, id] = min(score);
        sink = Vns(id); % looking for the sink of the current graph
        gns(hypo_cau{id},hypo_eff{id}) = 1;
        gns(hypo_eff{id},hypo_cau{id}) = 0;
        Vns = setdiff(Vns,sink);
    end
    % apply Meek's rule again
    gns = meeks(-gns,-adjMatrix);
    gns = -gns; % skeleton + identified directions with generalization of invariance & independent change principle
end


%% Phase 4: nonstationary driving force estimation
% visualize a low dimensional change of the causal module
% note that here we only concern the variability, but the scale can be
% different from the ground truth
if(Type==0)
    cmodule_id = find(g(n,:)==1); % variable id with changing modules
    for i = 1:length(cmodule_id)
        x_id = cmodule_id(i);
        pa_id = find(gns(1:n-1,x_id)==1)';
        
        x_label = dlabel{x_id};
        pa_label = [];
        for j = 1:length(pa_id)
           pa_label = [pa_label,dlabel{pa_id(j)}];
        end
        
        [Yg,Yl,Mg,Ml,D,eigValueg,eigValuel] = cd_non_con_fun(X(:,pa_label),X(:,x_label),c_indx,pars.width,pars.if_GP2);
        Yg_save{i} = Yg;
        Yl_save{i} = Yl;
        
        figure, subplot(211),plot(Yg(:,1),'b'); hold on; plot(Yg(:,2),'k--');
        title(strcat('Visualization of change in ',num2str(pa_id), '\rightarrow', num2str(x_id), ' (with Gaussian kernel)'))
        legend(['First component of \lambda_i (eigenvalue: ' num2str(eigValueg(1)) ')'],['Second component of \lambda_i (eigenvalue: ' num2str(eigValueg(2)) ')']);
        subplot(212),plot(abs(eigValueg),'.-');
        title('Absolute value of eigenvalues')
    end
end
