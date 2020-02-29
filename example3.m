% example 1: heterogeneous data (data from multiple domains)
clear all,clc,close all
addpath(genpath(pwd))
rng(10)
% generate data from the first domain
T_1 = 300;
x1_1 = randn(T_1,1);
x2_1 = 0.9*x1_1 + 0.6*randn(T_1,1);
x3_1 = 0.9*x2_1 + 0.6*randn(T_1,1); 
x4_1 = 0.9*x3_1 + 0.6*randn(T_1,1); 
Data_1 = [x1_1,x2_1,x3_1,x4_1];

% generate data from the second domain
T_2 = 300;
x1_2 = randn(T_2,1);
x2_2 = sin(x1_2) + 0.2*randn(T_2,1);
x3_2 = sin(x2_2) + 0.2*randn(T_2,1); 
x4_2 = sin(x3_2) + 0.2*randn(T_2,1); 
Data_2 = [x1_2,x2_2,x3_2,x4_2];

% concateneate data from the two domains
Data = [Data_1;Data_2];

% setting the parameters
alpha = 0.05; % signifcance level of independence test
maxFanIn = 2; % maximum number of conditional variables
cond_ind_test='indtest_new_t';
c_indx = [ones(1,T_1),2*ones(1,T_2)]'; % surrogate variable to capture the distribution shift; 
                 %here it is the doamin index, because the data is from multiple domains
[gns, g, SP] = nonsta_cd_new(Data,cond_ind_test,c_indx,maxFanIn,alpha)

% INPUT: 
%       Data: - T*n matrix. T is number of data points and n is the number
%               of observed variables
%       cond_ind_test: - function handle that computes p-values for X ind. Y given Z: 
%                 (p_val = cond_ind_test(X, Y, Z, pars))
%       c_indx: surrogate variable to capture the distribution shift. If
%               the data is nonstationary, then it is the time index. If the data
%               is from multiple domains, then it is the domain index
%       maxFanIn:  - maximum number of variables in the conditioning set 
%       alpha: - significance level of the independence test

% OUTPUT:
%       gns: - (n+1)*(n+1) matrix to represent recovered graph structure by
%       the methods for Markov equivalence class learning on augmented
%       causal graph & causal direction determination by making use of
%       independent changing causal modules
%            i->j: gns(i,j)=1; i-j: gns(i,j)=-1; i j: gns(i,j)=0
%          - the last row of gns indicates the connection of nonstationarity
%            indicator (C) with other observed variables
%       g: - (n+1)*(n+1) matrix to represent recovered graph structure only by
%       the methods for Markov equivalence class learning on augmented
%       causal graph
%            i->j: g(i,j)=1; i-j: g(i,j)=-1; i j: g(i,j)=0
%          - the last row of g indicates the connection of nonstationarity
%            indicator (C) with other observed variables
%       ("gns" should have more oriented edges than "g")
%       SP: - details of each independence test


% estimated gns: 
%  0  1  0  0  0 
%  0  0  1  0  0 
%  0  0  0  1  0 
%  0  0  0  0  0 
%  0  1  0  1  0 


% estimated g: 
%  0  1  0  0  0 
%  0  0  1  0  0 
%  0  0  0  1  0 
%  0  0  0  0  0 
%  0  1  0  1  0 
 
