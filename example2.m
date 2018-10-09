% example 2
clear all,clc,close all
addpath(genpath(pwd))

% x1->x2->x3, and the causal module of x1, x2, and x3 are nonstationary,
% and the causal modules change independently
load smooth_module
% R0 saves generated nonstatioanry driving force which are independent of each other
T = 500;
x1 = 0.5*randn(T,1) + 5*R0{1}(1:T); 
x2 = 0.8*x1 + 4*R0{2}(1:T) + 0.5*randn(T,1);
x3 = 6*R0{6}(1:T)+ 0.8*x2 + 0.3*randn(T,1); 
Data = [x1,x2,x3];

% setting the parameters
alpha = 0.05; % signifcance level of independence test
maxFanIn=2; % maximum number of conditional variables
cond_ind_test='indtest_new_t';
[gns, g, SP] = nonsta_cd_new(Data,cond_ind_test,maxFanIn,alpha)

% INPUT: 
%       Data: - T*N matrix. T is number of data points and N is the number
%               of observed variables
%       cond_ind_test: - function handle that computes p-values for X ind. Y given Z: 
%                 (p_val = cond_ind_test(X, Y, Z, pars))
%       maxFanIn:  - maximum number of variables in the conditioning set 
%       alpha: - significance level of the independence test

% OUTPUT:
%       gns: - (n+1)*(n+1) matrix to represent recovered graph structure by
%       the methods for Markov equivalence class learning on augmented
%       causal graph & causal direction determination by making use of
%       independently changing causal modules
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
%    0  1  0  0
%    0  0  1  0
%    0  0  0  0
%    1  1  1  0

% estimated g: 
%    0  -1  0  0
%    0  0  -1  0
%    0  0  0  0
%    1  1  1  0




