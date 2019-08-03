% example 1
clear all,clc,close all
addpath(genpath(pwd))

% x1->x2->x3->x4, and the causal modules of x2 and x4 is nonstationary, and
% their changes are related 
T = 500;
x1 = randn(T,1);
x2 = 0.8*x1 + 1.5*sin([1:T]'/50) + 0.5*randn(T,1);
x3 = 0.8*x2 + 0.5*randn(T,1); 
x4 = 0.8*x3 + (sin([1:T]'/50)+sin([1:T]'/20)) + 0.5*randn(T,1); 
Data = [x1,x2,x3,x4];

% setting the parameters
alpha = 0.05; % signifcance level of independence test
maxFanIn=2; % maximum number of conditional variables
cond_ind_test='indtest_new_t';
[gns, g, SP] = nonsta_cd_new(Data,cond_ind_test,maxFanIn,alpha)

% INPUT: 
%       Data: - T*n matrix. T is number of data points and n is the number
%               of observed variables
%       cond_ind_test: - function handle that computes p-values for X ind. Y given Z: 
%                 (p_val = cond_ind_test(X, Y, Z, pars))
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
 