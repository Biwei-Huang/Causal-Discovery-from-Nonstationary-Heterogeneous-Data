# Causal-Discovery-from-Nonstationary-Heterogeneous-Data

Causal Discovery from Nonstationary/Heterogeneous Data. Copyright (c) 2017-2018 Kun Zhang & Biwei Huang

### MAIN FUNCTIONS
function [gns, g, SP] = nonsta_cd_new(X,cond_ind_test,maxFanIn,alpha) 

INPUT: 
 *  Data: - T*n matrix. T is number of data points and n is the number of observed variables 
 *  cond_ind_test: - function handle that computes p-values for X ind. Y given Z: (p_val = cond_ind_test(X, Y, Z, pars))
 *  maxFanIn: - maximum number of variables in the conditioning set 
 *  alpha: - significance level of the independence test

OUTPUT: 
 * gns: - (n+1)_(n+1) matrix to represent recovered graph structure by the methods for Markov equivalence class learning on augmented causal graph & causal direction determination by making use of independently changing causal modules 
   * i->j: gns(i,j)=1; i-j: gns(i,j)=-1; i j: gns(i,j)=0 
   * the last row of gns indicates the connection of nonstationarity indicator (C) with other observed variables 
 * g: - (n+1)_(n+1) matrix to represent recovered graph structure only by % the methods for Markov equivalence class learning on augmented causal graph 
   * i->j: g(i,j)=1; i-j: g(i,j)=-1; i j: g(i,j)=0 
   * the last row of g indicates the connection of nonstationarity indicator (C) with other observed variables 
   * ("gns" should have more oriented edges than "g") 
 * SP: - details of each independence test
 
 ### EXAMPLE 
example1.m and example2.m give two examples in using this package.

### Notes
For large-scale systems, there are several ways to speed up the process:
- Fix the kernel width with in conditional independence tests: set IF_GP = 1 at ./KCI-test/algorithms/CInd_test_new_withGP_t.
- Approximate the kernel learning with random fourier feature, Nystrom Method, or other kernel approximation techniques.
- Do a pre-processing step to remove some spurious edges, e.g., first using partial correlation to remove some edges.

### CITATION
 If you use this code, please cite the following paper:

1.  "Zhang, K., Huang, B., Zhang, J., Glymour, C., Schölkopf, B.. Causal Discovery from Nonstationary/Heterogeneous Data: Skeleton Estimation and Orientation Determination. IJCAI 2017."
2.  "Huang, B., Zhang, K., Zhang, J., Glymour, C., Schölkopf, B. Behind Distribution Shift: Mining Driving Forces of Changes and Causal Arrows. ICDM 2017."

If you have problems or questions, do not hesitate to send an email to  [biweih@andrew.cmu.edu](mailto:biweih@andrew.cmu.edu)
