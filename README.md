# Causal-Discovery-from-Nonstationary-Heterogeneous-Data

Causal Discovery from Nonstationary/Heterogeneous Data. Copyright (c) 2017-2019 Biwei Huang & Kun Zhang

### MAIN FUNCTIONS
function [g_skeleton, g_inv, gns, SP] = nonsta_cd_new(X,cond_ind_test,c_indx,maxFanIn,alpha, Type, pars)

INPUT:
 *  X: - T*n matrix. T is number of data points and n is the number of observed variables
 *  cond_ind_test: - function handle that computes p-values for X ind. Y given Z: (p_val = cond_ind_test(X, Y, Z, pars))
 *  c_indx: surrogate variable to capture the distribution shift. If data is nonstationary, then it is the time index. If data is from multiple domains, then it is the domain index
 *  maxFanIn:  - maximum number of variables in the conditioning set
 *  alpha: - significance level of the independence test
 *  Type: - run corresponding phases of CD-NOD
   *  If Type=0, run all phases of CD-NOD (including phase 1: learning causal skeleton, phase 2: identifying causal directions with generalization of invariance, phase 3: identifying directions with independent change principle, and phase 4: recovering the nonstationarity driving force).
   *  If Type = 1, perform phase 1 + phase 2 + phase 3 
   *  If Type = 2, perform phase 1 + phase 2
   *  If Type = 3, only perform phase 1
 *  pars: - including pars.pairwise, pars.bonferroni, pars.if_GP1, pars.if_GP2, pars.width, and pars.widthT
   *  If pars.if_GP1 = 1, optimize the kernel width with GP in conditional independence tests; otherwise, use a fixed kernel width
   *  If pars.if_GP2 = 1, optimize the kernel width with GP in direction determination with independent change principle & nonstationary driving force visualization
   *  pars.width: kernel width on observational variables (except the time index). If it is 0, then use the default kernel width when IF_GP1 = 0
   *  pars.widthT: kernel width on the time index


OUTPUT:
 *  g_skeleton: - (n+1)*(n+1) matrix to represent recovered causal skeleton over augmented set of variables
   *  i-j: gns(i,j)=-1 & gns(j,i)=-1; i j: gns(i,j)=0 & gns(j,i)=0
   *  - the last row of gns indicates the connection of nonstationarity indicator (C) with other observed variables
 *  g_inv: - (n+1)*(n+1) matrix to represent recovered graph structure up to the Markov equivalence class learning on augmented causal graph, with directions inferred by generalization of invariance
   *  i->j: g_inv(i,j)=1; i-j: g_inv(i,j)=-1; i j: g(i,j)=0
   *  - the last row of g indicates the connection of nonstationarity indicator (C) with other observed variables
 *  gns: - (n+1)*(n+1) matrix to represent recovered graph structure, with directions inferred by generalization of invariance & independent change principle
   *  i->j: gns(i,j)=1; i-j: gns(i,j)=-1; i j: gns(i,j)=0
   *  - the last row of gns indicates the connection of nonstationarity indicator (C) with other observed variables
   *  ("gns" should have more oriented edges than "g_inv")
 *  SP: - details of each independence test


 ### EXAMPLE 
example1.m, example2.m, and example3.m give three example of using this package. 
Specifically, example1.m and example2.m are for nonstationary data, and example3.m is
for data from multiple domains.

### If there are multi-dimensional variables, use 
  function [g_skeleton, g_inv, gns, SP] = nonsta_cd_new_multi(X,dlabel,cond_ind_test,c_indx,maxFanIn,alpha, Type, pars)
  * dlabel: - In the case with multi-dimensional variables, we use dlable to indicat the index of each variable 
  Please see the example given in example4.m 

### Notes
For large-scale systems, there are several ways to speed up the process:

- Fix the kernel width in conditional independence tests: set "pars.if_GP1 = 0".
- Approximate the kernel learning with random fourier feature, by setting "cond_ind_test='indtest_new_t_RFF'".
- Do a pre-processing step to remove some spurious edges, e.g., first using partial correlation to remove some edges.

Note that if fixing the kernel width, you may need to tune the kernel width a bit to get better results,
especially for the kernel width on the time index ("pars.widthT"), "width" and "Wt" in "infer_nonsta_dir.m" and "cd_non_con_fun.m".

### CITATION
 If you use this code, please cite the following paper:

1.  Zhang, K., Huang, B., Zhang, J., Glymour, C., Scholkopf, B.. Causal Discovery from Nonstationary/Heterogeneous Data: Skeleton Estimation and Orientation Determination. IJCAI 2017.
2.  Huang, B., Zhang, K., Zhang, J., Glymour, C., Scholkopf, B. Behind Distribution Shift: Mining Driving Forces of Changes and Causal Arrows. ICDM 2017.
3.  Huang, B., Zhang, K., Zhang, J., Ramsey, J., Sanchez-Romero, R., Glymour, C., Scholkopf, B.. Causal Discovery from Heterogeneous/Nonstationary Data. JMLR, 21(89), 2020.

If you have problems or questions, do not hesitate to send an email to biweih@andrew.cmu.edu
