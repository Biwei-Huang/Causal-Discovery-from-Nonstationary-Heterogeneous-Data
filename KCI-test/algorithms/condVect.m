function[res] = condVect(S,T,U,M,ths,BG)
%calculates independence measure on vectors
%
%S,T,U: disjoint lists of of indizes
%M: covariance matrix of vectors
%ths: threshold for independence if ths=0, then condVect returns I(S:T|U),
%if ths > 0, then condVect returns 1, if I(S:T|U)<=ths and 0 otherwise
%BG: list of indizes of background variables (are always conditioned on but do not count as variables for PC)
%
%example: 
%condVect([1,2],[3],[4,5],M,0,[6]) returns I(1,2 : 3 | 4,5,6)
%condVect([1,2],[3],[],M,0.5,[]) returns 1, if  I(1,2 : 3)<=0.5, and 0
%otherwise
% Copyright (c) 2010-2011  ...
% All rights reserved.  See the file COPYING for license terms.

if (ths>=0)
fprintf('Calculating I(%s : %s | %s)  ',num2str(S,'%1.0d '),num2str(T,' %1.0d '),num2str(U,' %1.0d'));
end
U=union(U,BG);
I = union(S,T);I=union(I,U); %set of indices

res = entVect(M(union(S,U),union(S,U)))+entVect(M(union(T,U),union(T,U)))-entVect(M(U,U))-entVect(M(I,I));

if (ths >= 0)
    fprintf(' Result: %1.2d\n',res);
    res = (res<=ths);
end
end

%information from covariance matrix
function[res] = entVect(M)

res = 1/2*log(det(M));
end
