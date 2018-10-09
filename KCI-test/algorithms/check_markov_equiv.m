function [res]=check_markov_equiv(g1,g2)
% function res=check_markov_equiv(g1,g2)
% INPUT: two graphs g1, g2. 
%		g(i,j)=-1       if there is a directed arrow from i to j.
%		g(i,j)=g(j,i)=1 if there is an undirected edge between i and j 
%
% OUTPUT: res==1: the two graphs are markov equivalent
% 	  res==0: they are not
%
% Copyright (c) 2010-2011  ...
% All rights reserved.  See the file COPYING for license terms. 


res=1;
num_nodes=size(g1,1);


%check whether they have the same skeleton
%
skeleton1=g1+g1'; skeleton1(skeleton1~=0)=skeleton1(skeleton1~=0)./skeleton1(skeleton1~=0);
skeleton2=g2+g2'; skeleton2(skeleton2~=0)=skeleton2(skeleton2~=0)./skeleton2(skeleton2~=0);
if ~isequal(skeleton1,skeleton2)
    res=0;
    fprintf('not the same skeletons\n');
end

if res==1
%check whether they have the same set of immoralites
%
    for i=1:num_nodes
        i_parents=find(g1(:,i)==-1);
        for ii1=1:(length(i_parents)-1)
            for ii2=(ii1+1):length(i_parents)
                if g1(i_parents(ii2),i_parents(ii1))==0 & g1(i_parents(ii1),i_parents(ii2))==0
                    if g2(i,i_parents(ii1))~=0 | g2(i,i_parents(ii2))~=0
                        res=0;               
                        fprintf('there is an immorality in the 1st graph that is not in the 2nd graph\n');
                    end
                end
            end
        end
    end

    for i=1:num_nodes
        i_parents=find(g2(:,i)==-1);
        for ii1=1:(length(i_parents)-1)
            for ii2=(ii1+1):length(i_parents)
                if g2(i_parents(ii2),i_parents(ii1))==0 & g2(i_parents(ii1),i_parents(ii2))==0
                    if g1(i,i_parents(ii1))~=0 | g1(i,i_parents(ii2))~=0
                        res=0;               
                        fprintf('there is an immorality in the 2nd graph that is not in the 1st graph\n');
                    end
                end
            end
        end
    end
end
