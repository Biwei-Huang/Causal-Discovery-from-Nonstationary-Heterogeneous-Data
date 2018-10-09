%Copyright (C) 
%               1997-2002 Kevin Murphy
%		2010-2011 Robert Tillman 
%
%    This file is part of pc.
%
%    discrete_anm is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    discrete_anm is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with discrete_anm.  If not, see <http://www.gnu.org/licenses/>.    


% meeks rules - adapted from version in BNT
function pdag = meeks(pdag,G)

n = size(pdag,1);
old_pdag=zeros(n,n);
while ~isequal(pdag, old_pdag)
  old_pdag = pdag;
  % rule 1
  [A,B] = find(pdag==-1); % a -> b
  for i=1:length(A)
    a = A(i); b = B(i);
    C = find(pdag(b,:)==1 & G(a,:)==0); % all nodes adj to b but not a
    if ~isempty(C)
      pdag(b,C) = -1; pdag(C,b) = 0;
      %fprintf('rule 1: a=%d->b=%d and b=%d-c=%d implies %d->%d\n', a, b, b, C, b, C);
    end
  end
  % rule 2
  [A,B] = find(pdag==1); % unoriented a-b edge
  for i=1:length(A)
    a = A(i); b = B(i);
    if any( (pdag(a,:)==-1) & (pdag(:,b)==-1)' );
      pdag(a,b) = -1; pdag(b,a) = 0;
      %fprintf('rule 2: %d -> %d\n', a, b);
    end
  end
  % rule 3
  [A,B] = find(pdag==1); % a-b
  for i=1:length(A)
    a = A(i); b = B(i);
    C = find( (G(a,:)==1) & (pdag(:,b)==-1)' );
    % C contains nodes c s.t. a-c->ba
    G2 = setdiag(G(C, C), 1);
    if any(G2(:)==0) % there are 2 different non adjacent elements of C
      pdag(a,b) = -1; pdag(b,a) = 0;
      %fprintf('rule 3: %d -> %d\n', a, b);
    end
  end
  % rule 4
  [A, B] = find(pdag==1); % a-b
  for i=1:length(A)
    a = A(i); b = B(i);
    C = find((pdag(:,b)==-1) & (G(:,a)==1));
    for j=1:length(C)
      c = C(j); % c -> b and c - a
      D = find((pdag(:,c)==-1) & (pdag(:,a)==1)); % d -> c and d - a
      if (length(D)>0)
         pdag(a,b) = -1;
         %pdag(b,a) = -1; % It is a bug in the original version!!!!!!!!!!
         pdag(b,a) = 0;
      end
    end
  end
end
