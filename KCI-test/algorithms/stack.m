function v = stack(M)
% stack the matrix M into the vector v
% Copyright (c) 2010-2011  ...
% All rights reserved.  See the file COPYING for license terms.
[n,t] = size(M);
v = zeros(n*t,1);

for i=1:t
    v((i-1)*n+1:i*n) = M(:,i);
end
