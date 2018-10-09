function M = setdiag(M, v)
% SETDIAG Set the diagonal of a matrix to a specified scalar/vector.
% M = set_diag(M, v)
%Copyright (C) 
%		2010 Robert Tillman 
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

n = length(M);
if length(v)==1
  v = repmat(v, 1, n);
end

% e.g., for 3x3 matrix,  elements are numbered
% 1 4 7 
% 2 5 8 
% 3 6 9
% so diagnoal = [1 5 9]


J = 1:n+1:n^2;
M(J) = v;

%M = triu(M,1) + tril(M,-1) + diag(v);

