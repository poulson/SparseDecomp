%
%  Copyright (c) 2014, Jack Poulson
%
%  All rights reserved.
%
%  This file is part of SparseDecomp and is under the BSD 2-Clause License, 
%  which can be found in the LICENSE file in the root directory, or at 
%  http://opensource.org/licenses/BSD-2-Clause
%
function [lambda,SList] = hierPALM(A,SList,projPair,numIter)

% This algorithm is based upon the paper
%   Luc Le Magoarou and Remi Gribonval, 
%   "Learning computationally efficient dictionaries and their implementation
%   as fast transforms", arXiv:1406.5388v2, 30 Jun 2014.
%
% Question: What is the canonical initialization of the set of sparse factors?
%
% NOTE: projPair{1} should be a 'strong' projection, while projPair{2}
%       should be a 'weak' projection, which usually means that the former
%       should keep a smaller number of nonzero entries

numMats=length(SList);

% Ensure that A is square
[m,n] = size(A);
if( m ~= n ), error('Matrix must be square'); end

R=A;
for k=1:numMats-1,
  TPair{1}=SList{k};   % Proper initialization???
  TPair{2}=SList{k+1}; % Proper initialization???
  [lambda,TPair] = PALM(R,TPair,projPair,numIter);
  SList{k}=lambda*TPair{1};
  R=TPair{2};

  for j=1:k, SListLeft{j}=SList{j}; projListLeft{j}=projPair{1}; end
  SListLeft{k+1}=R; projListLeft{k+1}=projPair{2};
  [lambda,SListLeft] = PALM(A,SListLeft,projListLeft,numIter);
  for j=1:k+1, SList{j}=SListLeft{j}; end
end

return
