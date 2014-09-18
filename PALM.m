%
%  Copyright (c) 2014, Jack Poulson
%
%  All rights reserved.
%
%  This file is part of SparseDecomp and is under the BSD 2-Clause License, 
%  which can be found in the LICENSE file in the root directory, or at 
%  http://opensource.org/licenses/BSD-2-Clause
%
function [lambda,SList] = PALM(A,SList,projList,numIter)

% This algorithm is based upon the paper
%   Luc Le Magoarou and Remi Gribonval, 
%   "Learning computationally efficient dictionaries and their implementation
%   as fast transforms", arXiv:1406.5388v2, 30 Jun 2014.
%
% Question: What is the canonical initialization of the set of sparse factors?
%

numMats=length(SList);

% Ensure that A is square
[m,n] = size(A);
if( m ~= n ), error('Matrix must be square'); end

% Initialize the scaling factor
AHat = eye(n,n); for k=1:numMats, AHat=AHat*SList{k}; end
lambda = (vec(A)'*vec(AHat))/(vec(AHat)'*vec(AHat));

for i=1:numIter,
  for j=1:numMats,
    L = eye(n,n); for k=1:j-1,       L = L*SList{k}; end
    R = eye(n,n); for k=j+1:numMats, R = R*SList{k}; end

    % Select an inverse step length greater than the Lipschitz coefficient
    % -------------------------------------------------------------------- 
    % NOTE: Since frobNorm(L)=frobNorm(R)=1, norm(L), norm(R) <= 1
    %cj = 2*(lambda*norm(L)*norm(R))^2;
    cj = 1.5*lambda^2; 
    
    % Take a step in the negative gradient direction
    % ----------------------------------------------
    SList{j} = SList{j} - (lambda/cj)*L'*(lambda*L*SList{j}*R-A)*R';
    
    % Project the updated state
    % -------------------------
    SList{j} = projList{j}(SList{j});
  end
  AHat = eye(n,n); for k=1:numMats, AHat=AHat*SList{k}; end
  lambda = (vec(A)'*vec(AHat))/(vec(AHat)'*vec(AHat));
end

return
