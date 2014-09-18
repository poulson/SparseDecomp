%
%  Copyright (c) 2014, Jack Poulson
%
%  All rights reserved.
%
%  This file is part of SparseDecomp and is under the BSD 2-Clause License, 
%  which can be found in the LICENSE file in the root directory, or at 
%  http://opensource.org/licenses/BSD-2-Clause
%

randomInit=false;
numMats=2;
N=2^4;
numIter=10000;

strongSparse=3*N;
weakSparse=N^2/2;

A=hadamard(N); A=A/norm(A,'fro');

% NOTE: Initializing with the identity is problematic since
%       the identity matrix is orthogonal to the Hadamard matrix
%       with respect to the Hilbert-Schmidt inner product
if randomInit,
  for j=1:numMats,
    SList{j}=randn(N,N); SList{j}=SList{j}/norm(SList{j},'fro');
  end
else
  SList{1}=randn(N,N); SList{1}=SList{1}/norm(SList{1},'fro');
  SList{2}=A;
end

projStrong=@(A) keepTop(A,strongSparse);
projWeak=@(A) keepTop(A,weakSparse);
projPair{1}=projStrong;
projPair{2}=projWeak;

[lambda,SList]=PALM(A,SList,projPair,numIter);
ATild=lambda*eye(N,N); for j=1:numMats, ATild=ATild*SList{j}; end

figure(1); imagesc(A); title('A'); colorbar;
figure(2); imagesc(ATild); title('ATild'); colorbar;
figure(3); imagesc(A-ATild); title('A-ATild'); colorbar;
for j=1:numMats,
  figure(3+j); imagesc(SList{j}); title(sprintf('S %d',j)); colorbar;
end
