%
%  Copyright (c) 2014, Jack Poulson
%
%  All rights reserved.
%
%  This file is part of SparseDecomp and is under the BSD 2-Clause License, 
%  which can be found in the LICENSE file in the root directory, or at 
%  http://opensource.org/licenses/BSD-2-Clause
%
function A = keepTop(A,numEntries)
[m,n]=size(A);

% Find the magnitude of the numEntries+1'th largest entry
sortedSizes=sort(abs(vec(A)),'descend');
threshold=sortedSizes(min(numEntries+1,m*n));

% Threshold every value in A less than or equal to the cutoff, so that at 
% most numEntries values remain nonzero
A=arrayfun(@(alpha) hardThresh(alpha,threshold),A);

% Rescale A so that its Frobenius norm is one
A=A/norm(A,'fro');

return
