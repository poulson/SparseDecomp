%
%  Copyright (c) 2014, Jack Poulson
%
%  All rights reserved.
%
%  This file is part of SparseDecomp and is under the BSD 2-Clause License, 
%  which can be found in the LICENSE file in the root directory, or at 
%  http://opensource.org/licenses/BSD-2-Clause
%
function gamma = hardThresh(alpha,threshold)

if abs(alpha) > threshold,
  gamma = alpha;
else
  gamma = 0;
end

return
