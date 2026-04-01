function p = positive_inertia(F)
%POSITIVE_INERTIA Number of positive diagonal entries in D.

% Copyright (c) 2026 Jason H. Nicholson
% Copyright (c) Paul Goulart and QDLDL.jl contributors
% SPDX-License-Identifier: Apache-2.0
% Ported to MATLAB from QDLDL.jl (https://github.com/osqp/QDLDL.jl)
    p = F.positiveInertia();
end
