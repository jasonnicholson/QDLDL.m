function refactor_qdldl(F)
%REFACTOR_QDLDL Recompute the numerical factorization after data updates.

% Copyright (c) 2026 Jason H. Nicholson
% Copyright (c) Paul Goulart and QDLDL.jl contributors
% SPDX-License-Identifier: Apache-2.0
% Ported to MATLAB from QDLDL.jl (https://github.com/osqp/QDLDL.jl)
    F.refactor();
end
