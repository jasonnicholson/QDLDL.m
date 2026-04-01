function x = solve_inplace(F, b)
%SOLVE_INPLACE Solve Ax=b and return the overwritten vector.

% Copyright (c) 2026 Jason H. Nicholson
% Copyright (c) Paul Goulart and QDLDL.jl contributors
% SPDX-License-Identifier: Apache-2.0
% Ported to MATLAB from QDLDL.jl (https://github.com/osqp/QDLDL.jl)
    x = F.solveInPlace(b);
end
