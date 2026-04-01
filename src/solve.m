function x = solve(F, b)
%SOLVE Solve Ax=b using a QDLDLFactorization.

% Copyright (c) 2026 Jason H. Nicholson
% Copyright (c) Paul Goulart and QDLDL.jl contributors
% SPDX-License-Identifier: Apache-2.0
% Ported to MATLAB from QDLDL.jl (https://github.com/osqp/QDLDL.jl)
    x = F.solve(b);
end
