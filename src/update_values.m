function update_values(F, indices, values)
%UPDATE_VALUES Overwrite selected upper-triangular nonzeros in-place.

% Copyright (c) 2026 Jason H. Nicholson
% Copyright (c) Paul Goulart and QDLDL.jl contributors
% SPDX-License-Identifier: Apache-2.0
% Ported to MATLAB from QDLDL.jl (https://github.com/osqp/QDLDL.jl)
    F.updateValues(indices, values);
end
