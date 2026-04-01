% QDLDL  Quasi-definite LDL factorization.
%
%   qdldl               - Factor a sparse quasi-definite matrix using LDL.
%   solve               - Solve Ax=b using a QDLDLFactorization.
%   solve_inplace       - Solve Ax=b, overwriting b with the solution.
%   refactor_qdldl      - Recompute numerical factorization after data updates.
%   update_values       - Overwrite selected upper-triangular nonzeros in-place.
%   scale_values        - Scale selected upper-triangular nonzeros in-place.
%   positive_inertia    - Number of positive diagonal entries in D.
%   regularized_entries - Number of regularized diagonal entries in D.
%
%   The main entry point is QDLDL.  All other functions operate on the
%   QDLDLFactorization handle object returned by QDLDL.
%
%   Example:
%       F = qdldl(A);
%       x = F \ b;
%
%   See also QDLDLFactorization.
