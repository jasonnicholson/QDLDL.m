classdef QDLDLFactorization < handle
% QDLDLFactorization  Handle object for a sparse quasi-definite LDL factorization.

% Copyright (c) 2026 Jason H. Nicholson
% Copyright (c) Paul Goulart and QDLDL.jl contributors
% SPDX-License-Identifier: Apache-2.0
% Ported to MATLAB from QDLDL.jl (https://github.com/osqp/QDLDL.jl)
    properties
        perm
        iperm
        L
        L_unit   % unit lower triangular: L + I  (for fast MATLAB sparse \)
        Dinv_vec % Dinv as a plain vector   (avoids spdiags overhead per solve)
        Dinv
        workspace
        logical = false
    end

    methods
        function obj = QDLDLFactorization(perm, iperm, L, Dinv, workspace, logicalFlag)
            % QDLDLFactorization  Construct a QDLDLFactorization object.
            %
            %   OBJ = QDLDLFactorization(PERM, IPERM, L, DINV, WORKSPACE, LOGICALFLAG)
            %   stores the results of a sparse LDL^T factorization.  Typically created
            %   via QDLDLFactorization.fromMatrix rather than called directly.
            %
            %   Inputs:
            %     perm        - Fill-reducing permutation vector ([] if none applied)
            %     iperm       - Inverse of perm ([] if none applied)
            %     L           - Strict lower-triangular factor as a MATLAB sparse matrix
            %     Dinv        - Sparse diagonal matrix of inverse diagonal entries
            %     workspace   - Struct with CSC arrays and intermediate working storage
            %     logicalFlag - true if only a symbolic/logical factorization was performed
            obj.perm = perm;
            obj.iperm = iperm;
            obj.L = L;
            obj.L_unit   = L + speye(size(L, 1));
            obj.Dinv_vec = workspace.Dinv;
            obj.Dinv = Dinv;
            obj.workspace = workspace;
            obj.logical = logicalFlag;
        end

        function x = mldivide(obj, b)
            % mldivide  Overload the backslash operator (A \ b) for the factorization.
            %
            %   X = OBJ \ B  delegates to solve(B) to compute X = A^{-1} * B.
            x = obj.solve(b);
        end

        function x = solve(obj, b)
            % solve  Solve the linear system A*x = b using the LDL^T factorization.
            %
            %   X = OBJ.solve(B) applies the fill-reducing permutation, performs
            %   L*D*L^T back-substitution, and then undoes the permutation.
            x = obj.solveInPlace(b);
        end

        function x = solveInPlace(obj, b)
            % solveInPlace  Solve A*x = b with optional reordering and in-workspace storage.
            %
            %   X = OBJ.solveInPlace(B) computes the solution to A*x = b by:
            %     1. Permuting b as b_perm = b(perm)  (skipped if no permutation).
            %     2. Applying forward substitution   L  * y      = b_perm.
            %     3. Applying diagonal scaling        D^{-1} * z = y.
            %     4. Applying backward substitution  L' * x_perm = z.
            %     5. Unpermuting: x(perm) = x_perm.
            %
            %   Throws QDLDL:LogicalFactorization if only a symbolic factor is stored.
            if obj.logical
                error('QDLDL:LogicalFactorization', 'Cannot solve with logical factorization only.');
            end

            % Apply permutation
            if isempty(obj.perm)
                tmp = b;
            else
                tmp = b(obj.perm);
            end

            % L * D * L^T solve using MATLAB sparse triangular backsolve
            % (much faster than manual scalar loops for large problems)
            L1 = obj.L_unit;          % unit lower triangular (I + L_strict)
            tmp = L1 \ tmp;           % forward solve:  L * y = b
            tmp = tmp .* obj.Dinv_vec; % diagonal scale: y = Dinv .* y
            tmp = L1' \ tmp;          % backward solve: L' * x = y

            % Inverse permutation
            if isempty(obj.perm)
                % No reordering: result is already in the original ordering
                x = tmp;
            else
                x = b;             % allocate same size/type
                x(obj.perm) = tmp;
            end
        end

        function refactor(obj)
            % refactor  Recompute the numerical LDL^T factorization in place.
            %
            %   OBJ.refactor() reuses the existing symbolic structure (elimination tree,
            %   nonzero counts, column pointers) and performs a fresh numerical QDLDL
            %   factorization using whatever values are currently in workspace.triuA.
            %   The logical flag is cleared so the result can be used for solves.
            %   Call updateValues or scaleValues first to change matrix values.
            obj.logical = false;
            obj.workspace = QDLDLFactorization.factorWorkspace(obj.workspace, obj.logical);
            obj.refreshFactors();
        end

        function updateValues(obj, indices, values)
            % updateValues  Replace specific nonzero values in the stored upper-triangular A.
            %
            %   OBJ.updateValues(INDICES, VALUES) overwrites the nonzero entries at
            %   positions INDICES with VALUES.  If a fill-reducing permutation was applied,
            %   INDICES are mapped through AtoPAPt to the corresponding positions in the
            %   permuted matrix before writing.
            %
            %   Inputs:
            %     indices - Linear indices into the nonzero array of the original A
            %     values  - Replacement values (same length as indices)
            if isempty(obj.workspace.AtoPAPt)
                % No permutation: write directly into the CSC value array
                obj.workspace.triuA.Ax(indices) = values;
            else
                % Map original indices to permuted positions, then write
                obj.workspace.triuA.Ax(obj.workspace.AtoPAPt(indices)) = values;
            end
        end

        function scaleValues(obj, indices, scale)
            % scaleValues  Multiply specific nonzero values in the stored matrix by a scale factor.
            %
            %   OBJ.scaleValues(INDICES, SCALE) multiplies the nonzero entries at positions
            %   INDICES element-wise by SCALE.  If a fill-reducing permutation was applied,
            %   INDICES are mapped through AtoPAPt to the corresponding positions in the
            %   permuted matrix before scaling.
            %
            %   Inputs:
            %     indices - Linear indices into the nonzero array of the original A
            %     scale   - Multiplicative scale factors (same length as indices)
            if isempty(obj.workspace.AtoPAPt)
                % No permutation: scale directly in the CSC value array
                obj.workspace.triuA.Ax(indices) = obj.workspace.triuA.Ax(indices) .* scale;
            else
                % Map original indices to permuted positions before scaling
                mapped = obj.workspace.AtoPAPt(indices);
                obj.workspace.triuA.Ax(mapped) = obj.workspace.triuA.Ax(mapped) .* scale;
            end
        end

        function n = positiveInertia(obj)
            % positiveInertia  Return the number of positive entries in D.
            %
            %   N = OBJ.positiveInertia() returns the count of positive diagonal entries
            %   from the most recent numerical factorization.  By Sylvester's law of
            %   inertia this equals the number of positive eigenvalues of A.
            n = obj.workspace.positive_inertia;
        end

        function n = regularizedEntries(obj)
            % regularizedEntries  Return the number of diagonal entries that were regularized.
            %
            %   N = OBJ.regularizedEntries() returns the count of diagonal entries in D
            %   that were replaced by regularize_delta during the most recent numerical
            %   factorization (i.e., entries where Dsigns(i)*D(i) < regularize_eps).
            n = obj.workspace.regularize_count;
        end

        function refreshFactors(obj)
            % refreshFactors  Rebuild the public L and Dinv properties from workspace arrays.
            %
            %   OBJ.refreshFactors() converts the internal CSC representation of the
            %   strict lower-triangular factor L and the vector of inverse diagonal entries
            %   into MATLAB sparse matrices stored as OBJ.L and OBJ.Dinv.
            %   Called automatically by refactor() after each numerical factorization.

            % Convert CSC strict-lower arrays to a MATLAB sparse matrix
            obj.L = QDLDLFactorization.cscStrictLowerToSparse( ...
                obj.workspace.Ln, ...
                obj.workspace.Lp, ...
                obj.workspace.Li, ...
                obj.workspace.Lx);
            % Pre-compute unit lower triangular factor (I + L) for fast solve
            obj.L_unit   = obj.L + speye(obj.workspace.Ln);
            obj.Dinv_vec = obj.workspace.Dinv;   % plain vector, no spdiags
            obj.Dinv     = spdiags(obj.workspace.Dinv, 0, obj.workspace.Ln, obj.workspace.Ln);
        end
    end

    methods (Static)
        function obj = fromMatrix(A, options)
            % fromMatrix  Factory: compute the LDL^T factorization of a quasi-definite matrix.
            %
            %   OBJ = QDLDLFactorization.fromMatrix(A, OPTIONS) factorizes the symmetric
            %   quasi-definite matrix A using the QDLDL algorithm and returns a
            %   QDLDLFactorization object ready for solving or refactoring.
            %   Only the upper triangle of A is referenced.
            %
            %   Inputs:
            %     A       - Real symmetric (quasi-definite) matrix; dense or sparse.
            %     options - Struct with fields:
            %                 perm           : [] for no reordering, 'auto' for symamd
            %                                  (default), or an explicit permutation vector
            %                 Dsigns         : Length-n sign vector (+/-1) indicating the
            %                                  expected sign of each diagonal D entry
            %                                  (used for regularization)
            %                 logical        : true to perform only symbolic factorization
            %                 regularize_eps : Threshold triggering diagonal regularization
            %                 regularize_delta: Replacement value for regularized D entries
            %
            %   Output:
            %     obj - QDLDLFactorization object

            if ~issparse(A)
                % Promote dense input to sparse before extracting CSC data
                A = sparse(A);
            end

            [m, n] = size(A);
            if m ~= n
                error('QDLDL:DimensionMismatch', 'Input matrix must be square.');
            end

            % Ensure only the upper triangle is stored; triu strips lower entries
            if ~istriu(A)
                A = triu(A);
            else
                A = sparse(A);
            end

            % Determine the fill-reducing column permutation
            if ~isfield(options, 'perm') || (ischar(options.perm) && strcmp(options.perm, 'auto'))
                % Use symamd approximate minimum-degree ordering by default
                perm = symamd(A);
                perm = perm(:);
            elseif isempty(options.perm)
                % No reordering requested
                perm = [];
            else
                % Validate user-supplied permutation
                perm = options.perm(:);
                if numel(perm) ~= n || ~isequal(sort(double(perm(:))), (1:n)')
                    error('QDLDL:InvalidPermutation', 'perm must be a permutation of 1:n.');
                end
            end

            % Compute the inverse permutation (maps new indices back to original)
            iperm = [];
            if ~isempty(perm)
                iperm = QDLDLFactorization.inversePerm(perm);
            end

            % Convert the upper-triangular sparse matrix to an internal CSC struct
            cscA = QDLDLFactorization.sparseToCSCUpper(A);
            if ~isempty(perm)
                % Symmetrically permute A -> P*A*P^T and record the index map AtoPAPt
                [cscA, AtoPAPt] = QDLDLFactorization.permuteSymmetricCSC(cscA, iperm);
            else
                AtoPAPt = [];
            end

            % Build the permuted sign vector for diagonal regularization (if provided)
            mySigns = [];
            if isfield(options, 'Dsigns') && ~isempty(options.Dsigns)
                signs = options.Dsigns(:);
                if numel(signs) ~= n
                    error('QDLDL:DimensionMismatch', 'Dsigns must have one entry per row/column.');
                end
                if isempty(perm)
                    mySigns = signs;
                else
                    % Reorder signs to match the permuted row/column numbering
                    mySigns = signs(perm);
                end
            end

            % Allocate the workspace (etree, Lnz counts, and CSC arrays for L and D)
            workspace = QDLDLFactorization.initializeWorkspace( ...
                cscA, ...
                AtoPAPt, ...
                mySigns, ...
                options.regularize_eps, ...
                options.regularize_delta);

            % Run the numerical (or symbolic) LDL^T factorization
            workspace = QDLDLFactorization.factorWorkspace(workspace, options.logical);

            % Convert internal CSC arrays to public MATLAB sparse matrix objects
            L = QDLDLFactorization.cscStrictLowerToSparse(workspace.Ln, workspace.Lp, workspace.Li, workspace.Lx);
            Dinv = spdiags(workspace.Dinv, 0, workspace.Ln, workspace.Ln);

            obj = QDLDLFactorization(perm, iperm, L, Dinv, workspace, options.logical);
        end
    end

    methods (Static, Access = private)
        function ip = inversePerm(p)
            % inversePerm  Compute the inverse of a permutation vector.
            %
            %   IP = inversePerm(P) returns IP such that IP(P(i)) = i for all i.
            %   P must be a permutation of 1:n.
            ip = zeros(size(p));
            ip(double(p)) = 1:numel(p);
        end

        function csc = sparseToCSCUpper(S)
            % sparseToCSCUpper  Convert an upper-triangular MATLAB sparse matrix to a CSC struct.
            %
            %   CSC = sparseToCSCUpper(S) extracts the nonzero pattern and values from the
            %   upper-triangular sparse matrix S and returns a struct with fields:
            %     n  - Matrix dimension (S is assumed square)
            %     Ap - Column pointer array of length n+1 (1-based)
            %     Ai - Row index array for each nonzero entry
            %     Ax - Nonzero value array
            n = size(S, 1);
            [rows, cols, vals] = find(S);

            % Count nonzeros per column to build the CSC column pointer array
            counts = accumarray(cols, 1, [n, 1]);
            Ap = [1; cumsum(counts) + 1];

            csc = struct();
            csc.n = n;
            csc.Ap = Ap;
            csc.Ai = rows;
            csc.Ax = vals;
        end

        function [P, AtoPAPt] = permuteSymmetricCSC(A, iperm)
            % permuteSymmetricCSC  Symmetrically permute a CSC upper-triangular matrix.
            %
            %   [P, ATOPAP T] = permuteSymmetricCSC(A, IPERM) computes the upper-triangular
            %   part of P_mat * A_mat * P_mat^T, where P_mat is the permutation matrix
            %   derived from the inverse permutation IPERM.
            %
            %   Inputs:
            %     A     - CSC struct for the upper triangle of the original matrix
            %     iperm - Inverse permutation vector (iperm(old_index) = new_index)
            %
            %   Outputs:
            %     P       - CSC struct for the upper triangle of P_mat*A_mat*P_mat^T
            %     AtoPAPt - Index map: AtoPAPt(k) is the position in P.Ax corresponding
            %               to the k-th nonzero in A.Ax.  Used by updateValues/scaleValues.
            n = A.n;
            Ap = A.Ap;
            Ai = A.Ai;
            Ax = A.Ax;

            % --- Pass 1: count nonzeros per column in the permuted upper triangle ---
            numEntries = zeros(n, 1);
            for colA = 1:n
                colP = iperm(colA);
                for idx = Ap(colA):(Ap(colA + 1) - 1)
                    rowA = Ai(idx);
                    if rowA <= colA
                        % Map both row and column through the inverse permutation
                        rowP = iperm(rowA);
                        % The permuted entry lives in the upper triangle under the larger index
                        colIdx = max(rowP, colP);
                        numEntries(colIdx) = numEntries(colIdx) + 1;
                    end
                end
            end

            % Build the CSC column pointer array from per-column counts
            Pc = zeros(n + 1, 1);
            Pc(1) = 1;
            for k = 1:n
                Pc(k + 1) = Pc(k) + numEntries(k);
                % Reuse numEntries as a running insertion pointer for the fill pass
                numEntries(k) = Pc(k);
            end

            % --- Pass 2: fill row indices, values, and the AtoPAPt index map ---
            rowStarts = numEntries;  % Current insertion position per column
            nnzA = numel(Ax);
            Pr = zeros(nnzA, 1);
            Pv = zeros(nnzA, 1, class(Ax));
            AtoPAPt = zeros(nnzA, 1);

            for colA = 1:n
                colP = iperm(colA);
                for rowAIdx = Ap(colA):(Ap(colA + 1) - 1)
                    rowA = Ai(rowAIdx);
                    if rowA <= colA
                        rowP = iperm(rowA);
                        % Determine the permuted column (larger re-indexed coordinate)
                        colIdx = max(colP, rowP);
                        rowPIdx = rowStarts(colIdx);

                        % Store permuted row index and value
                        Pr(rowPIdx) = min(colP, rowP);
                        Pv(rowPIdx) = Ax(rowAIdx);
                        % Record where this original nonzero landed in the permuted array
                        AtoPAPt(rowAIdx) = rowPIdx;

                        rowStarts(colIdx) = rowStarts(colIdx) + 1;
                    end
                end
            end

            P = struct();
            P.n = n;
            P.Ap = Pc;
            P.Ai = Pr;
            P.Ax = Pv;
        end

        function workspace = initializeWorkspace(cscA, AtoPAPt, Dsigns, regularizeEps, regularizeDelta)
            % initializeWorkspace  Allocate and populate the factorization workspace struct.
            %
            %   WORKSPACE = initializeWorkspace(CSCA, ATOPAP T, DSIGNS, EPS, DELTA)
            %   computes the symbolic elimination tree and per-column nonzero counts for
            %   the L factor, then allocates all arrays needed by the numerical factorization.
            %
            %   Inputs:
            %     cscA           - CSC struct for the (permuted) upper-triangular matrix
            %     AtoPAPt        - Index map from original to permuted nonzeros ([] if none)
            %     Dsigns         - Sign vector for diagonal regularization ([] to disable)
            %     regularizeEps  - Threshold: regularize entry i if Dsigns(i)*D(i) < eps
            %     regularizeDelta- Replacement value for regularized diagonal entries
            %
            %   Output:
            %     workspace - Struct with all CSC arrays and scratch buffers required
            %                 by qdldlFactor and qdldlSolve
            n = cscA.n;

            % Pre-allocate symbolic working arrays
            etree = zeros(n, 1);       % Elimination tree parent pointers
            Lnz = zeros(n, 1);         % Nonzero count per column of L
            iwork = zeros(3 * n, 1);   % General integer scratch buffer (reserved)
            bwork = false(n, 1);       % Boolean marker array (reused each factor step)
            fwork = zeros(n, 1);       % Floating-point scratch buffer (reused during solve)

            % Compute the elimination tree and nonzero counts for L
            [sumLnz, Lnz, etree] = QDLDLFactorization.qdldlEtree(n, cscA.Ap, cscA.Ai, Lnz, etree);
            if sumLnz < 0
                error('QDLDL:InvalidInput', 'Input matrix is not upper triangular or has an empty column.');
            end

            workspace = struct();

            % Symbolic structure arrays
            workspace.etree = etree;   % Parent of each node in the elimination tree
            workspace.Lnz = Lnz;       % Number of nonzeros in each column of L
            workspace.iwork = iwork;   % Integer scratch (reserved for future use)
            workspace.bwork = bwork;   % Boolean scratch (reused each factor step)
            workspace.fwork = fwork;   % Float scratch (reused during solves)

            % CSC representation of the L factor (filled during factorization)
            workspace.Ln = n;
            workspace.Lp = zeros(n + 1, 1);                      % Column pointers for L
            workspace.Li = zeros(sumLnz, 1);                     % Row indices for L nonzeros
            workspace.Lx = zeros(sumLnz, 1, class(cscA.Ax));    % Values of L nonzeros

            % Diagonal arrays
            workspace.D    = zeros(n, 1, class(cscA.Ax));  % Diagonal of D
            workspace.Dinv = zeros(n, 1, class(cscA.Ax));  % Element-wise inverse of D
            workspace.positive_inertia = -1;  % Filled after numerical factorization

            % References to the permuted input matrix and index map
            workspace.triuA   = cscA;
            workspace.AtoPAPt = AtoPAPt;

            % Regularization parameters
            workspace.Dsigns          = Dsigns;
            workspace.regularize_eps  = regularizeEps;
            workspace.regularize_delta = regularizeDelta;
            workspace.regularize_count = 0;  % Reset regularization counter
        end

        function workspace = factorWorkspace(workspace, logicalFactor)
            % factorWorkspace  Run the QDLDL numerical (or logical) factorization.
            %
            %   WORKSPACE = factorWorkspace(WORKSPACE, LOGICALFACTOR) performs the LDL^T
            %   factorization of the matrix stored in workspace.triuA.  If LOGICALFACTOR
            %   is true, all L values and D entries are set to 1 (symbolic-only mode).
            %   Otherwise, the full numerical factorization is executed and regularization
            %   may adjust diagonal entries that violate the quasi-definite sign conditions.
            %
            %   Throws QDLDL:NonQuasiDefinite if a zero appears on the D diagonal.
            if logicalFactor
                % Symbolic mode: fill L and D with trivial unit values
                workspace.Lx(:) = 1;
                workspace.D(:) = 1;
                workspace.Dinv(:) = 1;
            end

            % Execute the core QDLDL factorization and count positive D entries
            [posDCount, workspace] = QDLDLFactorization.qdldlFactor(workspace, logicalFactor);
            if posDCount < 0
                error('QDLDL:NonQuasiDefinite', 'Zero entry in D (matrix is not quasidefinite).');
            end

            % Store inertia count for later inspection via positiveInertia()
            workspace.positive_inertia = posDCount;
        end

        function [sumLnz, Lnz, etree] = qdldlEtree(n, Ap, Ai, Lnz, etree)
            % qdldlEtree  Compute the elimination tree and per-column L nonzero counts.
            %
            %   [SUMLNZ, LNZ, ETREE] = qdldlEtree(N, AP, AI, LNZ, ETREE) computes the
            %   elimination tree of the symmetric matrix whose upper triangle is given in
            %   CSC form (Ap, Ai), and counts the number of nonzeros in each column of L.
            %
            %   Inputs:
            %     n     - Matrix dimension
            %     Ap    - CSC column pointer array (length n+1, 1-based)
            %     Ai    - CSC row index array
            %     Lnz   - Pre-allocated output array of length n (overwritten)
            %     etree - Pre-allocated parent array of length n (overwritten)
            %
            %   Outputs:
            %     sumLnz - Total number of nonzeros in L (-1 if input is invalid)
            %     Lnz    - Number of nonzeros in column i of L
            %     etree  - etree(i) is the parent of node i in the elimination tree
            %              (-1 if i is a root node)

            % work(i) records the most recent column j that visited row i during the walk
            work = zeros(n, 1);

            % Initialize node markers, nonzero counts, and parent pointers
            for i = 1:n
                work(i) = 0;
                Lnz(i) = 0;
                etree(i) = -1;

                % An empty diagonal column makes the matrix structurally singular; abort
                if Ap(i) == Ap(i + 1)
                    sumLnz = -1;
                    return;
                end
            end

            % Walk each column j and climb the elimination tree to count fill-in
            for j = 1:n
                work(j) = j;  % Mark j as "claimed" by column j
                for p = Ap(j):(Ap(j + 1) - 1)
                    i = Ai(p);
                    if i > j
                        % Row index exceeds column index: matrix is not upper triangular
                        sumLnz = -1;
                        return;
                    end

                    % Walk ancestors of row i until we reach a node already claimed by j
                    while work(i) ~= j
                        if etree(i) == -1
                            % Discovered a new parent edge in the elimination tree
                            etree(i) = j;
                        end
                        Lnz(i) = Lnz(i) + 1;  % Column i gains one more nonzero in L
                        work(i) = j;            % Mark i as visited during column j's pass
                        i = etree(i);           % Climb to the parent node
                    end
                end
            end

            sumLnz = sum(Lnz);
        end

        function [positiveValuesInD, workspace] = qdldlFactor(workspace, logicalFactor)
            % qdldlFactor  Core QDLDL numerical LDL^T factorization.
            %
            %   [POSDCOUNT, WORKSPACE] = qdldlFactor(WORKSPACE, LOGICALFACTOR) computes
            %   the LDL^T factorization of the quasi-definite matrix stored in
            %   workspace.triuA (CSC upper triangle).  The factor L is stored back into
            %   workspace.Lp / workspace.Li / workspace.Lx, and the diagonal D is stored
            %   in workspace.D / workspace.Dinv.
            %
            %   The algorithm processes columns left-to-right, maintaining a sparse
            %   accumulator yVals (indexed by yIdx) for the pivot column, and walks the
            %   elimination tree to identify contributing ancestor columns.
            %
            %   Inputs:
            %     workspace     - Initialized workspace struct (from initializeWorkspace)
            %     logicalFactor - true to skip numerical fill and keep unit L/D values
            %
            %   Outputs:
            %     positiveValuesInD - Number of positive D entries (-1 if a zero is found)
            %     workspace         - Updated struct with L and D arrays filled in

            % Unpack input matrix arrays for readability
            n  = workspace.triuA.n;
            Ap = workspace.triuA.Ap;
            Ai = workspace.triuA.Ai;
            Ax = workspace.triuA.Ax;

            % Unpack L and D working arrays from workspace
            Lp     = workspace.Lp;
            Li     = workspace.Li;
            Lx     = workspace.Lx;
            D      = workspace.D;
            Dinv   = workspace.Dinv;
            Lnz    = workspace.Lnz;
            etree  = workspace.etree;
            yMarkers = workspace.bwork;  % Boolean: is column c in the active set?
            yVals    = workspace.fwork;  % Sparse accumulator values for current column

            % Scratch buffers for the sparse column elimination loop
            yIdx            = zeros(n, 1);  % Ordered list of active column indices
            elimBuffer      = zeros(n, 1);  % Stack for etree ancestor traversal
            LNextSpaceInCol = zeros(n, 1);  % Next available slot in Li/Lx per column

            positiveValuesInD = 0;
            workspace.regularize_count = 0;

            % Build CSC column pointers for L and initialize per-column scratch state
            Lp(1) = 1;
            for i = 1:n
                Lp(i + 1) = Lp(i) + Lnz(i);  % Column i starts at offset Lp(i)
                yMarkers(i) = false;
                yVals(i) = 0.0;
                D(i) = 0.0;
                LNextSpaceInCol(i) = Lp(i);   % Writing cursor for column i of L
            end

            if logicalFactor
                % Symbolic mode: skip numerical work, keep unit values set by caller
                D(:) = 1;
                Dinv(:) = 1;
            else
                % Bootstrap column 1: D(1) is the (1,1) diagonal entry of A
                D(1) = Ax(1);

                % Apply regularization if D(1) has the wrong sign or is too small
                if ~isempty(workspace.Dsigns) && workspace.Dsigns(1) * D(1) < workspace.regularize_eps
                    D(1) = workspace.regularize_delta * workspace.Dsigns(1);
                    workspace.regularize_count = workspace.regularize_count + 1;
                end

                if D(1) == 0.0
                    positiveValuesInD = -1;  % Zero pivot: matrix is not quasi-definite
                    return;
                end
                if D(1) > 0.0
                    positiveValuesInD = positiveValuesInD + 1;
                end
                Dinv(1) = 1 / D(1);
            end

            % --- Main column factorization loop (columns 2 through n) ---
            for k = 2:n
                nnzY = 0;  % Number of active entries in the sparse accumulator

                % --- Scatter column k of A into the sparse accumulator yVals ---
                for i = Ap(k):(Ap(k + 1) - 1)
                    bidx = Ai(i);  % Row index of this nonzero in column k

                    if bidx == k
                        % Diagonal entry: initialize D(k) from A(k,k)
                        D(k) = Ax(i);
                        continue;
                    end

                    % Off-diagonal upper entry: load into the sparse accumulator
                    yVals(bidx) = Ax(i);
                    nextIdx = bidx;

                    if ~yMarkers(nextIdx)
                        % First visit to this column: walk up the elimination tree
                        % to collect all ancestors that contribute to column k
                        yMarkers(nextIdx) = true;
                        elimBuffer(1) = nextIdx;
                        nnzE = 1;

                        nextIdx = etree(bidx);
                        while nextIdx ~= -1 && nextIdx < k
                            if yMarkers(nextIdx)
                                break;  % Already queued; stop climbing
                            end

                            yMarkers(nextIdx) = true;
                            nnzE = nnzE + 1;
                            elimBuffer(nnzE) = nextIdx;
                            nextIdx = etree(nextIdx);
                        end

                        % Reverse the ancestor stack into yIdx (root-to-leaf order)
                        while nnzE ~= 0
                            nnzY = nnzY + 1;
                            yIdx(nnzY) = elimBuffer(nnzE);
                            nnzE = nnzE - 1;
                        end
                    end
                end

                % --- Elimination loop: update D(k) and fill column k of L ---
                % Process active columns in reverse order (deepest ancestor first)
                for i = nnzY:-1:1
                    cidx   = yIdx(i);                % Contributing column index
                    tmpIdx = LNextSpaceInCol(cidx);  % Next free slot in column cidx of L

                    if ~logicalFactor
                        yVals_cidx = yVals(cidx);
                        % Subtract contributions of already-computed L entries:
                        % y <- y - L(cidx, :) * y(cidx)
                        for j = Lp(cidx):(tmpIdx - 1)
                            yVals(Li(j)) = yVals(Li(j)) - Lx(j) * yVals_cidx;
                        end

                        % L(k, cidx) = y(cidx) / D(cidx)
                        Lx(tmpIdx) = yVals_cidx * Dinv(cidx);
                        % Schur complement update: D(k) -= y(cidx) * L(k, cidx)
                        D(k) = D(k) - yVals_cidx * Lx(tmpIdx);
                    end

                    % Store the row index of this L entry and advance the write cursor
                    Li(tmpIdx) = k;
                    LNextSpaceInCol(cidx) = LNextSpaceInCol(cidx) + 1;

                    % Clear the accumulator entry and marker for column cidx
                    yVals(cidx) = 0.0;
                    yMarkers(cidx) = false;
                end

                if ~logicalFactor
                    % Apply regularization to D(k) if it has the wrong sign or is too small
                    if ~isempty(workspace.Dsigns) && workspace.Dsigns(k) * D(k) < workspace.regularize_eps
                        D(k) = workspace.regularize_delta * workspace.Dsigns(k);
                        workspace.regularize_count = workspace.regularize_count + 1;
                    end

                    if D(k) == 0.0
                        positiveValuesInD = -1;  % Zero pivot: abort
                        return;
                    end
                    if D(k) > 0.0
                        positiveValuesInD = positiveValuesInD + 1;
                    end

                    Dinv(k) = 1 / D(k);
                end
            end

            % Write updated arrays back to the workspace struct
            workspace.Lp    = Lp;
            workspace.Li    = Li;
            workspace.Lx    = Lx;
            workspace.D     = D;
            workspace.Dinv  = Dinv;
            workspace.bwork = yMarkers;
            workspace.fwork = yVals;
        end

        function x = qdldlLsolve(n, Lp, Li, Lx, x)
            % qdldlLsolve  Forward substitution: solve L*x = b in place.
            %
            %   X = qdldlLsolve(N, LP, LI, LX, X) solves the lower-triangular system
            %   L * x = b, where L is the unit lower-triangular LDL factor stored in
            %   CSC form (Lp, Li, Lx) with strictly lower-triangular structure.
            %   The solution overwrites the input vector X (which enters as b).
            for i = 1:n
                % Propagate the effect of x(i) onto all rows below it in column i
                for j = Lp(i):(Lp(i + 1) - 1)
                    x(Li(j)) = x(Li(j)) - Lx(j) * x(i);
                end
            end
        end

        function x = qdldlLtsolve(n, Lp, Li, Lx, x)
            % qdldlLtsolve  Backward substitution: solve L'*x = b in place.
            %
            %   X = qdldlLtsolve(N, LP, LI, LX, X) solves the upper-triangular system
            %   L' * x = b, where L is the unit lower-triangular LDL factor stored in
            %   CSC form (Lp, Li, Lx).  The solution overwrites the input vector X.
            for i = n:-1:1
                % Gather contributions from rows below i back into x(i)
                for j = Lp(i):(Lp(i + 1) - 1)
                    x(i) = x(i) - Lx(j) * x(Li(j));
                end
            end
        end

        function b = qdldlSolve(n, Lp, Li, Lx, Dinv, b)
            % qdldlSolve  Solve the system (L*D*L') * x = b.
            %
            %   B = qdldlSolve(N, LP, LI, LX, DINV, B) solves the factored system by
            %   applying the three-phase LDL^T back-substitution:
            %     1. Forward  substitution: solve L   * y = b  (overwrites b)
            %     2. Diagonal scaling:             D^{-1}*z = y  (element-wise multiply)
            %     3. Backward substitution: solve L'  * x = z  (overwrites b)
            b = QDLDLFactorization.qdldlLsolve(n, Lp, Li, Lx, b);   % Step 1: L * y = b
            b = b .* Dinv;                                            % Step 2: z = D^{-1} * y
            b = QDLDLFactorization.qdldlLtsolve(n, Lp, Li, Lx, b);  % Step 3: L' * x = z
        end

        function L = cscStrictLowerToSparse(n, Lp, Li, Lx)
            % cscStrictLowerToSparse  Convert a CSC strict-lower-triangular factor to MATLAB sparse.
            %
            %   L = cscStrictLowerToSparse(N, LP, LI, LX) converts the strictly
            %   lower-triangular L factor stored in CSC form (Lp, Li, Lx) into a standard
            %   MATLAB sparse matrix of size n-by-n.
            %
            %   Inputs:
            %     n  - Matrix dimension
            %     Lp - CSC column pointer array (length n+1, 1-based)
            %     Li - Row index array for L nonzeros
            %     Lx - Value array for L nonzeros
            %
            %   Output:
            %     L - n-by-n MATLAB sparse lower-triangular matrix
            nnzL = Lp(end) - 1;  % Total number of nonzeros in L
            if nnzL == 0
                L = sparse(n, n);
                return;
            end
            % Vectorised: build column-index vector from CSC pointer differences
            lens = diff(Lp);   % number of nonzeros per column
            cols = repelem((1:n)', lens);
            L = sparse(Li, cols, Lx, n, n);
        end
    end
end
