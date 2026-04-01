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
        Dinv
        workspace
        logical = false
    end

    methods
        function obj = QDLDLFactorization(perm, iperm, L, Dinv, workspace, logicalFlag)
            obj.perm = perm;
            obj.iperm = iperm;
            obj.L = L;
            obj.Dinv = Dinv;
            obj.workspace = workspace;
            obj.logical = logicalFlag;
        end

        function x = mldivide(obj, b)
            x = obj.solve(b);
        end

        function x = solve(obj, b)
            x = obj.solveInPlace(b);
        end

        function x = solveInPlace(obj, b)
            if obj.logical
                error('QDLDL:LogicalFactorization', 'Cannot solve with logical factorization only.');
            end

            x = b;
            if isempty(obj.perm)
                tmp = x;
            else
                obj.workspace.fwork(:) = x(obj.perm);
                tmp = obj.workspace.fwork;
            end

            tmp = QDLDLFactorization.qdldlSolve( ...
                obj.workspace.Ln, ...
                obj.workspace.Lp, ...
                obj.workspace.Li, ...
                obj.workspace.Lx, ...
                obj.workspace.Dinv, ...
                tmp);

            if isempty(obj.perm)
                x = tmp;
            else
                x(obj.perm) = tmp;
            end
        end

        function refactor(obj)
            obj.logical = false;
            obj.workspace = QDLDLFactorization.factorWorkspace(obj.workspace, obj.logical);
            obj.refreshFactors();
        end

        function updateValues(obj, indices, values)
            if isempty(obj.workspace.AtoPAPt)
                obj.workspace.triuA.Ax(indices) = values;
            else
                obj.workspace.triuA.Ax(obj.workspace.AtoPAPt(indices)) = values;
            end
        end

        function scaleValues(obj, indices, scale)
            if isempty(obj.workspace.AtoPAPt)
                obj.workspace.triuA.Ax(indices) = obj.workspace.triuA.Ax(indices) .* scale;
            else
                mapped = obj.workspace.AtoPAPt(indices);
                obj.workspace.triuA.Ax(mapped) = obj.workspace.triuA.Ax(mapped) .* scale;
            end
        end

        function n = positiveInertia(obj)
            n = obj.workspace.positive_inertia;
        end

        function n = regularizedEntries(obj)
            n = obj.workspace.regularize_count;
        end

        function refreshFactors(obj)
            obj.L = QDLDLFactorization.cscStrictLowerToSparse( ...
                obj.workspace.Ln, ...
                obj.workspace.Lp, ...
                obj.workspace.Li, ...
                obj.workspace.Lx);
            obj.Dinv = spdiags(obj.workspace.Dinv, 0, obj.workspace.Ln, obj.workspace.Ln);
        end
    end

    methods (Static)
        function obj = fromMatrix(A, options)
            if ~issparse(A)
                A = sparse(A);
            end

            [m, n] = size(A);
            if m ~= n
                error('QDLDL:DimensionMismatch', 'Input matrix must be square.');
            end

            if ~istriu(A)
                A = triu(A);
            else
                A = sparse(A);
            end

            if ~isfield(options, 'perm') || (ischar(options.perm) && strcmp(options.perm, 'auto'))
                perm = symamd(A);
                perm = perm(:);
            elseif isempty(options.perm)
                perm = [];
            else
                perm = options.perm(:);
                if numel(perm) ~= n || ~isequal(sort(double(perm(:))), (1:n)')
                    error('QDLDL:InvalidPermutation', 'perm must be a permutation of 1:n.');
                end
            end

            iperm = [];
            if ~isempty(perm)
                iperm = QDLDLFactorization.inversePerm(perm);
            end

            cscA = QDLDLFactorization.sparseToCSCUpper(A);
            if ~isempty(perm)
                [cscA, AtoPAPt] = QDLDLFactorization.permuteSymmetricCSC(cscA, iperm);
            else
                AtoPAPt = [];
            end

            mySigns = [];
            if isfield(options, 'Dsigns') && ~isempty(options.Dsigns)
                signs = options.Dsigns(:);
                if numel(signs) ~= n
                    error('QDLDL:DimensionMismatch', 'Dsigns must have one entry per row/column.');
                end
                if isempty(perm)
                    mySigns = signs;
                else
                    mySigns = signs(perm);
                end
            end

            workspace = QDLDLFactorization.initializeWorkspace( ...
                cscA, ...
                AtoPAPt, ...
                mySigns, ...
                options.regularize_eps, ...
                options.regularize_delta);

            workspace = QDLDLFactorization.factorWorkspace(workspace, options.logical);

            L = QDLDLFactorization.cscStrictLowerToSparse(workspace.Ln, workspace.Lp, workspace.Li, workspace.Lx);
            Dinv = spdiags(workspace.Dinv, 0, workspace.Ln, workspace.Ln);

            obj = QDLDLFactorization(perm, iperm, L, Dinv, workspace, options.logical);
        end
    end

    methods (Static, Access = private)
        function ip = inversePerm(p)
            ip = zeros(size(p));
            ip(double(p)) = 1:numel(p);
        end

        function csc = sparseToCSCUpper(S)
            n = size(S, 1);
            [rows, cols, vals] = find(S);

            counts = accumarray(cols, 1, [n, 1]);
            Ap = [1; cumsum(counts) + 1];

            csc = struct();
            csc.n = n;
            csc.Ap = Ap;
            csc.Ai = rows;
            csc.Ax = vals;
        end

        function [P, AtoPAPt] = permuteSymmetricCSC(A, iperm)
            n = A.n;
            Ap = A.Ap;
            Ai = A.Ai;
            Ax = A.Ax;

            numEntries = zeros(n, 1);
            for colA = 1:n
                colP = iperm(colA);
                for idx = Ap(colA):(Ap(colA + 1) - 1)
                    rowA = Ai(idx);
                    if rowA <= colA
                        rowP = iperm(rowA);
                        colIdx = max(rowP, colP);
                        numEntries(colIdx) = numEntries(colIdx) + 1;
                    end
                end
            end

            Pc = zeros(n + 1, 1);
            Pc(1) = 1;
            for k = 1:n
                Pc(k + 1) = Pc(k) + numEntries(k);
                numEntries(k) = Pc(k);
            end

            rowStarts = numEntries;
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
                        colIdx = max(colP, rowP);
                        rowPIdx = rowStarts(colIdx);

                        Pr(rowPIdx) = min(colP, rowP);
                        Pv(rowPIdx) = Ax(rowAIdx);
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
            n = cscA.n;

            etree = zeros(n, 1);
            Lnz = zeros(n, 1);
            iwork = zeros(3 * n, 1);
            bwork = false(n, 1);
            fwork = zeros(n, 1);

            [sumLnz, Lnz, etree] = QDLDLFactorization.qdldlEtree(n, cscA.Ap, cscA.Ai, Lnz, etree);
            if sumLnz < 0
                error('QDLDL:InvalidInput', 'Input matrix is not upper triangular or has an empty column.');
            end

            workspace = struct();
            workspace.etree = etree;
            workspace.Lnz = Lnz;
            workspace.iwork = iwork;
            workspace.bwork = bwork;
            workspace.fwork = fwork;

            workspace.Ln = n;
            workspace.Lp = zeros(n + 1, 1);
            workspace.Li = zeros(sumLnz, 1);
            workspace.Lx = zeros(sumLnz, 1, class(cscA.Ax));

            workspace.D = zeros(n, 1, class(cscA.Ax));
            workspace.Dinv = zeros(n, 1, class(cscA.Ax));
            workspace.positive_inertia = -1;

            workspace.triuA = cscA;
            workspace.AtoPAPt = AtoPAPt;

            workspace.Dsigns = Dsigns;
            workspace.regularize_eps = regularizeEps;
            workspace.regularize_delta = regularizeDelta;
            workspace.regularize_count = 0;
        end

        function workspace = factorWorkspace(workspace, logicalFactor)
            if logicalFactor
                workspace.Lx(:) = 1;
                workspace.D(:) = 1;
                workspace.Dinv(:) = 1;
            end

            [posDCount, workspace] = QDLDLFactorization.qdldlFactor(workspace, logicalFactor);
            if posDCount < 0
                error('QDLDL:NonQuasiDefinite', 'Zero entry in D (matrix is not quasidefinite).');
            end
            workspace.positive_inertia = posDCount;
        end

        function [sumLnz, Lnz, etree] = qdldlEtree(n, Ap, Ai, Lnz, etree)
            work = zeros(n, 1);

            for i = 1:n
                work(i) = 0;
                Lnz(i) = 0;
                etree(i) = -1;

                if Ap(i) == Ap(i + 1)
                    sumLnz = -1;
                    return;
                end
            end

            for j = 1:n
                work(j) = j;
                for p = Ap(j):(Ap(j + 1) - 1)
                    i = Ai(p);
                    if i > j
                        sumLnz = -1;
                        return;
                    end

                    while work(i) ~= j
                        if etree(i) == -1
                            etree(i) = j;
                        end
                        Lnz(i) = Lnz(i) + 1;
                        work(i) = j;
                        i = etree(i);
                    end
                end
            end

            sumLnz = sum(Lnz);
        end

        function [positiveValuesInD, workspace] = qdldlFactor(workspace, logicalFactor)
            n = workspace.triuA.n;
            Ap = workspace.triuA.Ap;
            Ai = workspace.triuA.Ai;
            Ax = workspace.triuA.Ax;

            Lp = workspace.Lp;
            Li = workspace.Li;
            Lx = workspace.Lx;
            D = workspace.D;
            Dinv = workspace.Dinv;
            Lnz = workspace.Lnz;
            etree = workspace.etree;
            yMarkers = workspace.bwork;
            yVals = workspace.fwork;

            yIdx = zeros(n, 1);
            elimBuffer = zeros(n, 1);
            LNextSpaceInCol = zeros(n, 1);

            positiveValuesInD = 0;
            workspace.regularize_count = 0;

            Lp(1) = 1;
            for i = 1:n
                Lp(i + 1) = Lp(i) + Lnz(i);
                yMarkers(i) = false;
                yVals(i) = 0.0;
                D(i) = 0.0;
                LNextSpaceInCol(i) = Lp(i);
            end

            if logicalFactor
                D(:) = 1;
                Dinv(:) = 1;
            else
                D(1) = Ax(1);

                if ~isempty(workspace.Dsigns) && workspace.Dsigns(1) * D(1) < workspace.regularize_eps
                    D(1) = workspace.regularize_delta * workspace.Dsigns(1);
                    workspace.regularize_count = workspace.regularize_count + 1;
                end

                if D(1) == 0.0
                    positiveValuesInD = -1;
                    return;
                end
                if D(1) > 0.0
                    positiveValuesInD = positiveValuesInD + 1;
                end
                Dinv(1) = 1 / D(1);
            end

            for k = 2:n
                nnzY = 0;

                for i = Ap(k):(Ap(k + 1) - 1)
                    bidx = Ai(i);

                    if bidx == k
                        D(k) = Ax(i);
                        continue;
                    end

                    yVals(bidx) = Ax(i);
                    nextIdx = bidx;

                    if ~yMarkers(nextIdx)
                        yMarkers(nextIdx) = true;
                        elimBuffer(1) = nextIdx;
                        nnzE = 1;

                        nextIdx = etree(bidx);
                        while nextIdx ~= -1 && nextIdx < k
                            if yMarkers(nextIdx)
                                break;
                            end

                            yMarkers(nextIdx) = true;
                            nnzE = nnzE + 1;
                            elimBuffer(nnzE) = nextIdx;
                            nextIdx = etree(nextIdx);
                        end

                        while nnzE ~= 0
                            nnzY = nnzY + 1;
                            yIdx(nnzY) = elimBuffer(nnzE);
                            nnzE = nnzE - 1;
                        end
                    end
                end

                for i = nnzY:-1:1
                    cidx = yIdx(i);
                    tmpIdx = LNextSpaceInCol(cidx);

                    if ~logicalFactor
                        yVals_cidx = yVals(cidx);
                        for j = Lp(cidx):(tmpIdx - 1)
                            yVals(Li(j)) = yVals(Li(j)) - Lx(j) * yVals_cidx;
                        end

                        Lx(tmpIdx) = yVals_cidx * Dinv(cidx);
                        D(k) = D(k) - yVals_cidx * Lx(tmpIdx);
                    end

                    Li(tmpIdx) = k;
                    LNextSpaceInCol(cidx) = LNextSpaceInCol(cidx) + 1;

                    yVals(cidx) = 0.0;
                    yMarkers(cidx) = false;
                end

                if ~logicalFactor
                    if ~isempty(workspace.Dsigns) && workspace.Dsigns(k) * D(k) < workspace.regularize_eps
                        D(k) = workspace.regularize_delta * workspace.Dsigns(k);
                        workspace.regularize_count = workspace.regularize_count + 1;
                    end

                    if D(k) == 0.0
                        positiveValuesInD = -1;
                        return;
                    end
                    if D(k) > 0.0
                        positiveValuesInD = positiveValuesInD + 1;
                    end

                    Dinv(k) = 1 / D(k);
                end
            end

            workspace.Lp = Lp;
            workspace.Li = Li;
            workspace.Lx = Lx;
            workspace.D = D;
            workspace.Dinv = Dinv;
            workspace.bwork = yMarkers;
            workspace.fwork = yVals;
        end

        function x = qdldlLsolve(n, Lp, Li, Lx, x)
            for i = 1:n
                for j = Lp(i):(Lp(i + 1) - 1)
                    x(Li(j)) = x(Li(j)) - Lx(j) * x(i);
                end
            end
        end

        function x = qdldlLtsolve(n, Lp, Li, Lx, x)
            for i = n:-1:1
                for j = Lp(i):(Lp(i + 1) - 1)
                    x(i) = x(i) - Lx(j) * x(Li(j));
                end
            end
        end

        function b = qdldlSolve(n, Lp, Li, Lx, Dinv, b)
            b = QDLDLFactorization.qdldlLsolve(n, Lp, Li, Lx, b);
            b = b .* Dinv;
            b = QDLDLFactorization.qdldlLtsolve(n, Lp, Li, Lx, b);
        end

        function L = cscStrictLowerToSparse(n, Lp, Li, Lx)
            nnzL = Lp(end) - 1;
            if nnzL == 0
                L = sparse(n, n);
                return;
            end

            cols = zeros(nnzL, 1);
            for c = 1:n
                first = Lp(c);
                last = Lp(c + 1) - 1;
                if first <= last
                    cols(first:last) = c;
                end
            end

            L = sparse(Li, cols, Lx, n, n);
        end
    end
end
