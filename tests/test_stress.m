function tests = test_stress
    tests = functiontests(localfunctions);
end

function testLargeScale(testCase)
    % Test a large-scale system (e.g. 2000x2000)
    rng(1);
    n = 2000;
    
    A = speye(n) + sprandsym(n, 5/n, 0.1, 1);
    % Make sure it's quasidefinite (positive definite in this case)
    A = A + n * speye(n); 
    
    b = randn(n, 1);
    
    F = qdldl(A);
    x = F \ b;
    
    res = norm(A*x - b, inf) / norm(b, inf);
    testCase.verifyLessThanOrEqual(res, 1e-8, sprintf('Large scale residual failed: %e', res));
end

function testIllConditioned(testCase)
    % Test highly ill-conditioned quasidefinite matrix
    rng(2);
    n = 100;
    % Create eigenvalues spanning a large but realistic range
    d = logspace(-5, 5, n)';
    Q = qr(randn(n));
    A = Q * diag(d) * Q';
    % sparse it up a bit, though Q*D*Q' is dense.
    A = sparse(A);
    
    b = randn(n, 1);
    F = qdldl(A, 'regularize_eps', 1e-12, 'regularize_delta', 1e-7); 
    x = F \ b;
    res = norm(A*x - b, inf) / norm(b, inf);
    
    % cond(A) is ~1e10. we expect res < 1e-2. Double precision eps is 1e-16.
    testCase.verifyLessThanOrEqual(res, 1e-2, sprintf('Ill-conditioned residual failed: %e', res));
end

function testRepeatedRefactor(testCase)
    rng(3);
    n = 200;
    m = 100;
    A = speye(n);
    B = sprandn(m, n, 0.1);
    C = -speye(m);
    M = [A, B'; B, C];
    
    b = randn(n + m, 1);
    triuM = triu(M);
    F = qdldl(triuM);
    
    % repeated update & refactor
    for i = 1:100
        M_new = M + 1e-3 * speye(n+m);
        triuM_new = triu(M_new);
        [~, ~, vals] = find(triuM_new);
        
        update_values(F, (1:numel(vals))', vals);
        refactor_qdldl(F);
        
        x = solve(F, b);
        res = norm(M_new*x - b, inf) / norm(b, inf);
        testCase.verifyLessThanOrEqual(res, 1e-8);
        M = M_new;
    end
end

function testMultipleRHS(testCase)
    % Verify that solving with a matrix of right-hand sides produces the
    % correct columnwise solution (regression test for permutation-aware
    % solveInPlace handling).
    rng(7);
    n = 60; m = 40;
    A = sprandsym(n, 0.3) + 5*speye(n);
    B = sprandn(m, n, 0.3);
    C = -sprandsym(m, 0.3) - 5*speye(m);
    K = [A, B'; B, C];
    F = qdldl(triu(K));

    Brhs = randn(n + m, 8);
    X = F \ Brhs;

    res = norm(K*X - Brhs, 'fro') / norm(Brhs, 'fro');
    testCase.verifyLessThanOrEqual(res, 1e-10, ...
        sprintf('Multi-RHS solve residual: %e', res));

    % Per-column equivalence with single-RHS solves
    for k = 1:size(Brhs, 2)
        xk = F \ Brhs(:, k);
        testCase.verifyLessThanOrEqual(norm(X(:,k) - xk, inf), 1e-12);
    end
end

function testInertiaIndefinite(testCase)
    % Validate that positive_inertia matches the true number of positive
    % eigenvalues for an indefinite quasi-definite KKT-like matrix.
    rng(11);
    n = 50; m = 30;
    H = sprandsym(n, 0.2) + 5*speye(n);
    G = sprandn(m, n, 0.2);
    Cblk = -sprandsym(m, 0.2) - 5*speye(m);
    K = [H, G'; G, Cblk];

    F = qdldl(triu(K));
    pos_qdldl = positive_inertia(F);
    pos_true  = sum(eig(full(K)) > 0);
    testCase.verifyEqual(pos_qdldl, pos_true, ...
        sprintf('Inertia mismatch: qdldl=%d, true=%d', pos_qdldl, pos_true));
end

function testAMDFillReduction(testCase)
    % Arrowhead matrices have catastrophic fill without a fill-reducing
    % permutation but only O(n) fill with AMD. This is a regression test
    % verifying that the default permutation actually reorders.
    rng(13);
    n = 200;
    offdiag = randn(n - 1, 1);
    A = speye(n) ...
        + sparse(1, 2:n, offdiag, n, n) ...
        + sparse(2:n, 1, offdiag, n, n);
    A = A + n * speye(n);

    F_no   = qdldl(triu(A), 'perm', []);
    F_amd  = qdldl(triu(A));

    % AMD must dramatically reduce fill on this pattern
    testCase.verifyLessThan(nnz(F_amd.L), nnz(F_no.L) / 10, ...
        sprintf('AMD ineffective: nnz no perm=%d, AMD=%d', ...
                nnz(F_no.L), nnz(F_amd.L)));

    % Both should still solve correctly
    b = randn(n, 1);
    res_no  = norm(A * (F_no  \ b) - b) / norm(b);
    res_amd = norm(A * (F_amd \ b) - b) / norm(b);
    testCase.verifyLessThanOrEqual(res_no,  1e-10);
    testCase.verifyLessThanOrEqual(res_amd, 1e-10);
end
