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
