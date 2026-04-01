function tests = test_regularization
    tests = functiontests(localfunctions);
end

function testRegularization(testCase)
    a = 1e-20;
    A = sparse([a, 0; 0, a]);

    signs = [1; 1];
    f = qdldl(A, 'Dsigns', signs);
    testCase.verifyEqual(f.workspace.D, [1e-7; 1e-7], 'AbsTol', 0);
    testCase.verifyEqual(regularized_entries(f), 2);

    signs = [1; -1];
    f = qdldl(A, 'Dsigns', signs);
    testCase.verifyEqual(f.workspace.D, [1e-7; -1e-7], 'AbsTol', 0);
    testCase.verifyEqual(regularized_entries(f), 2);

    delta = 1e-8;
    f = qdldl(A, 'Dsigns', signs, 'regularize_delta', delta);
    testCase.verifyEqual(f.workspace.D, [delta; -delta], 'AbsTol', 0);
    testCase.verifyEqual(regularized_entries(f), 2);

    epsValue = 1e-30;
    f = qdldl(A, 'Dsigns', signs, 'regularize_eps', epsValue, 'regularize_delta', delta);
    testCase.verifyEqual(f.workspace.D, [a; -delta], 'RelTol', 1e-12);
    testCase.verifyEqual(regularized_entries(f), 1);

    % MATLAB sparse assignment removes explicit zeros, so exercise the same
    % structural-zero case by zeroing an existing stored value in-place.
    f = qdldl(A, 'Dsigns', signs);
    update_values(f, 1, 0);
    refactor_qdldl(f);
    testCase.verifyEqual(f.workspace.D, [1e-7; -1e-7], 'AbsTol', 0);

    g = qdldl(A);
    update_values(g, 1, 0);
    threw = false;
    try
        refactor_qdldl(g);
    catch
        threw = true;
    end
    testCase.verifyTrue(threw);

    A = sparse([0, 0; 0, a]);
    threw = false;
    try
        qdldl(A);
    catch
        threw = true;
    end
    testCase.verifyTrue(threw);

    threw = false;
    try
        qdldl(A, 'Dsigns', signs);
    catch
        threw = true;
    end
    testCase.verifyTrue(threw);

    A = speye(3) * a;
    signs = [1; 1; -1];
    perm = [2; 3; 1];
    f = qdldl(A, 'perm', perm, 'Dsigns', signs);
    testCase.verifyEqual(f.workspace.D, [1e-7; -1e-7; 1e-7], 'AbsTol', 0);
    testCase.verifyEqual(regularized_entries(f), 3);

    rng(2401);
    m = 20;
    n = 30;
    A = random_psd(n) + spdiags(rand(n, 1), 0, n, n);
    B = sprandn(m, n, 0.2);
    C = -spdiags(rand(m, 1), 0, m, m);
    M = [A, B'; B, C];
    b = randn(m + n, 1);
    s = [ones(n, 1); -ones(m, 1)];
    p = randperm(m + n)';
    F = qdldl(M, 'perm', p, 'Dsigns', s);
    testCase.verifyLessThanOrEqual(norm(F \ b - M \ b, inf), 1e-10);
end
