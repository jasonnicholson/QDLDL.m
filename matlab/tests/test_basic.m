function tests = test_basic
    tests = functiontests(localfunctions);
end

function testLinearSolves(testCase)
    rng(706);

    m = 20;
    n = 30;

    A = random_psd(n);
    B = sprandn(m, n, 0.2);
    C = -random_psd(m);
    M = [A, B'; B, C];

    b = randn(m + n, 1);
    F = qdldl(M);

    testCase.verifyLessThanOrEqual(norm(F \ b - M \ b, inf), 1e-10);
    testCase.verifyLessThanOrEqual(norm(solve(F, b) - M \ b, inf), 1e-10);

    x = b;
    x = solve_inplace(F, x);
    testCase.verifyLessThanOrEqual(norm(x - M \ b, inf), 1e-10);

    As = sparse(randn(1, 1));
    bs = randn(1, 1);
    testCase.verifyLessThanOrEqual(norm((qdldl(As) \ bs) - (As \ bs), inf), 1e-10);

    Ms = sparse(single(M));
    bs = single(b);
    testCase.verifyLessThanOrEqual(norm(double(qdldl(Ms) \ bs) - (M \ b), inf), 1e-5);
end
