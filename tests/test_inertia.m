function tests = test_inertia
    tests = functiontests(localfunctions);
end

function testInertiaChecks(testCase)
    rng(2701);

    m = 20;
    n = 30;

    A = random_psd(n) + spdiags(rand(n, 1), 0, n, n);
    B = sprandn(m, n, 0.2);
    C = -spdiags(rand(m, 1), 0, m, m);
    M = [A, B'; B, C];

    F = qdldl(M);
    testCase.verifyEqual(positive_inertia(F), n);
end
