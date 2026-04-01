function tests = test_scale_values
    tests = functiontests(localfunctions);
end

function testScaleValues(testCase)
    rng(1312);

    nz = 100;
    nc = 70;

    H = sprand(nz, nz, 0.05);
    H = H' * H + speye(nz);
    b = randn(nz + nc, 1);

    A1 = sprand(nc, nz, 0.8);
    K1 = [H, A1'; A1, -1e-3 * speye(nc)];

    triuK1 = triu(K1);
    F = qdldl(triuK1);

    testCase.verifyLessThan(norm(K1 \ b - F \ b, inf), 1e-12);

    % Scale all nonzeros by 2 and verify the factorization of 2*K1
    [~, ~, vals] = find(triuK1);
    scale_values(F, (1:numel(vals))', 2.0);
    refactor_qdldl(F);

    testCase.verifyLessThan(norm((2 * K1) \ b - F \ b, inf), 1e-12);

    % Scale back by 0.5 to recover K1 and verify again
    scale_values(F, (1:numel(vals))', 0.5);
    refactor_qdldl(F);

    testCase.verifyLessThan(norm(K1 \ b - F \ b, inf), 1e-12);
end
