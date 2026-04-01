function tests = test_update_values
    tests = functiontests(localfunctions);
end

function testUpdateValues(testCase)
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

    A2 = A1;
    [r, c, ~] = find(A2);
    A2 = sparse(r, c, randn(numel(r), 1), size(A2, 1), size(A2, 2));

    K2 = [H, A2'; A2, -1e-7 * speye(nc)];
    triuK2 = triu(K2);

    [~, ~, vals] = find(triuK2);
    update_values(F, (1:numel(vals))', vals);
    refactor_qdldl(F);

    testCase.verifyLessThan(norm(K2 \ b - F \ b, inf), 1e-12);
end
