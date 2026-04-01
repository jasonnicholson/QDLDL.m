function tests = test_non_quasidef
    tests = functiontests(localfunctions);
end

function testNonQuasidefinite(testCase)
    rng(1312);

    m = 20;
    A = random_psd(m);
    A(:, 10) = 0;
    A(10, :) = 0;

    threw = false;
    try
        qdldl(A);
    catch
        threw = true;
    end
    testCase.verifyTrue(threw);
end
