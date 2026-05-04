function tests = bench_solve
    tests = functiontests(localfunctions);
end
function testBenchSolve(testCase)
    n = 2000;
    A = speye(n) + sprandsym(n, 5/n, 0.1, 1);
    A = A + n * speye(n); 
    b = randn(n, 1);
    F = qdldl(A);
    tic;
    for i=1:1000
        x = solve(F, b);
    end
    t = toc;
    fprintf('1000 solves took %f s\n', t);
end
