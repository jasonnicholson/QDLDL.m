function tests = bench
    tests = functiontests(localfunctions);
end
function testBench(testCase)
    n = 2000;
    A = speye(n) + sprandsym(n, 5/n, 0.1, 1);
    A = A + n * speye(n); 
    F = qdldl(A);
    tic;
    for i=1:100
        refactor_qdldl(F);
    end
    t = toc;
    fprintf('100 refactors took %f s\n', t);
end
