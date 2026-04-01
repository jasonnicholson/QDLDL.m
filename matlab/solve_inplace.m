function x = solve_inplace(F, b)
%SOLVE_INPLACE Solve Ax=b and return the overwritten vector.
    x = F.solveInPlace(b);
end
