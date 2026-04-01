# QDLDL.m — A free LDL factorisation routine for MATLAB

A MATLAB port of [QDLDL.jl](https://github.com/osqp/QDLDL.jl), which is itself a Julia implementation of the C language [QDLDL solver](https://github.com/osqp/qdldl).

QDLDL is a factorisation routine for quasi-definite linear systems `Ax = b`. It includes additional functionality for refactorisations, making it well-suited for iterative solvers such as OSQP.

## Getting Started

Add the `matlab/` folder (and subfolders) to your MATLAB path:

```matlab
addpath('matlab')
```

## Using QDLDL

Given a quasidefinite sparse matrix `A` and right-hand side vector `b`, compute the factorisation `F` with:

```matlab
F = qdldl(A)
```

Solve the linear system for `x` with:

```matlab
x = solve(F, b)
```

To solve and overwrite `b` in-place:

```matlab
x = solve_inplace(F, b)
```

`F` also supports the backslash operator directly:

```matlab
x = F \ b
```

### Options

`qdldl` accepts optional name-value pairs:

| Option | Default | Description |
|---|---|---|
| `'perm'` | `'auto'` | Permutation vector, `[]` to disable reordering, or `'auto'` to use AMD |
| `'amd_dense_scale'` | `1.0` | Dense row/column scale factor for AMD ordering |
| `'logical'` | `false` | Perform symbolic (logical) factorisation only |
| `'Dsigns'` | `[]` | Expected signs of diagonal entries of D |
| `'regularize_eps'` | `1e-12` | Regularisation threshold |
| `'regularize_delta'` | `1e-7` | Regularisation perturbation magnitude |

Example:

```matlab
F = qdldl(A, 'perm', [], 'regularize_delta', 1e-6)
```

### Refactorisation

After updating matrix values, refactor numerically without re-running the symbolic phase:

```matlab
update_values(F, indices, values)   % overwrite selected nonzeros
scale_values(F, indices, scale)     % scale selected nonzeros
refactor_qdldl(F)                   % recompute L and D in-place
```

### Inertia

Get the number of positive diagonal entries in D:

```matlab
p = positive_inertia(F)
```

## Running Tests

```matlab
cd matlab/tests
run_all_tests
```

## Credits

- Original C solver: [osqp/qdldl](https://github.com/osqp/qdldl)
- Julia implementation: [osqp/QDLDL.jl](https://github.com/osqp/QDLDL.jl) by [Paul Goulart](http://users.ox.ac.uk/~engs1373/)

## License

This project is licensed under the Apache License — see the [LICENSE](LICENSE) file for details.
