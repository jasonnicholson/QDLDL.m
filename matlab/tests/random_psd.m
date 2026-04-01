function A = random_psd(n, density, seed)
%RANDOM_PSD Create a sparse, symmetric, diagonally dominant PSD matrix.
    if nargin < 2
        density = 0.2;
    end
    if nargin >= 3
        rng(seed);
    end

    A = sprandn(n, n, density);
    A = A + A';
    A = A + spdiags(sum(abs(A), 2), 0, n, n);
end
