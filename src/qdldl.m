function F = qdldl(A, varargin)
%QDLLDL Factor a sparse quasi-definite matrix using LDL.
%
% F = qdldl(A)
% F = qdldl(A, 'perm', p)
% F = qdldl(A, 'perm', [])         % disable reordering
% F = qdldl(A, 'logical', true)
% F = qdldl(A, 'Dsigns', s, 'regularize_eps', eps, 'regularize_delta', delta)

    options = struct();
    options.amd_dense_scale = 1.0;
    options.perm = 'auto';
    options.logical = false;
    options.Dsigns = [];
    options.regularize_eps = 1e-12;
    options.regularize_delta = 1e-7;

    if mod(numel(varargin), 2) ~= 0
        error('QDLDL:InvalidOptions', 'Name/value options must be provided in pairs.');
    end

    for k = 1:2:numel(varargin)
        key = varargin{k};
        val = varargin{k + 1};
        switch lower(key)
            case 'amd_dense_scale'
                options.amd_dense_scale = val;
            case 'perm'
                options.perm = val;
            case 'logical'
                options.logical = logical(val);
            case 'dsigns'
                options.Dsigns = val;
            case 'regularize_eps'
                options.regularize_eps = val;
            case 'regularize_delta'
                options.regularize_delta = val;
            otherwise
                error('QDLDL:InvalidOptions', 'Unknown option: %s', key);
        end
    end

    F = QDLDLFactorization.fromMatrix(A, options);
end
