function T = tuckerNonneg(X, R, maxitr, tol)
% Solves nonnegative Tucker model by nonnegative alternating least squares
%
% Inputs:
%   X: tensor to be decomposed (tensor, sptensor)
%   R: core tensor dimension
%   maxitr: maximum number of iterations (default 50)
%   tol: convergence criterion when changes in fits < tol between succussive itertions
%
% Outputs:
%   T: ttensor
%
% This function uses the tensor toolbox by Kolda etc and the nonnegative
% least squares solver (NNLS) pfls from N-way toolbox by Bro

%% Set default parameter values
if ~exist('maxitr', 'var'), maxitr = 50; end
if ~exist('tol', 'var'), tol = 1E-4; end

%% Extract number of dimensions and norm of X.
N = ndims(X);
normX = norm(X);

%% Initialize factor and core matrices
U = cell(N,1);
for n = 1:N
    U{n} = rand(size(X,n),R(n));
end
V = 1;
for i = N:-1:1
    V = kron(V,U{i});
end
core = tensor(reshape(fastnnls(V'*V, V'*X(:)), R));

fit = 0;

%% Main ALS loop
for iter = 1:maxitr

    fitold = fit;

    % Iterate over all N modes of the tensor
    for n = 1:N
        % compute consecutive Kronecker products of the factor matrices in reverse order
        Z = 1;
        for i = N:-1:1
            if i == n, continue; end
            Z = kron(Z, U{i});
        end
        
        Z = Z * double(tenmat(core, n))';
        Xn = double(tenmat(X, n));         
        ZtZ = Z'*Z;
        ZtX = Z'*Xn';

        % solve Xn = UnZ' with NNLS by Bro
        constraint = 2; % nonnegativity
        Unew = pfls(ZtZ, ZtX, size(Xn, 1), constraint, U{n}, 0, []);

        % normalize columns of factor matrix
        lambda = sqrt(sum(Unew.^2,1))';
        Unew = Unew * spdiags(1./lambda, 0, R(n), R(n));
        U{n} = full(Unew);
        core = ttm(core, diag(lambda), n);
    end

    % Update core tensor
    V = 1;
    for i = N:-1:1
        V = kron(V,U{i});
    end
    core = tensor(reshape(fastnnls(V'*V, V'*X(:)), R));

    % Compute fit
    T = ttensor(core, U);
    normresidual = sqrt( normX^2 + norm(T)^2 - 2 * innerprod(X, T));
    fit = 1 - (normresidual / normX); %fraction explained by model
    fitchange = abs(fitold - fit);

    fprintf('Iteration %2d: fit = %e\tfit change = %7.1e\n', iter, fit, fitchange);
    
    % Check for convergence
    if (iter > 1) && (fitchange < tol)
        break;
    end

end
end
