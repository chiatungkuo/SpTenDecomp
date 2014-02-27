function T = cpNonnegSp(X, numf, spmode, tol, maxitr, spiter)
% Solves nonnegative Parafac model by nonnegative alternating least squares
% with sparse projection. 
%
% Inputs:
%   X: tensor to be decomposed (tensor, sptensor)
%   numf: number of factors
%   spmode: vector of size the same as number of dimensions of X,
%       indicating desired sparsity [0, 1) on each mode; 0 if sparsity not desired
%   tol: convergence criterion when changes in fits is smaller than tol between 
%       succussive itertions (default 1e-4)
%   maxitr: maximum number of iterations (default 50)
%   spiter: maximum number of iterations with sparse projection (<=
%       maxitr); an estimate is used if not provided by the user.
%
% Outputs:
%   T: ttensor
%
% This function uses the tensor toolbox by Kolda etc ,the nonnegative
% least squares solver (NNLS) pfls from N-way toolbox by Bro, and the
% sparse projection implementation from Hoyer's nmfpack
%

%% Set default parameter values
if ~exist('spmode', 'var'), spmode = zeros(ndims(X), 1); end
if ~exist('tol', 'var'), tol = 1E-4; end
if ~exist('maxitr', 'var'), maxitr = 50; end
if ~exist('spiter', 'var'), spiter = round(0.6*maxitr); end

%% check for potential errors
if ~isempty(find(spmode > 1-tol, 1))
    fprintf('Setting desired sparity (very close) to 1 may cause instability in the algorithm');
end

%% Extract number of dimensions and norm of X.
N = ndims(X);
normX = norm(X);

%% Initialize factor and core matrices
U = cell(N, 1);
for n = 1:N
    U{n} = rand(size(X, n), numf);
end

fit = 0;
constraint = 2; % nonnegativity

%% Main ALS loop
for iter = 1:maxitr

    fitold = fit;
    
    for n = 1:N
        
        %% solve NNLS subproblems row by row
        Unew = U{n};
        Uexceptn = [U(1:n-1); U(n+1:end)];
        Z = khatrirao(Uexceptn, 'r');
        Xn = double(tenmat(X, n));
        ZtZ = Z'*Z;
        
        % do not flip over zeros if sparity is desired;
        % only work on currently remaining positive entries
        if spmode(n) > 0           
            for i = 1:size(Xn,1)
                old_a = Unew(i,:);
                gtz_ind = old_a > 0;
                if isempty(gtz_ind)
                    fprintf('Factor matrix %d row %d has no positive entries.\n', n, i)
                    continue;
                end

                xn = Xn(i,:);
                ztx = Z(:,gtz_ind)'*xn';
                ztz = ZtZ(gtz_ind,gtz_ind);
                a = pfls(ztz,ztx,1,2,old_a(gtz_ind),0,[]);
                Unew(i,gtz_ind) = a;
            end
        else
            ZtX = Z'*Xn';
            Unew = pfls(ZtZ, ZtX, size(Xn, 1), constraint, Unew, 0, []);
        end      

        %% normalize columns of factor matrix
        lambda = sqrt(sum(Unew.^2, 1))';
        Unew = Unew * spdiags(1./lambda, 0, numf, numf);
        
        %% project each factor to the nearest point with the increasing sparsity
        if spmode(n) > 0 && iter <= spiter
            len = size(Unew, 1);
            sqrtlen = sqrt(len);
            for i = 1:numf
                f = Unew(:,i);
                fsp = sparseness(f, len);
                targetsp = fsp + iter * (spmode(n) - fsp) / spiter;
                if fsp >= targetsp, continue; end
                k1 = sqrtlen - (sqrtlen - 1) * targetsp;
                f = projfunc(f, k1, 1, 1);
                Unew(:, i) = f;
                fprintf('Project mode %d component %d from sparseness %4f to %4f\n', ...
                    n, i, fsp, sparseness(f, len));
            end
        end
        U{n} = Unew;

        %% abort if the resulting factors have nan or Inf
        if hasInfNaN(Unew)
            error('Factor matrix %d has Inf or NaN\n',n);
        end
        
    end

    %% Compute fit and output 
    T = ktensor(lambda, U);
    normresidual = sqrt(normX^2 + norm(T)^2 - 2 * innerprod(X, T));
    fit = 1 - (normresidual / normX); % fraction explained by model
    fitchange = abs(fitold - fit);
    fprintf('Iteration %2d: fit = %e fit change = %7.1e\n', iter, fit, fitchange);
    
    %% Check for convergence
    if (iter > spiter) && (fitchange < tol)
        break;
    end
end
%% arrange the output so that weights are ordered in magnitude
T = arrange(T);

end
