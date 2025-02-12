clear
rng(2024);

% matrix size
n = 100;
dim = n;
% randomly generate eigenvlaues, eigenvectors recover matrix
[A, evals, evecs] = rand_mat_real(n);

% initial state: a certain combination of the eigen-vectors
x0 = evecs * (1:n)';

% lower rank for dmd, and the steps for the simulation snapshots
r = 10; steps = 31;
or_evals = evals(1:r); or_evecs = evecs(:, 1:r);

% simulate snapshots observation
d = zeros([n, steps]);
for i = 1:steps
    d(:, i) = A^(i-1) * x0;
end
X = d(:, 1:end-1); Y = d(:, 2:end);

% best for ref
A = Y * pinv(X);
[evecs, evals] = main_eig(A, r);

% A = Y * V * S * U';
[P, ~] = qr(evecs, 'econ');
R = P' * A * P;
A_approx = P * R * P';
[right_r, lambda] = main_eig(R, r);
[left, ~] = main_eig_left(R, r);
right = P * right_r; left = P * left;

U = right; D = diag(lambda);

err_vec = vecs_distance(right, evecs);
err_val = abs(lambda - evals);

err_ed = norm(Y * pinv(X) * U - U * D, 'fro')

B = pinv(X) * U / D;
err_u = norm(U - Y*B)

Yp = P'*Y; Xp = P'*X;
err_proj = norm( (P'*(Y*pinv(X))*P - Yp*pinv(Xp)) , 'fro')


function [Vs, Ds] = sort_by_val(V, D)
    %%% simply sort the [eigenvectors, eigenvalues] by the module of
    %%% eigenvalues
    eigenvalues = abs(diag(D));
    [~, indices] = sort(eigenvalues, 'descend');
    Vs = V(:, indices);
    D = diag(D);
    Ds = diag(D(indices));
end


function [Vs, Ds] = main_eig(A, r)
    %%% simply sort the [eigenvectors, eigenvalues] by the module of
    %%% eigenvalues
    [V, D] = eig(A);
    [Vs, Ds] = sort_by_val(V, D);
    Vs = Vs(:, 1:r);
    Ds = diag(Ds(1:r, 1:r));
end

function [Vs, Ds] = main_eig_left(A, r)
    %%% simply sort the [eigenvectors, eigenvalues] by the module of
    %%% eigenvalues
    A = A';
    [V, D] = eig(A);
    [Vs, Ds] = sort_by_val(V, D);
    Vs = Vs(:, 1:r);
    Ds = diag(Ds(1:r, 1:r));
end

% random matrix with complex eigenvalues and vectors
% (same as the test case of diff-svd)
function [A, evals, evecs] = rand_mat(n)
    mat = normrnd(0, 1, [n, n]);
    [u, ~] = qr(mat);
    
    mat = normrnd(0, 1, [n, n]);
    [v, ~] = qr(mat);
    
    s = logspace(1, -1, n);
    
    A = u * diag(s) * v';
    [evecs, evals] = main_eig(A, n);
end


% random matrix with real eigenvalues and vectors
function [A, evals, evecs] = rand_mat_real(n)
    evals = logspace(0, -2, n) .* (1 + 0.1*randn(1, n));
    evecs = rand_col(n, n);
    A = evecs * diag(evals) / evecs;
end

function vecs = rand_col(dim, num)
% generate a matrix of column vecters, randomly genereted with directions
% input:
% dim: dimension of the vectors, row of vecs
% num: number of the random vectors. column of vecs
% output:
% vecs: dim*num matrix. each column is a random vector
vecs = unifrnd(0, 2*pi, [dim, num]);
vecs = cos(vecs);
len = vecnorm(vecs, 2, 1);
vecs = vecs ./ len;
end


function projected_term = ortho_projection(basis, ex_basis)
    projected_term = ex_basis - basis * basis' * ex_basis;
end

function err = vecs_distance(vecs1, vecs2)
    prod = vecs1' * vecs2;
    diag_err = diag(abs(prod) - eye(size(prod)));
    err = norm(diag_err, 'fro');
end

function err = vals_distance(vals1, vals2)
    err = norm(vals1 - vals2, 'fro');
end
