clear all
close all
rng(2024);

% matrix size
n = 400;
dim = n;
% randomly generate eigenvlaues, eigenvectors recover matrix
[A, evals, evecs] = rand_mat_real(n);
A_org = A;

% initial state: a certain combination of the eigen-vectors
x0 = evecs * (1:n)';

% lower rank for dmd, and the steps for the simulation snapshots
r = 10; steps = 201;

% simulate snapshots observation
d = zeros([n, steps]);
for i = 1:steps
    d(:, i) = A^(i-1) * x0;
end
X = d(:, 1:end-1); Y = d(:, 2:end);

% best for ref
% A = Y * pinv(X);
[evecs, evals] = main_eig(A, r);

iter = 200;
rate = 1;
init = 21;

% dmd result
X = d(:, 1:init-1);
Y = d(:, 2:init);
[U, S, V] = svds(X, r);
R = U' * Y * V * S;
[right_r, lambda] = main_eig(R, r);
right = U * right_r;
P = right; % P is in span(X)

err_vec = vecs_distance(evecs, right);
err_val = vals_distance(lambda, evals);

for i = init: steps
    X = d(:, 1:i-1);
    Y = d(:, 2:i);
    A = Y * pinv(X);
    R = pinv(P) * A * P;
    [right_r, lambda] = main_eig(R, r);
    right = P * right_r;

    for n = 1:r
        u = right(:, n);

        term1 = A * u / lambda(n) - u;
        term3 = right * pinv( diag(lambda(n) - lambda) )/right_r* pinv(P) * A * u;

        du(:, n) = term1 + term3;
    end

    P = right + du;
    norm_P = vecnorm(P);
    P = P ./ norm_P;
    [P, ~] = qr(P + du / right_r, 'econ');
    

    % record error
    err_vec = [err_vec, vecs_distance(evecs, right)];
    err_val = [err_val, vals_distance(lambda, evals)];

    i
end


figure('Position', [100, 100, 800, 600])
subplot(211)
plot(err_vec)
yscale log
subtitle('error of eigenvectors')

subplot(212)
box on
hold on
plot(err_val)
yscale log
legend()
subtitle('error of eigenvalues')




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


function [A, evals, evecs] = case1(n)
    cut = 5;
    evals = [logspace(0.05, -0.05, cut), logspace(-1, -2, n-cut)];
    evecs = rand_col(n, n);
    A = evecs * diag(evals) / evecs;
end




