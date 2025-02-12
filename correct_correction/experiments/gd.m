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

% dmd result
[U, S, V] = svds(X, r);
R = U' * Y * V / S;
A_approx = U * R * U';
[right, lambda] = main_eig(A_approx, r);
[left, ~] = main_eig_left(A_approx, r);

iter = 1000;
rate = 1;

err_vec = vecs_distance(evecs, right);
err_val = abs(lambda - evals);

for i = 1:iter
    du = zeros(dim, r);
    dlambda = zeros(r, 1);
    for n = 1:r
        u = right(:, n);
        v = left(:, n);
        term_1 = pinv(lambda(n) * eye(dim) - A_approx) * (eye(dim) - u*v' / (v'*u));
        term_2 = (A - A_approx) * u;
        du(:, n) = term_1 * term_2;
        dlambda(n) = v' * term_2 / (v'*u);
    end

    % update eigens
    [q, ~] = qr(right, 'econ');
    % du = ortho_projection(q, du);   % !!!: NaN because of this step
    right = right + rate*du;
    right = right ./ vecnorm(right);

    lambda = vecnorm(A * right) ./ vecnorm(right);
    % lambda = lambda + dlambda;
    lambda = lambda';

    % update core matrix
    A_approx = right * diag(lambda) / right;
    [left, ~] = main_eig_left(A_approx, r);

    % record error
    err_vec = [err_vec, vecs_distance(evecs, right)];
    err_val = [err_val, abs(lambda - evals)];
end

figure('Position', [100, 100, 800, 600])
subplot(211)
plot(err_vec)
yscale log
subtitle('error of eigenvectors')

subplot(212)
box on
hold on
for i = 1:r
    plot(err_val(i, :), DisplayName=[num2str(i) '-th \lambda'])
end
yscale log
legend()
subtitle('error of each eigenvalue')




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
function A = rand_mat(n)
    mat = normrnd(0, 1, [n, n]);
    [u, ~] = qr(mat);
    
    mat = normrnd(0, 1, [n, n]);
    [v, ~] = qr(mat);
    
    s = logspace(1, -1, n);
    
    A = u * diag(s) * v';
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
