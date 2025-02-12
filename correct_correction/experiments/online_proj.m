clear
rng(2024);

% matrix size
n = 100;
dim = n;
% randomly generate eigenvlaues, eigenvectors recover matrix
[A, evals, evecs] = rand_mat_real(n);

% initial state: a certain combination of the eigen-vectors
x0 = evecs * (1:n)';

% steps for the simulation snapshots
steps = 51;
% simulate snapshots observation
D = zeros([n, steps]);
D(:, 1) = x0;
for i = 2:steps
    D(:, i) = A * D(:, i-1);
end
data = D;
X = D(:, 1:end-1); Y = D(:, 2:end);

% best for ref
r = 5; 
A = Y * pinv(X);
[V, D] = main_eig(A, r);


% init for online subspace iteration (not practical, just for reference)
init = r+1;
Dt = data(:, 1:init);
X = Dt(:, 1:end-1);
Y = Dt(:, 2:end);
[V_on, D_on, Q_on] = simple_dmd(X, Y, r);

% current sota
[V_sota, D_sota] = direct_solve(X, Y, r);

% dmd
[V_dmd, D_dmd, Q_dmd] = simple_dmd(X, Y, r);

err_vec_on = vecs_distance(V_on, V);
err_val_on = vals_distance(D_on, D);

err_vec_sota = vecs_distance(V_sota, V);
err_val_sota = vals_distance(D_sota, D);

err_vec_dmd = vecs_distance(V_dmd, V);
err_val_dmd = vals_distance(D_dmd, D);

rate = 1;
U = V_on; u = Q_on' * U; P = Q_on;
% online iteration
for i = 1: steps-init

    % data preparation
    Dt = data(:, 1:init+i);
    Xt = Dt(:, 1:end-1);
    Yt = Dt(:, 2:end);
    A = Yt*pinv(Xt);

    % online preparation
    D_on = diag(U' * A * U);  % diag(U'*A*U / U'*U)
    [U, D_on] = sort_by_abs(U, D_on);
    
    for j = 1:r
        Ur = U(:, j);
        term1 = A * Ur / D_on(j);
        term3 = U * pinv( diag(D_on(j) - D_on) ) * pinv(U) * A * Ur;

        dU(:, j) = term1 + term3;
    end


    U = U + rate * dU;
    P = P + rate * dU / u;
    [P, ~] = qr(P, 'econ');
    u = P' * U;
    length = vecnorm(u);
    u = u ./ length;
    U = P * u;

    % online variables
    Q_on = P; V_on = U;

    err_vec_on = [err_vec_on, vecs_distance(V_on, V)];
    err_val_on = [err_val_on, vals_distance(D_on, D)];

    % other 2 ref
    % current sota
    [V_sota, D_sota] = direct_solve(Xt, Yt, r);
    
    % dmd
    [V_dmd, D_dmd, Q_dmd] = simple_dmd(Xt, Yt, r);
    
    err_vec_sota = [err_vec_sota, vecs_distance(V_sota, V)];
    err_val_sota = [err_val_sota, vals_distance(D_sota, D)];
    
    err_vec_dmd = [err_vec_dmd, vecs_distance(V_dmd, V)];
    err_val_dmd = [err_val_dmd, vals_distance(D_dmd, D)];
end


figure(4)
subplot(121)
hold on
plot(err_vec_on)
plot(err_vec_sota)
plot(err_vec_dmd)
% yscale log
legend('online', 'sota', 'dmd')

subplot(122)
hold on
plot(err_val_on)
plot(err_val_sota)
plot(err_val_dmd)
% yscale log
legend('online', 'sota', 'dmd')






function [V_sort, D_sort] = sort_by_abs(V, D)
    [~, idx] = sort(abs(D), 'descend');
    D_sort = D(idx);
    V_sort = V(:, idx);
end


function [Vs, Ds] = main_eig(A, r)
    %%% simply sort the [eigenvectors, eigenvalues] by the module of
    %%% eigenvalues
    [V, D] = eig(A);
    D = diag(D);
    [Vs, Ds] = sort_by_abs(V, D);
    Vs = Vs(:, 1:r);
    Ds = Ds(1:r);
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

function [V, D, Q] = simple_dmd(X, Y, r)
    [u, s, v] = svds(X, r);
    R = u' * Y * v / s;
    A_dmd = u * R * u';
    [V, D] = main_eig(R, r);
    Q = u;
    V = Q * V;
end

function [V, D] = direct_solve(X, Y, r)
    A = Y * pinv(X);
    [V, D] = main_eig(A, r);
end
    
