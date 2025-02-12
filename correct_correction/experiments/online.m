clear
rng(2024);

% matrix size
n = 100;
dim = n;
% randomly generate eigenvlaues, eigenvectors recover matrix
[A, evals, evecs] = rand_mat_real(n);
A_real = A;

% initial state: a certain combination of the eigen-vectors
x0 = evecs * (1:n)';

% lower rank for dmd, and the steps for the simulation snapshots
r = 5; steps = 51;

[or_evecs, or_evals] = main_eig(A, r);
evals = or_evals; evecs = or_evecs;

% simulate snapshots observation
D = zeros([n, steps]);
D(:, 1) = x0;
for i = 2:steps
    D(:, i) = A * D(:, i-1);
end
X = D(:, 1:end-1); Y = D(:, 2:end);

% best for ref
A = Y * pinv(X);
[evecs, evals] = main_eig(A, r);
%%% for reference
A_real = A;

init = 6;
Dt = D(:, 1:init);
X = Dt(:, 1:end-1);
Y = Dt(:, 2:end);
A = Y * pinv(X);

%%% suppose we can begin with an exactly accurate initial ed
% [U, lambda] = main_eig(A, r);

%%% inaccurate initial ed
[U, S, V] = svds(X, r);
R = U' * Y * V / S;
A_approx = U * R * U';
[U, lambda] = main_eig(A_approx, r);

B = pinv(Dt) * U;
U = Dt * B;

rate = 1;
% online iteration
for i = 1:steps - init

    Dt = D(:, 1:init+i);
    D_old = Dt(:, 1:end-1);
    Dp = Dt(:, 2:end);

    lambda = diag(U' * Dp * B);  % diag(U'*A*U / U'*U)

    Uinv = pinv(U);
    new_B = [B; zeros(1, r)];
    for j = 1:r        
        Bj = B(:, j);

        temp1 = Bj / lambda(j);
        term1 = [0; temp1];

        temp2 = B * pinv(diag(lambda(j) - lambda)) * Uinv * Dp * Bj;
        term2 = [temp2; 0];

        new_B(:, j) = term1 + term2;

        %%% for reference
        u = U(:, j);
        ref1 = A_real * D_old * Bj / lambda(j);
        ref2 = U * pinv(diag(lambda(j) - lambda)) * pinv(U) * A_real * D_old * Bj;
        %%% ref1 = Dp * temp1,  ref2 = D_old * temp2
        ref_new_u(:, j) = ref1 + ref2;

        test_u = Dt * new_B(:, j);
        z1 = ref_new_u - test_u;
        norm(z1, 'fro');
        a = 1;
    end

    new_U = Dt * new_B;
    err_u(i) = norm(diag(new_U' * ref_new_u - eye(r)), 'fro');

    unorm = vecnorm(new_U);
    new_B = new_B ./ unorm;

    B = new_B;
    U = Dt * B;
    err_unit(i) = norm(diag(U'*U - eye(r)), 'fro');

    %%% for reference
    ref_new_u = ref_new_u ./ vecnorm(ref_new_u);

    X = Dt(:, 1:end-1);
    Y = Dt(:, 2:end);
    A = Y * pinv(X);
    [ref_evecs, ref_evals] = main_eig(A, r);

    
    % error of the results from the best of current data
    err_vec_data(i) = vecs_distance(evecs, ref_evecs);
    err_val_data(:, i) = abs(evals - ref_evals);

    % error of the results from the updating prediction
    err_vec_pred(i) = vecs_distance(evecs, U);
    err_val_pred(:, i) = abs(evals - lambda);

    % relative error of the prediction to the current best
    err_vec_rel(i) = vecs_distance(ref_evecs, U);
    err_val_rel(:, i) = abs(ref_evals - lambda);
end

dmd_start = 41;
X = D(:, dmd_start: end-1);
Y = D(:, dmd_start+1: end);
[U, S, V] = svds(X, r);
R = U' * Y * V / S;
A_approx = U * R * U';
[dmd_evecs, dmd_evals] = main_eig(R, r);
dmd_evecs = U * dmd_evecs;
err_vec_dmd = vecs_distance(evecs, dmd_evecs);
err_val_dmd = abs(evals - dmd_evals);

figure('Position', [100, 100, 1200, 600])
sgtitle('errors with increasing data. eigenvectors and eigenvalues')
s1 = subplot(211);
plot(err_vec_pred)
yline(err_vec_dmd, 'LineWidth', 2., 'LineStyle', '--')
subtitle('eigenvectors')
yscale log

s2 = subplot(212);
box on
hold on
colors = lines(r);
for i = 1:r
    plot(err_val_pred(i, :), 'Color', colors(i, :)); % plot 使用颜色
    yline(err_val_dmd(i), 'Color', colors(i, :), 'LineWidth', 2., 'LineStyle', '--');    % yline 使用相同颜色
end
subtitle('eigenvalues')
yscale log




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
