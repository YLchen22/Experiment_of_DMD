clear all
rng(2024);

% matrix size
n = 200;
dim = n;
% randomly generate eigenvlaues, eigenvectors recover matrix
[A, evals, evecs] = rand_mat_real(n);

% initial state: a certain combination of the eigen-vectors
x0 = evecs * ones(n, 1);

% steps for the simulation snapshots
steps = 200;

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
U = V_on; UR = Q_on' * U; P = Q_on; 
B = pinv(Dt) * P;
% P == Dt*B, U == P*UR == Dt*B*UR

% online iteration
for i = 1: steps-init

    % data preparation
    Dt = data(:, 1:init+i);
    D_old = Dt(:, 1:end-1);
    D_new = Dt(:, 2:end);
    % A * D_old == D_new

    % preparation
    A_r = B' * D_old' * D_new * B;
    [UR, lambda] = main_eig(A_r, r); 
    U = P*UR;
    BU = B * UR;

    Uinv = UR \ P';
    new_BU = [BU; zeros(1, r)];
    for j = 1:r
        Bj = BU(:, j);

        temp1 = Bj / lambda(j);
        term1 = [0; temp1];

        temp2 = BU * pinv(diag(lambda(j) - lambda)) * Uinv * D_new * Bj;
        term2 = [temp2; 0];

        new_BU(:, j) = term1 + term2;   % dU == Dt*(new_BU-BU)
    end

    BU = [BU; zeros(1, r)]; % U == D_old*BU ---> U == Dt*BU
    delta_BU = new_BU - BU; % dU == Dt*delta_BU
    B_ex = [B; zeros(1, r)]; % P == D_old*B ---> P == Dt*B_ex

    %%% route 1, closed-form
    new_B = delta_BU / UR + B_ex; % new_P == dP + P == Dt*new_B
    new_P = Dt * new_B;
    err1(i) = norm(new_P - Dt * pinv(Dt) * new_P, 'fro');

    %%% route 2, directly adding up
    % dU = Dt * delta_BU;
    % dP = dU / UR;
    % new_P = P + dP;

    % [P, R] = qr(new_P, 'econ'); Rinv = inv(R);
    % [P, s, v] = svd(new_P, 'econ'); Rinv = v * pinv(s);
    [P, R, T] = qr(new_P, 'econ'); Rinv = T / R;

    B1 = new_B * Rinv;     % closed-form parameter

    B2 = pinv(Dt) * P;  % directly calculated (not practical) parameter
    B = B1;
    err(i) = norm(Dt*(B1 - B2), 'fro');    % why???????????????????????
    err_ana(i) = norm(P - Dt*B1, 'fro');
    err_ref(i) = norm(P - Dt*B2, 'fro');   % why???????????????????????
    err_qr(i) = norm(P - new_P * Rinv, 'fro');

    V_on = U; D_on = lambda;

    err_vec_on = [err_vec_on, vecs_distance(V_on, V)];
    err_val_on = [err_val_on, vals_distance(D_on, D)];

    % other 2 ref
    % current sota
    Xt = D_old; Yt = D_new;
    [V_sota, D_sota] = direct_solve(Xt, Yt, r);
    
    % dmd
    [V_dmd, D_dmd, Q_dmd] = simple_dmd(Xt, Yt, r);
    
    err_vec_sota = [err_vec_sota, vecs_distance(V_sota, V)];
    err_val_sota = [err_val_sota, vals_distance(D_sota, D)];
    
    err_vec_dmd = [err_vec_dmd, vecs_distance(V_dmd, V)];
    err_val_dmd = [err_val_dmd, vals_distance(D_dmd, D)];
end

figure()
hold on
plot(err_ana)
plot(err_ref)
plot(err_qr)
yscale log
legend('Analytic', 'Reference', 'location', 'best')
title('Error caused by parameter matrix - reference form')

figure()
hold on
plot(err1)
yscale log

figure('Position', [100, 100, 1200, 600])
subplot(121)
box on
hold on
plot(err_vec_on)
plot(err_vec_sota)
plot(err_vec_dmd)
% yscale log
legend('Online algorithm', 'SOTA', 'DMD', 'location', 'best')
subtitle('Error of eigenvector')

subplot(122)
box on
hold on
plot(err_val_on)
plot(err_val_sota)
plot(err_val_dmd)
% yscale log
legend('Online algorithm', 'SOTA', 'DMD', 'location', 'best')
subtitle('Error of eigenvalue')





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

function [V, D, Q] = simple_dmd(X, Y, r)
    [u, s, v] = svds(X, r);
    R = u' * Y * v * pinv(s);
    A_dmd = u * R * u';
    [V, D] = main_eig(R, r);
    Q = u;
    V = Q * V;
end

function [V, D] = direct_solve(X, Y, r)
    A = Y * pinv(X);
    [V, D] = main_eig(A, r);
end
    

function [V_sort, D_sort] = sort_by_abs(V, D)
    [~, idx] = sort(abs(D), 'descend');
    D_sort = D(idx);
    V_sort = V(:, idx);
end
