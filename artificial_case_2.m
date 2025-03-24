clear all
close all
rng(2024)

%% hyper parameters
n = 400;    % dimension
steps = 200;    % data size
r = 10;     % lower-rank
init = steps/2;  % initial step for the alg
w = init;

%% generate matrix and data
[A_org, evals, evecs] = rand_mat_real(n);
% [A_org, evals, evecs] = case2(n);

distri = ones(n, 1);
x0 = evecs * distri;
% control the beginning distribution of eigenvectors,
% and simulate snapshots observation

data = zeros([n, steps]);
data(:, 1) = x0;
for i = 2:steps
    data(:, i) = A_org * data(:, i-1);
end

% data = data + 1e-8 * eye(size(data));
% data = data + 1e-6 * randn(size(data));     % Robust!

% reference result
V = evecs(:, 1:r); D = evals(1:r);
[P_real, ~] = qr(V, 'econ');
b = 1e0;
[P_ref, ~] = qr(V + b*eye(size(V)), 'econ');
[P_svd, ~] = svds(data, r);

X = data(:, 1:end-1); Y = data(:, 2:end);
A = Y * pinv(X);
% Ar_real = P_real' * A * P_real;
% Ar_svd = P_svd' * A * P_svd;
% Ar_ref = P_ref' * A * P_ref;
Ar_real = P_real' * Y * pinv(P_real' * X);
Ar_svd = P_svd' * Y * pinv(P_svd' * X);
Ar_ref = P_ref' * Y * pinv(P_ref' * X);

err_real = norm((eye(n) - P_real * P_real') * V);
err_svd = norm((eye(n) - P_svd * P_svd') * V);
err_ref = norm((eye(n) - P_ref * P_ref') * V);
figure()
bar(["Real", "SVD", "Biased ref"], [err_real, err_svd, err_ref])
title('Error of projected evecs for different subspaces')


t = 100000;
distri = randn(n, 1)*3;
pred_x0 = evecs * distri;
pred_x0 = data(:, end);
pred_data = zeros([n, t]);
pred_data(:, 1) = pred_x0;
for i = 2:t
    pred_data(:, i) = A_org * pred_data(:, i-1);
end

pred_real = projected_reconstruction(pred_x0, P_real, Ar_real, t);
pred_svd = projected_reconstruction(pred_x0, P_svd, Ar_svd, t);
pred_ref = projected_reconstruction(pred_x0, P_ref, Ar_ref, t);

err_real = vecnorm(pred_data - pred_real);
err_svd = vecnorm(pred_data - pred_svd);
err_ref = vecnorm(pred_data - pred_ref);
figure()
hold on
plot(err_real)
plot(err_svd)
plot(err_ref)
legend('Real', 'SVD', 'Biased ref')
xlabel('Steps')
ylabel('Vecnorm')
yscale log
title('Error of low-rank prediction for different subspaces')






% figure()
% clear bnorm
% for i = init: steps
%     bnorm(i) = cond(B_on{i});
%     bnormex(i) = cond(B_ex{i});
% end
% hold on
% plot(bnorm)
% plot(bnormex)
% title('Conditional number of the parameter matrix B')
% legend('online', 'expensive online')
% yscale log





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

    scalar = logspace(0., -2, n);
    evals = evals ./ abs(evals) .* scalar';
    A = evecs * diag(evals) / evecs;
end


function [A, evals, evecs] = case1(n)
    k = 50;
    evals = logspace(0., -2, n);
    evals(1: k) = logspace(0., -0.01, k);
    evecs = rand_col(n, n);
    A = evecs * diag(evals) / evecs;
end


function [A, evals, evecs] = case2(n)
    k = 5;
    evals = logspace(0., -2, n);
    evals(1: k) = logspace(0., -0.01, k);
    evecs = rand_col(n, n);
    A = evecs * diag(evals) / evecs;
end


function [A, evals, evecs] = case3(n)
    k = 15;
    evals = logspace(0., -2, n);
    evals(1: k) = logspace(0., -0.01, k);
    evecs = rand_col(n, n);
    A = evecs * diag(evals) / evecs;
end


% random matrix with real eigenvalues and vectors
function [A, evals, evecs] = rand_mat_real(n)
    k = n/2;
    evals = logspace(0., -2, n);
    evals(1: k) = logspace(0., -1, k);
    evecs = rand_col(n, n);
    A = evecs * diag(evals) / evecs;
end


% random matrix with real eigenvalues and vectors
function [A, evals, evecs] = rand_mat_sym(n)
    k = n/2;
    evals = logspace(0., -2, n);
    evals(1: k) = logspace(0., -1, k);
    [evecs, ~] = qr(randn(n));
    A = evecs * diag(evals) * evecs';
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

function pred = projected_reconstruction(x0, P, Ar, steps)
    
    r = size(Ar, 1);
    [v, d] = main_eig(Ar, r);
    
    px0 = P'*x0;
    ampl = diag(v \ px0);
    evol = d .^ (0: steps-1);
    pred = P * v * ampl * evol;
end


