clear all
rng(2024);

% matrix size
n = 400;
dim = n;
% randomly generate eigenvlaues, eigenvectors recover matrix
[A_org, evals, evecs] = rand_mat_real(n);
% [A_org, evals, evecs] = case1(n);

% n = 100;
% 
% k=-1;
% R = normrnd(0,1,n,n);
% [U,~] = qr(R,'econ');
% 
% Sigma = diag([1.00001 logspace(0,k,n-1)]);
% eig_A = diag(Sigma);
% A_org = U*Sigma*U';
% evecs = U; evals = diag(Sigma);

% initial state: a certain combination of the eigen-vectors
x0 = evecs * ones(n, 1);

% steps for the simulation snapshots
steps = 201;

% simulate snapshots observation
data = zeros([n, steps]);
data(:, 1) = x0;
for i = 2:steps
    data(:, i) = A_org * data(:, i-1);
end

% best for ref
r = 10; 
% X = data(:, 1:end-1); Y = data(:, 2:end);
% A = Y * pinv(X);

A = A_org;
[V, D] = main_eig(A, r);

% init for online subspace iteration (not practical, just for reference)
init = 2*r;
Dt = data(:, 1:init);
X = Dt(:, 1:end-1);
Y = Dt(:, 2:end);
[V_on, D_on, Q_on] = simple_dmd(X, Y, r);
B = pinv(Dt) * Q_on;
% P == Dt*B, V == P*vr == Dt*B*vr

% online iteration
for i = 1: steps-init
    % data preparation
    Dt = data(:, 1:init+i);
    D_old = Dt(:, 1:end-1);
    D_new = Dt(:, 2:end);
    % A * D_old == D_new

    % alg: online
    % preparation
    A_r = Q_on' * D_new * B;

    [vr, D_on] = main_eig(A_r, r);
    V_on = Q_on*vr;
    Vinv = vr \ Q_on';

    new_BU = [B; zeros(1, r)];
    for j = 1:r

        temp1 = B * vr(:, j) / D_on(j);
        term1 = [0; temp1];
        
        lambda_minus = diag(D_on(j) - D_on);
        temp2 = B * vr * pinv(lambda_minus) * Vinv * D_new * B * vr(:, j);
        term2 = [temp2; 0];

        new_BU(:, j) = term1 + term2;   % dU == Dt*(new_BU-BU)
    end

    new_V = Dt * new_BU;
    [Q_on, R] = qr(new_V, 'econ');
    B = new_BU / R; % theoretically, Q_on = Dt * B

    for j = 1:r
        err_vec(j, i) = vecs_distance(V_on(:, j), V(:, j));
        err_val(j, i) = vals_distance(D_on(j), D(j));
    end

    i
end

% figure()
% plot(err_ons)
% yscale log
% 
% 
% figure()
% hold on
% plot(err_span1)
% plot(err_span2)
% plot(err_span3)
% plot(err_span4)
% yscale log
% legend()

figure('Position', [100, 100, 1500, 800])
sgtitle('Eigenvalue prediction of online algorithm')
timeline = 1: steps-init;
for i = 1:9
    subplot(3,3,i)
    hold on
    plot(err_val(i, :))
    xlabel('Iteration')
    ylabel('Error')
    yscale log
    title(['Online prediction of real eigenvalue:', num2str(D(i))])
    box on
end

figure('Position', [100, 100, 1500, 800])
sgtitle('Eigenvector prediction of online algorithm')
timeline = 1: steps-init;
for i = 1:9
    subplot(3,3,i)
    hold on
    plot(err_vec(i, :))
    xlabel('Iteration')
    ylabel('Error')
    yscale log
    box on
end

% figure()
% hold on
% plot(err_ana)
% plot(err_ref)
% plot(err_qr)
% yscale log
% legend('Analytic', 'Reference', 'QR', 'location', 'best')
% title('Error caused by parameter matrix - reference form')



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
    k = n/2;
    evals = logspace(0., -2, n) .* (1 + 0.01*randn(1, n));
    evals(1: k) = logspace(0., -1, k);
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


function mat = filt_min(A)
    threshold = 1e-2;
    mat = A;
    mat(abs(A) < threshold) = 0;
end


function [A, evals, evecs] = case1(n)
    cut = 5;
    evals = [logspace(0.05, -0.05, cut), logspace(-1, -2, n-cut)];
    evecs = rand_col(n, n);
    A = evecs * diag(evals) / evecs;
end



