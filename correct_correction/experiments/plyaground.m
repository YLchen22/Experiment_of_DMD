clear
rng(2024);

% matrix size
dim = 400;
step = 50;
period = 200;

[Xdata, Ydata, data] = circle_norm_data_generator(dim, step, period);
dmin = min(data, [], 'all'); dmax = max(data, [], 'all');
X = data(:, 1:end-1); Y = data(:, 2:end);

% best for ref
A = Y * pinv(X);
r = 10;
[V, D] = main_eig(A, r);
A_real = V * diag(D) * pinv(V);

% dmd
[U, S, V_] = svds(X, r);
R = U' * Y * V_ / S;
A_dmd = U * R * U';
[V_dmd, D_dmd] = main_eig(R, r);
V_dmd = U * V_dmd;

recover_data = [];
state = data(:, 1);
figure(1)
for t = 1: 200
    recover_data(:, t) = state;

    Zdata = reshape(state, dim^0.5, dim^0.5);

    contourf(Xdata, Ydata, real(Zdata));
    colormap(jet);
    colorbar;
    clim([dmin, dmax]);

    pause(0.02)
    state = A_real * state;
end

approx_data = [];
state = data(:, 1);
figure(2)
for t = 1: 200
    approx_data(:, t) = state;

    Zdata = reshape(state, dim^0.5, dim^0.5);

    contourf(Xdata, Ydata, real(Zdata));
    colormap(jet);
    colorbar;
    clim([dmin, dmax]);

    pause(0.02)
    state = A_dmd * state;
end


for i = 1:r
    vec_err(i) = vecs_distance(V(:, i), V_dmd(:, i));
    val_err(i) = vals_distance(D(i), D_dmd(i));
end
figure(3)
subplot(121)
plot(vec_err)
yscale log
subplot(122)
plot(abs(D))
hold on
plot(abs(D_dmd))
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
