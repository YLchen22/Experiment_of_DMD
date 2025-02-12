clear all
rng(2024);

% matrix size
n = 300;
dim = n;
% randomly generate eigenvlaues, eigenvectors recover matrix
[A_org, evals, evecs] = rand_mat_real(n);
% [A_org, evals, evecs] = case1(n);

x0 = evecs * ones(n, 1);

% steps for the simulation snapshots
steps = 101;

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

Dt = data;
X = Dt(:, 1:end-1);
Y = Dt(:, 2:end);

% current sota
[V_sota, D_sota] = direct_solve(X, Y, r);

% dmd
[V_dmd, D_dmd, Q_dmd] = simple_dmd(X, Y, r);

k = lsqr(X, Y(:, end), [], steps);
norm(Y(:, end) - X*k, 'fro')

[V_on, D_on, Q_on] = simple_dmd(X, Y, r);
vr = Q_on' * V_on;
B = pinv(X) * Q_on;
% P == X*B

err_vec_on = vecs_distance(V_on, V);
err_val_on = vals_distance(D_on, D);

err_vec_onsub = err_vec_on;
err_val_onsub = err_val_on;

err_vec_sota = vecs_distance(V_sota, V);
err_val_sota = vals_distance(D_sota, D);

err_vec_dmd = vecs_distance(V_dmd, V);
err_val_dmd = vals_distance(D_dmd, D);

iter = 400;
% online iteration
for i = 1: iter

    % alg: online

    % preparation
    % A_r = Q_on' * D_new * pinv(D_old) * Q_on;
    A_r = Q_on' * Y * B;
    r1 = A_r;
    % N = Q_on' * D_old; A_r = Q_on' * D_new * pinv(N);

    [vr, D_on] = main_eig(A_r, r);
    V_on = Q_on*vr;

    Uinv = vr \ Q_on';  % V_on = Dt * B * vr
    BU = B * vr;
    new_BU = [B; zeros(1, r)];
    for j = 1:r

        temp1 = B * vr(:, j) / D_on(j);
        term1 = [0; temp1];
        
        lambda_minus = diag(D_on(j) - D_on);
        temp2 = B * vr * pinv(lambda_minus) * Uinv * Y * B * vr(:, j);
        term2 = [temp2; 0];

        new_BU(:, j) = term1 + term2;   % dU == Dt*(new_BU-BU)
    end

    block1 = new_BU(1:end-1, :); block2 = new_BU(end, :);
    new_BU = block1 + k * block2;

    new_U = X * new_BU;

    [Q_on, R] = qr(new_U, 'econ'); % Q_on = Dt * new_BU / R
    % B = pinv(Dt) * Q_on;
    B = new_BU / R;

    % [Q_on, R, T] = qr(new_P, 'econ'); Rinv = T / R;   % Q_on = new_P * Rinv

    % record error
    % err(i) = norm(Dt*(B1 - B2), 'fro');    % why???????????????????????
    % err_ana(i) = norm(Q_on - Dt*B1, 'fro');
    % err_ref(i) = norm(Q_on - Dt*B2, 'fro');   % why???????????????????????
    % err_qr(i) = norm(Q_on - new_P * Rinv, 'fro');

    % ==========online subspace - reference version. (directly!)===========
    V_ons = V_on; D_ons = D_on;

    err_vec_onsub = [err_vec_onsub, vecs_distance(V_ons, V)];
    err_val_onsub = [err_val_onsub, vals_distance(D_ons, D)];

    % save methods error===================================================

    err_vec_on = [err_vec_on, vecs_distance(V_on, V)];
    err_val_on = [err_val_on, vals_distance(D_on, D)];

    % other 2 ref
    % current sota
    [V_sota, D_sota] = main_eig(Y * pinv(X), r);
    
    % dmd
    [V_dmd, D_dmd, Q_dmd] = simple_dmd(Y, X, r);
    
    err_vec_sota = [err_vec_sota, vecs_distance(V_sota, V)];
    err_val_sota = [err_val_sota, vals_distance(D_sota, D)];
    
    err_vec_dmd = [err_vec_dmd, vecs_distance(V_dmd, V)];
    err_val_dmd = [err_val_dmd, vals_distance(D_dmd, D)];

    % error by index
    for j = 1:r
        recval_on(j, i) = abs(D_on(j) - D(j));
        recval_ons(j, i) = abs(D_ons(j) - D(j));
        recval_sota(j, i) = abs(D_sota(j) - D(j));
        recval_dmd(j, i) = abs(D_dmd(j) - D(j));

        recvec_on(j, i) = vecs_distance(V_on(:, j), V(:, j));
        recvec_ons(j, i) = vecs_distance(V_ons(:, j), V(:, j));
        recvec_sota(j, i) = vecs_distance(V_sota(:, j), V(:, j));
        recvec_dmd(j, i) = vecs_distance(V_dmd(:, j), V(:, j));
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
timeline = 1: iter;
for i = 1:9
    subplot(3,3,i)
    hold on
    plot(recval_on(i, :))
    plot(recval_ons(i, :))
    plot(recval_dmd(i, :))
    plot(recval_sota(i, :))
    xlabel('Iteration')
    ylabel('Error')
    yscale log
    legend('Online', 'Online-reference', 'DMD', 'Direct regression', 'Location', 'southwest')
    title(['Online prediction of real eigenvalue:', num2str(D(i))])
    box on
end

figure('Position', [100, 100, 1500, 800])
sgtitle('Eigenvector prediction of online algorithm')
timeline = 1: iter;
for i = 1:9
    subplot(3,3,i)
    hold on
    plot(recvec_on(i, :))
    plot(recvec_ons(i, :))
    plot(recvec_dmd(i, :))
    plot(recvec_sota(i, :))
    xlabel('Iteration')
    ylabel('Error')
    yscale log
    legend('Online', 'Online-reference', 'DMD', 'Direct regression', 'Location', 'southwest')
    title(['Online prediction of real eigenvalue:', num2str(D(i))])
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


figure('Position', [100, 100, 1200, 600])
subplot(121)
box on
hold on
plot(err_vec_on)
plot(err_vec_onsub)
plot(err_vec_dmd)
plot(err_vec_sota)
yscale log
legend('Online algorithm', 'Computational online algorithm', 'DMD', 'SOTA', 'location', 'best')
subtitle('Error of eigenvector')

subplot(122)
box on
hold on
plot(err_val_on)
plot(err_val_onsub)
plot(err_val_dmd)
plot(err_val_sota)
yscale log
legend('Online algorithm', 'Computational online algorithm', 'DMD', 'SOTA', 'location', 'best')
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
    evals = logspace(0., -2, n);
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
    [V, D] = main_eig(R, r);
    Q = u;
    V = Q * V;
end


function [V, D, Q] = extend_dmd(X, Y, r, ex)
    [u, s, v] = svds(X, r + ex);
    R = u' * Y * v * pinv(s);
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
    cut = 10;
    evals = [logspace(0.001, -0.001, cut), logspace(-1, -2, n-cut)];
    evecs = rand_col(n, n);
    A = evecs * diag(evals) / evecs;
end


function [U, S, V] = randomized_svd(A, k, p, q)
    % 随机化 SVD 计算 A ≈ U S V'
    % A: 输入矩阵 (m x n)
    % k: 目标奇异值数
    % p: 过采样参数 (默认为 10)
    % q: 迭代次数 (默认为 1)
    
    if nargin < 3
        p = 10; % 过采样参数
    end
    if nargin < 4
        q = 1; % 默认迭代次数
    end
    
    [m, n] = size(A);
    l = k + p;  % 增加过采样，提高精度
    
    % Step 1: 生成随机高斯矩阵
    Omega = randn(n, l);
    
    % Step 2: 计算投影矩阵 Y = A * Omega
    Y = A * Omega;
    
    % Step 3: 迭代增强 (Power Iteration, 适用于低秩矩阵)
    for i = 1:q
        [Y, ~] = qr(Y, 0); % 正交化
        Y = A' * Y;
        [Y, ~] = qr(Y, 0);
        Y = A * Y;
    end
    
    % Step 4: 计算正交基 Q
    [Q, ~] = qr(Y, 0);
    
    % Step 5: 计算小矩阵 B
    B = Q' * A;
    
    % Step 6: 对 B 进行 SVD
    [U_tilde, S, V] = svd(B, 'econ');
    
    % Step 7: 计算最终 U
    U = Q * U_tilde;
    
    % 只返回前 k 个奇异值
    U = U(:, 1:k);
    S = S(1:k, 1:k);
    V = V(:, 1:k);
end


function [U,S,V] = bksvd(A, k, iter, bsize, center)
%--------------------------------------------------------------------------
% Randomized block Krylov iteration for truncated singular value decomposition
% Computes approximate top singular vectors and corresponding values
% Described in Musco, Musco, 2015 (http://arxiv.org/abs/1504.05477)
%
% usage : 
%
%  input:
%  * A : matrix to decompose
%  * k : number of singular vectors to compute, default = 6
%  * iter : number of iterations, default = 3
%  * bsize : block size, must be >= k, default = k
%  * center : set to true if A's rows should be mean centered before the
%  singular value decomposition (e.g. when performing principal component 
%  analysis), default = false
%
%
%  output:
%  k singular vector/value pairs. 
%  * U : a matrix whose columns are approximate top left singular vectors for A
%  * S : a diagonal matrix whose entries are A's approximate top singular values
%  * V : a matrix whose columns are approximate top right singular vectors for A
%
%  U*S*V' is a near optimal rank-k approximation for A
%--------------------------------------------------------------------------

% Check input arguments and set defaults.
if nargin > 5
    error('bksvd:TooManyInputs','requires at most 5 input arguments');
end
if nargin < 1
    error('bksvd:TooFewInputs','requires at least 1 input argument');
end
if nargin < 2
    k = 6;
end
k = min(k,min(size(A)));

if nargin < 3
    iter = 3;
end
if nargin < 4
    bsize = k;
end
if nargin < 5
    center = false;
end
if(k < 1 || iter < 1 || bsize < k)
    error('bksvd:BadInput','one or more inputs outside required range');
end

% Calculate row mean if rows should be centered.
u = zeros(1,size(A,2));
if(center)
    u = mean(A);
end
l = ones(size(A,1),1);


% We want to iterate on the smaller dimension of A.
[n, ind] = min(size(A));
tpose = false;
if(ind == 1) 
    tpose = true;
    l = u'; u = ones(1,size(A,1));
    A = A';
end

% Allocate space for Krylov subspace.
K = zeros(size(A,2),bsize*iter);
% Random block initialization.
block = randn(size(A,2),bsize);
[block,R] = qr(block,0);

% Preallocate space for temporary products.
T = zeros(size(A,2),bsize);

% Construct and orthonormalize Krlov Subspace. 
% Orthogonalize at each step using economy size QR decomposition.
for i=1:iter
    T = A*block - l*(u*block);
    block = A'*T - u'*(l'*T);
    [block,R] = qr(block,0);
    K(:,(i-1)*bsize+1:i*bsize) = block;
end
[Q,R] = qr(K,0);

% Rayleigh-Ritz postprocessing with economy size dense SVD.
T = A*Q - l*(u*Q);

[Ut,St,Vt] = svd(T,0);
S = St(1:k,1:k);
if(~tpose)
    U = Ut(:,1:k);
    V = Q*Vt(:,1:k);
else
    V = Ut(:,1:k);
    U = Q*Vt(:,1:k);
end

end


function A_filt = filt_cond(A, k)
    [u, s, v] = bksvd(A, k);
    A_filt = u * s * v';
end


