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
[V_on, D_on, Q_on] = extend_dmd(X, Y, r, r);
V_ons = V_on; D_ons = D_on; Q_ons = Q_on;

% current sota
[V_sota, D_sota] = direct_solve(X, Y, r);

% dmd
[V_dmd, D_dmd, Q_dmd] = extend_dmd(X, Y, r, r);

rate = 1;
B = pinv(X) * V_on;
% V == Dt*B,

% online iteration
for i = 1: steps-init
    % data preparation
    Dt = data(:, 1:init-1+i);
    D_old = Dt(:, 1:end-1);
    D_new = Dt(:, 2:end);
    % [Ux, Sx, Vx] = bksvd(D_old, 2*r);
    % A * D_old == D_new

    % alg: online
    err_subspace(i) = cond(V_ons);
    % err_subspace(i) = norm(Q_on - Q_ons, 'fro');
    % preparation
    % A_r = pinv(V_on) * D_new * pinv(D_old) * V_on;
    A_r = pinv(V_on) * D_new * B;
    % r1 = A_r;

    [vr, D_on] = main_eig(A_r, r);
    new_V_on = V_on*vr; % V_on = Dold * B, new_V_on = Dold * B * vr
    B = B * vr;
    err_vcond(i) = cond(V_ons);

    Uinv = pinv(new_V_on);

    new_BU = [B; zeros(1, r)];
    for j = 1:r

        temp1 = B * vr(:, j) / D_on(j);
        term1 = [0; temp1];
        
        lambda_minus = diag(D_on(j) - D_on);
        temp2 = B * vr * pinv(lambda_minus) * Uinv * D_new * B * vr(:, j);
        term2 = [temp2; 0];

        new_BU(:, j) = term1 + term2;   % dU == Dt*(new_BU-BU)
    end
    vecnorm_U = vecnorm(Dt * new_BU);
    B = new_BU ./ vecnorm_U;
    V_on = Dt * B;

    % [Q_on, R] = qr(V_on, 'econ');
    % err_rcond(i) = cond(R);
    % err_newu(i) = norm(V_on - Q_on * R, 'fro');
    % err_q(i) = norm(Q_on - Dt*B, 'fro');
    % err_qpinv(i) = norm(Q_on - Dt*pinv(Dt)*Q_on, 'fro');
    % 
    % err_condu(i) = cond(V_on);
    % err_uinq(i) = norm(V_on - Q_on * Q_on' * V_on, 'fro');
    % [Q_on, R, T] = qr(new_P, 'econ'); Rinv = T / R;   % Q_on = new_P * Rinv

    % record error
    % err(i) = norm(Dt*(B1 - B2), 'fro');    % why???????????????????????
    % err_ana(i) = norm(Q_on - Dt*B1, 'fro');
    % err_ref(i) = norm(Q_on - Dt*B2, 'fro');   % why???????????????????????
    % err_qr(i) = norm(Q_on - new_P * Rinv, 'fro');

    % ==========online subspace - reference version. (directly!)===========
    % At = A_org;
    % At = D_new * pinv(D_old);

    B_ons = pinv(D_old)*Q_ons;
    A_r = Q_ons' * D_new * B_ons;

    [vr, D_ons] = main_eig(A_r, r);
    V_ons = Q_ons * vr;

    % [u0, s0, v0] = svd(Dt);
    % % Now, V_ons should be in the span(D_old). is it?
    % err_span1(i) = norm(V_ons - u0 * u0' * V_ons, 'fro');
    
    for j = 1:r
        u = V_ons(:, j);

        % Q_ons == D_old * B_ons, u == D_old * B_ons * vr(:, j)
        % ======
        term1 = D_new * B_ons * vr(:, j) / D_ons(j);
        lambda_minus = diag(D_ons(j) - D_ons);
        term2 = V_ons * pinv(lambda_minus) / vr * Q_ons' * D_new * B_ons * vr(:, j);

        new_u(:, j) = term1 + term2;
    end
    [Q_ons, ~] = qr(new_u, 'econ');

    % % Now, P and Q_ons should be in the span(Dt). is it?
    % err_span3(i) = norm(P - Dt * pinv(Dt) * P, 'fro');
    % err_span4(i) = norm(Q_ons - Dt * pinv(Dt) * Q_ons, 'fro');

    % % is pinv(Dt) * Dt still almost identity?
    % err_norm_i(i) = norm(pinv(Dt) * Dt - eye(init+i), 'fro');

    % the condition number of Dt is related!
    cond_dt(i) = cond(Dt);

    % save methods error===================================================

    err_vec_on(i) = vecs_distance(V_on, V);
    err_val_on(i) = vals_distance(D_on, D);

    err_vec_onsub(i) = vecs_distance(V_ons, V);
    err_val_onsub(i) = vals_distance(D_ons, D);

    % other 2 ref
    % current sota
    [V_sota, D_sota] = main_eig(D_new * pinv(D_old), r);
    
    % dmd
    [V_dmd, D_dmd, Q_dmd] = extend_dmd(D_old, D_new, r, r);
    
    err_vec_sota(i) = vecs_distance(V_sota, V);
    err_val_sota(i) = vals_distance(D_sota, D);
    
    err_vec_dmd(i) = vecs_distance(V_dmd, V);
    err_val_dmd(i) = vals_distance(D_dmd, D);

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
timeline = 1: steps-init;
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
timeline = 1: steps-init;
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

% figure()
% hold on
% plot(err_uinq)
% yscale log
% title('QR decomposition of updated subspace error')
% 
% 
% figure()
% plot(cond_dt)
% yscale log


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


function [A, evals, evecs] = case1(n)
    cut = 10;
    evals = [logspace(0., -0.01, cut), logspace(-1, -2, n-cut)];
    evecs = rand_col(n, n);
    A = evecs * diag(evals) / evecs;
end


function [A, evals, evecs] = case2(n)
    cut = n/2;
    evals = logspace(0, -2, n);
    evals(1:cut) = logspace(0, -1, cut);
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


