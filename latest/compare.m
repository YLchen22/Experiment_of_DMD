clear all
close all

%% hyper parameters
n = 300;    % dimension
steps = 201;    % data size
init = 21;  % initial step for the alg
r = 10;     % lower-rank


%% generate matrix and data
[A_org, evals, evecs] = rand_mat_real(n);
% [A_org, evals, evecs] = case1(n);

x0 = evecs * ones(n, 1);
% simulate snapshots observation
data = zeros([n, steps]);
data(:, 1) = x0;
for i = 2:steps
    data(:, i) = A_org * data(:, i-1);
end

% reference result
[V, D] = main_eig(A_org, r);

[evals_on, vr_on, P_on, B_on] = online_iteration(data, init, r, 'cheap');
[evals_ex, vr_ex, P_ex, B_ex] = online_iteration(data, init, r, 'expensive');

[evals_dmd, evecs_dmd, P_dmd] = ref_by_steps(data, init, r, 'dmd');
[evals_dr, evecs_dr, P_dr] = ref_by_steps(data, init, r, 'direct');


%% separated visualization
figure('Position', [0, 100, 1500, 800])
sgtitle('Eigenvalue prediction of online algorithm')
timeline = init: steps;
for i = 1:9
    subplot(3,3,i)
    hold on

    % calculate used error
    eval_ref = D(i);
    for j = timeline
        err_on(j) = vals_distance(eval_ref, evals_on{j}(i));
        err_ex(j) = vals_distance(eval_ref, evals_ex{j}(i));
        err_dmd(j) = vals_distance(eval_ref, evals_dmd{j}(i));
        err_dr(j) = vals_distance(eval_ref, evals_dr{j}(i));
    end
    
    plot(err_on(timeline))
    plot(err_ex(timeline))
    plot(err_dmd(timeline))
    plot(err_dr(timeline))
    xlabel('Iteration')
    ylabel('Error')
    yscale log
    legend('Online', 'Online-expensive', 'DMD', 'Direct regression', 'Location', 'southwest')
    title(['Online prediction of real eigenvalue:', num2str(D(i))])
    box on
end


figure('Position', [0, 100, 1500, 800])
sgtitle('Eigenvector prediction of online algorithm')
timeline = init: steps;
for i = 1:9
    subplot(3,3,i)
    hold on

    % calculate used error
    evec_ref = V(:, i);
    for j = timeline
        err_on(j) = vecs_distance(evec_ref, P_on{j}*vr_on{j}(:, i));
        err_ex(j) = vecs_distance(evec_ref, P_ex{j}*vr_ex{j}(:, i));
        err_dmd(j) = vecs_distance(evec_ref, evecs_dmd{j}(:, i));
        err_dr(j) = vecs_distance(evec_ref, evecs_dr{j}(:, i));
    end
    
    plot(err_on(timeline))
    plot(err_ex(timeline))
    plot(err_dmd(timeline))
    plot(err_dr(timeline))
    xlabel('Iteration')
    ylabel('Error')
    yscale log
    legend('Online', 'Online-expensive', 'DMD', 'Direct regression', 'Location', 'southwest')
    box on
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


% random matrix with real eigenvalues and vectors
function [A, evals, evecs] = rand_mat_real(n)
    k = n/2;
    evals = logspace(0., -2, n);
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

function err = vecs_distance(vecs1, vecs2)
    prod = vecs1' * vecs2;
    diag_err = diag(abs(prod) - eye(size(prod)));
    err = norm(diag_err, 'fro');
end

function err = vals_distance(vals1, vals2)
    err = norm(vals1 - vals2, 'fro');
end