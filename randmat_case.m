clear all
close all
rng(2024)

%% hyper parameters
n = 400;    % dimension
steps = 200;    % data size
r = 9;     % lower-rank
init = steps/2;  % initial step for the alg
w = init;

%% generate matrix and data
[A_org, evals, evecs] = rand_mat_real(n);
% [A_org, evals, evecs] = rand_mat_sym(n);
% [A_org, evals, evecs] = case3(n);
% [A_org, evals, evecs] = case2(n);
% [A_org, evals, evecs] = rand_mat(n);

evecs_ = evecs; evals_ = evals;


distri = ones(n, 1);
x0 = evecs * distri;
% control the beginning distribution of eigenvectors,
% and simulate snapshots observation

data = zeros([n, steps]);
data(:, 1) = x0;
for i = 2:steps
    data(:, i) = A_org * data(:, i-1);
end

% for i = 1:steps
%     fcond(i) = cond(data(:, 1:i));
%     bcond(i) = cond(data(:, i:end));
% end
% figure()
% hold on
% plot(fcond)
% plot(bcond)
% yscale log
% 
% for i = 1:steps-w
%     wcond(i) = cond(data(:, i:i+w));
% end
% figure()
% hold on
% plot(wcond)
% yscale log


% data = data + 1e-8 * eye(size(data));
data = data + 1e-6 * randn(size(data));     % Robust!

% [qall, ~] = qr(evecs);
% [qsub, ~] = qr(evecs(:, 1:r));
% for i = 1:steps
%     dt = data(:, 1:i);
%     errall(i) = norm(dt - qall * qall' * dt, 'fro');
%     errsub(i) = norm(dt - qsub * qsub' * dt, 'fro');
% end
% figure()
% hold on
% plot(errall)
% plot(errsub)
% yscale log


% Z = data;
% for i = 2:steps
%     x = lsqr(Z(:, 1:i-1), Z(:, i), 1e-8, n);
%     Z(:, i) = Z(:, i) - Z(:, 1:i-1) * x;
%     c(i) = cond(Z(:, 1:i));
%     cd(i) = cond(data(:, 1:i));
% end
% figure()
% hold on
% plot(c)
% plot(cd)
% yscale log
% 

% win = 10;
% [q, ~] = qr(evecs(:, 1:win), 'econ');
% clear c
% for i = win: steps
%     Dt = data(:, i-win+1: i);
%     c(i) = cond(q * q' * Dt);
% end
% figure()
% plot(c)

% reference result
V = evecs(:, 1:r); D = evals(1:r);
% [V, D] = main_eig(A_org, r);
[Q, Q2V] = qr(evecs);
cond(Q2V(:, 1:r))

[evals_on, vr_on, P_on, B_on] = online_iteration_windowed(data, init, r, w);
% [evals_on, vr_on, P_on, B_on] = subspace_iteration(data, r, init, steps);

[evals_ex, vr_ex, P_ex, B_ex] = online_iteration(data, init, r, 'cheap');

[evals_dmd, vr_dmd, P_dmd] = baseline_simple(data, init, r, 'dmd');
% [evals_dr, evecs_dr, P_dr] = baseline_simple(data, init, r, 'direct');

% [evals_dmd, evecs_dmd, P_dmd] = baseline_windowed(data, init, r, w, 'dmd');
[evals_dr, vr_dr, P_dr] = baseline_windowed(data, init, r, w, 'direct');


new_x0 = evecs * (rand(n, 1)*2 + 1);
% new_x0 = data(:, 1);
% new_x0 = data(:, end);
t = 20000;
p_data = zeros([n, t]);
p_data(:, 1) = new_x0; 
pred_real = p_data; pred_on = p_data; pred_dmd = p_data;
A_real = V * diag(D) * pinv(V); A_real = data(:, 2:end) * pinv(data(:, 1: end-1));
A_on = P_on{end} * vr_on{end} * diag(evals_on{end}) / vr_on{end} * P_on{end}';
A_dmd = P_dmd{end} * vr_dmd{end} * diag(evals_dmd{end}) / vr_dmd{end} * P_dmd{end}';
for i = 2:t
    p_data(:, i) = A_org * p_data(:, i-1);
    pred_real(:, i) = A_real * pred_real(:, i-1);
    pred_on(:, i) = A_on * pred_on(:, i-1);
    pred_dmd(:, i) = A_dmd * pred_dmd(:, i-1);
end

% P = P_on; vr = vr_on; evals = evals_on;
% ampl = diag(vr{end} \ P{end}' * new_x0);      % r*r
% evol = evals{end} .^ (0: t-1);            % a trick to define vandermonde matrix
% pred_on = real(P{end} * vr{end} * ampl * evol);
% 
% P = P_dmd; vr = vr_dmd; evals = evals_dmd;
% ampl = diag(vr{end} \ P{end}' * new_x0);      % r*r
% evol = evals{end} .^ (0: t-1);            % a trick to define vandermonde matrix
% pred_dmd = real(P{end} * vr{end} * ampl * evol);
% 
% [P, ~] = qr(V, 'econ'); vr = P'*V; evals = evals_';
% P = {P,}; vr = {vr,}; evals = {evals(1:r)};
% ampl = diag(vr{end} \ P{end}' * new_x0);      % r*r
% evol = evals{end} .^ (0: t-1);            % a trick to define vandermonde matrix
% pred_real = real(P{end} * vr{end} * ampl * evol);

err_real = vecnorm(p_data - pred_real);
err_on = vecnorm(p_data - pred_on);
err_dmd = vecnorm(p_data - pred_dmd);
figure()
hold on
plot(err_real)
plot(err_on)
plot(err_dmd)
yscale log
legend('Recover error of real evecs', 'Recover error of online', 'Recover error of DMD')




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
        err_dmd(j) = vecs_distance(evec_ref, P_dmd{j}*vr_dmd{j}(:, i));
        err_dr(j) = vecs_distance(evec_ref, P_dr{j}*vr_dr{j}(:, i));
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


figure('Position', [0, 100, 1500, 800])
hold on
sgtitle('Difference of Subspace')
[Q, ~] = qr(evecs(:, 1:r));
v = evecs(:, 1:r);
timeline = init: steps;
for j = timeline
    % err_on(j) = norm((eye(n) - Q*Q')*P_on{j}, 'fro');
    % err_ex(j) = norm((eye(n) - Q*Q')*P_ex{j}, 'fro');
    % err_dmd(j) = norm((eye(n) - Q*Q')*P_dmd{j}, 'fro');
    % err_dr(j) = norm((eye(n) - Q*Q')*P_dr{j}, 'fro');

    err_on(j) = norm((eye(n) - P_on{j} * P_on{j}') * v, 'fro');
    err_ex(j) = norm((eye(n) - P_ex{j} * P_ex{j}') * v, 'fro');
    err_dmd(j) = norm((eye(n) - P_dmd{j} * P_dmd{j}') * v, 'fro');
    err_dr(j) = norm((eye(n) - P_dr{j} * P_dr{j}') * v, 'fro');
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


figure('Position', [0, 100, 1500, 800])
hold on
sgtitle('Difference of Each Projected Real Eig-vector')
timeline = init: steps;
for i = 1:r
    subplot(3,3,i)
    hold on
    for j = timeline
        v = evecs(:, i);
        err_on(i, j) = norm((eye(n) - P_on{j} * P_on{j}') * v, 'fro');
        err_ex(i, j) = norm((eye(n) - P_ex{j} * P_ex{j}') * v, 'fro');
        err_dmd(i, j) = norm((eye(n) - P_dmd{j} * P_dmd{j}') * v, 'fro');
        err_dr(i, j) = norm((eye(n) - P_dr{j} * P_dr{j}') * v, 'fro');
    end
    plot(err_on(i, timeline))
    plot(err_ex(i, timeline))
    plot(err_dmd(i, timeline))
    plot(err_dr(i, timeline))
    legend('online', 'expensive online', 'dmd', 'direct', 'Location', 'southwest')
    xlabel('Iteration')
    ylabel('Error')
    yscale log
    subtitle(['error of the ', num2str(i), '-th vector'])
    box on
end

evals_on{end}
evals_dr{end}


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

