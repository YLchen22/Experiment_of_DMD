clear
rng(2024)

%%% real case
decays = sort(0.1 * rand(1, 5), 'descend');
eig_mat = diag([1.2, 1, 0.9, 0.8, 0.3, decays]);
vec_mat = randn(10); vec_mat = vec_mat ./ vecnorm(vec_mat);
sys_mat = vec_mat * eig_mat / vec_mat;

%%% explicit eigval as eig_mat, and eigvec as vec_mat
init_state = unifrnd(2, 3, [10, 1]);

steps = 50;

r = 3;

observe = [init_state, zeros(10, steps)];
for i = 1: steps
    observe(:, i+1) = sys_mat * observe(:, i);
end

%%% simple DMD
X = observe(:, 1: end-1); Y = observe(:, 2: end);
[mode, eigenvalue] = dmd_decom(X, Y, r);
recon = dmd_recon(mode, eigenvalue, init_state, 0, steps);

%%% try to fix eigs
cut = 30;   %%% set this parameter
X = observe(:, 1: cut); Y = observe(:, 2: cut+1);
[mode1, eigenvalue1] = dmd_decom(X, Y, r);

for i = cut+1: steps
    state = observe(:, i); next_state = observe(:, i+1);
    pred_state = dmd_recon(mode1, eigenvalue1, state, 1, 1);
    delta_state = next_state - pred_state;

    %%% approx mode
    delta_mode = delta_state * pinv(state) * mode1 * diag(1 ./ eigenvalue1);
    %%% approx value
    delta_value = pinv(mode1) * delta_state * pinv(state) * mode1;
    delta_value = diag(delta_value);

    %%% add together
    mode1 = mode1 + delta_mode;
    eigenvalue1 = eigenvalue1 + delta_value;
end

recon1 = dmd_recon(mode1, eigenvalue1, init_state, 0, steps);

x_real = observe(1, :); y_real = observe(2, :);
x_rec = recon(1, :); y_rec = recon(2, :);
x_rec1 = recon1(1, :); y_rec1 = recon1(2, :);

