clear

sys_mat = diag([0.9, 1.1]);
%%% explicit eigval==0.9 and 1.1, eigvec==(1,0) and (0,1)
init_state = [5; 1];

% steps = 10;
% steps = 15;
steps = 18;
% steps = 20;
% steps = 30;
observe = [init_state, zeros(2, steps)];
for i = 1: steps
    observe(:, i+1) = sys_mat * observe(:, i);
end

%%% simple DMD
X = observe(:, 1: end-1); Y = observe(:, 2: end);
[mode, eigenvalue] = dmd_decom(X, Y, 1);
recon = dmd_recon(mode, eigenvalue, init_state, 0, steps);

%%% try to fix eigs
cut = 10;   %%% set this parameter
X = observe(:, 1: cut); Y = observe(:, 2: cut+1);
[mode1, eigenvalue1] = dmd_decom(X, Y, 1);

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

mode1 = mode1 ./ vecnorm(mode1);
recon1 = dmd_recon(mode1, eigenvalue1, init_state, 0, steps);

x_real = observe(1, :); y_real = observe(2, :);
x_rec = recon(1, :); y_rec = recon(2, :);
x_rec1 = recon1(1, :); y_rec1 = recon1(2, :);

figure()
sz = 25; c = linspace(1, 10, length(x_rec));
hold on
scatter(x_real, y_real, sz, c, "filled");
scatter(x_rec, y_rec, sz, c, "o");
scatter(x_rec1, y_rec1, sz, c, "*")
legend('real',...
    'recon',...
    'recon1')

cut = 10;   %%% set this parameter
X = observe(:, 1: cut); Y = observe(:, 2: cut+1);
[mode2, eigenvalue2] = dmd_decom(X, Y, 1);

for i = cut+1: steps
    state = observe(:, steps); next_state = observe(:, steps+1);
    pred_state = dmd_recon(mode2, eigenvalue2, state, 1, 1);
    delta_state = next_state - pred_state;

    %%% approx mode
    delta_mode = delta_state * pinv(state) * mode2 * diag(1 ./ eigenvalue2);
    %%% approx value
    delta_value = pinv(mode2) * delta_state * pinv(state) * mode2;
    delta_value = diag(delta_value);

    %%% add together
    mode2 = mode2 + delta_mode;
    eigenvalue2 = eigenvalue2 + delta_value;
end

mode2 = mode2 ./ vecnorm(mode2);
recon2 = dmd_recon(mode2, eigenvalue2, init_state, 0, steps);

x_real = observe(1, :); y_real = observe(2, :);
x_rec = recon(1, :); y_rec = recon(2, :);
x_rec2 = recon2(1, :); y_rec2 = recon2(2, :);

figure()
sz = 25; c = linspace(1, 10, length(x_rec));
hold on
scatter(x_real, y_real, sz, c, "filled");
scatter(x_rec, y_rec, sz, c, "o");
scatter(x_rec2, y_rec2, sz, c, "*")
legend('real',...
    'recon',...
    'recon1')
