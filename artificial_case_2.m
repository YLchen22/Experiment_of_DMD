clear all
close all
rng(2024)

%% algorithm hyper-parameters
n = 400;    % dimension
steps = 110;    % data size
veri_steps = steps*2; % preserve some data for verification
r = 10;     % lower-rank
init = steps-r;  % initial step for the alg
w = init;

case_name = 'artificial';

%% test case settings
noise = 1e-2;
split = 15;
dt = 0.01;
k = 2;
k=k-1;

exp_re = [zeros(1, k), linspace(-0., -0.5, split-k)/dt];
exp_im = linspace(0.02*pi, (0+2/3)*pi, split) / dt;

% %% test case settings
% noise = 0.01;
% split = 30;
% dt = 0.01;
% k = 3;
% k=k-1;
% 
% exp_re = [zeros(1, k), linspace(0, -0.5, split-k)/dt];
% exp_im = linspace(0.05*pi, (0+2/3)*pi, split) / dt;

%% generate data!
% artificial system and mat
A1 = [];
for ii = 1: split
    A2 = [exp_re(ii), exp_im(ii); -exp_im(ii) exp_re(ii)];
    A1 = [A1, A2];
end
Alowrank = [];
for ii = 1: split
    Alowrank = blkdiag(Alowrank,A1(:,(ii-1)*2+1: ii*2));
end

% handle to low-dimensional operator for simulations
dynsys = @(t,x) Alowrank*x;  % handle to operator

% construct map from low-rank to full state-dimension (Q: X^r --> X^n)
[Q, ~] = qr(randn(n, 2*split), 'econ');
Qrand = randn(n, 2*split);
Q = Qrand ./ vecnorm(Qrand, 2, 1);

% record (if analytic) sys-mat====================
Ak = expm(Alowrank*dt);
[Vall, Dall] = eig(Ak); Dall = diag(Dall);
[V, Dk] = main_eig(Ak, min([split*2, r]));
Vk = Q*V;
% ================================================

x0 = ones(2*split, 1);
t = (0: (steps+veri_steps)-1) * dt;
[~, rdata] = ode45(dynsys,t,x0);
rdata = rdata';

ex_data = Q * rdata;
data = ex_data(:, 1: steps) + noise*randn(n, steps); 
clean_data = ex_data(:, 1: steps) + 1e-6*randn(n, steps); 
veri_data = ex_data(:, steps+1: steps+veri_steps);
X = data(:, 1:end-1); Y = data(:, 2:end);

disp(['Mean absolute value of clean data:', num2str(mean(abs(clean_data), 'all'))])

%% implement our methods!
[evals_on, vr_on, P_on, B_on] = online_iteration_windowed(data, init, r, w);

[evals_dmd, vr_dmd, P_dmd] = baseline_simple(data, init, r, 'dmd');
[evals_td, vr_td, P_td] = baseline_simple(data, init, r, 'tdmd');
% [evals_dmd, vr_dmd, P_dmd] = baseline_windowed(data, init, r, w, 'dmd');

thetas = 0: 0.01: 2*pi;

%% data-based visualizations ================================================================

%% 1. residual map
[num_lin, res_map] = field_residual(evals_dmd{end}, X, Y);
[re_map, im_map] = meshgrid(num_lin);
max_res = max(res_map, [], 'all');
min_res = min(res_map, [], 'all');
%
id = init;
Pdmd = P_dmd{id}; vrdmd = vr_dmd{id}; Ddmd = evals_dmd{id};
Ptd = P_td{id}; vrtd = vr_td{id}; Dtd = evals_td{id};
Pon = P_on{id}; vron = vr_on{id}; Don = evals_on{id};

figure()
hold on
pcolor(re_map, im_map, res_map)
colormap sky
hold on
shading interp
cb = colorbar(); cb.Label.String = 'Residual';
contour(re_map, im_map, res_map, linspace(min_res, max_res, 8), '-k', LineWidth=0.5)
% scref = scatter(real(Dk), imag(Dk), 64, 'blue', 'o', 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
scdmd = scatter(real(Ddmd), imag(Ddmd), 'square', 'MarkerFaceAlpha', 0., 'LineWidth', 1.);
sctd = scatter(real(Dtd), imag(Dtd), 'green', 'diamond', 'MarkerFaceAlpha', 0., 'LineWidth', 1.);
scon= scatter(real(Don), imag(Don), 'red', '+', 'LineWidth', 1.5);

uc = plot(cos(thetas), sin(thetas), 'k--', LineWidth=1.5);
legend([scdmd, sctd, scon, uc], {'DMD', 'TDMD', 'FODMD', 'Unit circle'}, Location='southwest')
title('Residual map and the spectrum results - Before')
xlabel('Re')
ylabel('Im')
axis equal

% ================================
id = steps;
Pdmd = P_dmd{id}; vrdmd = vr_dmd{id}; Ddmd = evals_dmd{id};
Ptd = P_td{id}; vrtd = vr_td{id}; Dtd = evals_td{id};
Pon = P_on{id}; vron = vr_on{id}; Don = evals_on{id};

figure()
hold on
pcolor(re_map, im_map, res_map)
colormap sky
hold on
shading interp
cb = colorbar(); cb.Label.String = 'Residual';
contour(re_map, im_map, res_map, linspace(min_res, max_res, 8), '-k', LineWidth=0.5)
% scref = scatter(real(Dk), imag(Dk), 64, 'blue', 'o', 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
scdmd = scatter(real(Ddmd), imag(Ddmd), 'square', 'MarkerFaceAlpha', 0., 'LineWidth', 1.);
sctd = scatter(real(Dtd), imag(Dtd), 'green', 'diamond', 'MarkerFaceAlpha', 0., 'LineWidth', 1.);
scon= scatter(real(Don), imag(Don), 'red', '+', 'LineWidth', 1.5);

uc = plot(cos(thetas), sin(thetas), 'k--', LineWidth=1.5);
legend([scdmd, sctd, scon, uc], {'DMD', 'TDMD', 'FODMD', 'Unit circle'}, Location='southwest')
title('Residual map and the spectrum results - After')
xlabel('Re')
ylabel('Im')
axis equal

%% 2. phase - modulus
id = steps;
Pdmd = P_dmd{id}; vrdmd = vr_dmd{id}; Ddmd = evals_dmd{id};
Ptd = P_td{id}; vrtd = vr_td{id}; Dtd = evals_td{id};
Pon = P_on{id}; vron = vr_on{id}; Don = evals_on{id};
figure()
hold on
% scref = scatter(real(Dk), imag(Dk), 64, 'blue', 'o', 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
scdmd = scatter(angle(Ddmd), abs(Ddmd), 'square', 'MarkerFaceAlpha', 0., 'LineWidth', 1.);
sctd = scatter(angle(Dtd), abs(Dtd), 'green', 'diamond', 'MarkerFaceAlpha', 0., 'LineWidth', 1.);
scon= scatter(angle(Don), abs(Don), 'red', '+', 'LineWidth', 1.5);
legend('DMD', 'TDMD', 'FODMD', Location='best')
title('Phase-modulus diagram of the final spectrums')
xlabel('Phase')
ylabel('Modulus')
yscale log
axis padded
box on

%% 3. prediction error
id = steps;
Pdmd = P_dmd{id}; vrdmd = vr_dmd{id}; Ddmd = evals_dmd{id};
Ptd = P_td{id}; vrtd = vr_td{id}; Dtd = evals_td{id};
Pon = P_on{id}; vron = vr_on{id}; Don = evals_on{id};

x0 = veri_data(:, 1);
pred_on = make_prediction(Pon, vron, Don, veri_steps, x0);
pred_dmd = make_prediction(Pdmd, vrdmd, Ddmd, veri_steps, x0);
pred_td = make_prediction(Ptd, vrtd, Dtd, veri_steps, x0);

err_on = vecnorm(veri_data - pred_on, 2, 1);
err_td = vecnorm(veri_data - pred_td, 2, 1);
err_dmd = vecnorm(veri_data - pred_dmd, 2, 1);
figure()
hold on
% plot(err_dmd, Color='green')
% plot(err_td, Color='blue')
% plot(err_on, LineWidth=1.5, Color='Red')
plot(err_dmd)
plot(err_td)
plot(err_on, LineWidth=1.5)
legend('DMD', 'TDMD', 'FODMD', Location='northwest')
title('Squared Error of Future States Prediction')
xlabel('Future steps')
ylabel('Squared error')
xlim tight

figname = 'prederr';
saveas(gcf, [case_name, '_', figname, '.png']);

%% 4. calculate db-times residual
X = clean_data(:, 1:end-1); Y = clean_data(:, 2:end);
timeline = init:steps;
for i = timeline
    res_dmd(i) = sum_residual(evals_dmd{i}, P_dmd{i}*vr_dmd{i}, X, Y);
    res_td(i) = sum_residual(evals_td{i}, P_td{i}*vr_td{i}, X, Y);
    res_on(i) = sum_residual(evals_on{i}, P_on{i}*vr_on{i}, X, Y);
end
figure()
hold on
plot(timeline-init, res_dmd(timeline), '-*')
plot(timeline-init, res_td(timeline), '-*')
plot(timeline-init, res_on(timeline), '-o', LineWidth=1.5)
legend('DMD', 'TDMD', 'FODMD', Location='best')
title('Residual by iterations')
xlabel('Iteration')
ylabel('Residual')
xlim tight

figname = 'iterres';
saveas(gcf, [case_name, '_', figname, '.png']);

%% 5. mode - field



%% Reference based comparison ================================================================
%% 1. residual map - ref ver.
id = init;
Pdmd = P_dmd{id}; vrdmd = vr_dmd{id}; Ddmd = evals_dmd{id};
Ptd = P_td{id}; vrtd = vr_td{id}; Dtd = evals_td{id};
Pon = P_on{id}; vron = vr_on{id}; Don = evals_on{id};

figure()
hold on
pcolor(re_map, im_map, res_map)
colormap sky
hold on
shading interp
cb = colorbar(); cb.Label.String = 'Residual';
contour(re_map, im_map, res_map, linspace(min_res, max_res, 8), '-k', LineWidth=0.5)
scref = scatter(real(Dk), imag(Dk), 64, 'blue', 'o', 'filled', 'MarkerFaceAlpha', 1., 'MarkerEdgeAlpha', 1.);
scdmd = scatter(real(Ddmd), imag(Ddmd), 'square', 'MarkerFaceAlpha', 0., 'LineWidth', 1.);
sctd = scatter(real(Dtd), imag(Dtd), 'green', 'diamond', 'MarkerFaceAlpha', 0., 'LineWidth', 1.);
scon= scatter(real(Don), imag(Don), 'red', '+', 'LineWidth', 1.5);

uc = plot(cos(thetas), sin(thetas), 'k--', LineWidth=1.5);
legend([scref, scdmd, sctd, scon, uc], {'Real', 'DMD', 'TDMD', 'FODMD', 'Unit circle'}, Location='southwest')
title('Residual map and the spectrum results - Before')
xlabel('Re')
ylabel('Im')
axis equal
xlim([0.09, 1])
ylim([0., 0.6])

figname = 'mapbefore';
saveas(gcf, [case_name, '_', figname, '.png']);

% ================================

id = steps;
Pdmd = P_dmd{id}; vrdmd = vr_dmd{id}; Ddmd = evals_dmd{id};
Ptd = P_td{id}; vrtd = vr_td{id}; Dtd = evals_td{id};
Pon = P_on{id}; vron = vr_on{id}; Don = evals_on{id};

figure()
hold on
pcolor(re_map, im_map, res_map)
colormap sky
hold on
shading interp
cb = colorbar(); cb.Label.String = 'Residual';
contour(re_map, im_map, res_map, linspace(min_res, max_res, 8), '-k', LineWidth=0.5)
scref = scatter(real(Dk), imag(Dk), 64, 'blue', 'o', 'filled', 'MarkerFaceAlpha', 1., 'MarkerEdgeAlpha', 1.);
scdmd = scatter(real(Ddmd), imag(Ddmd), 'square', 'MarkerFaceAlpha', 0., 'LineWidth', 1.);
sctd = scatter(real(Dtd), imag(Dtd), 'green', 'diamond', 'MarkerFaceAlpha', 0., 'LineWidth', 1.);
scon= scatter(real(Don), imag(Don), 'red', '+', 'LineWidth', 1.5);

uc = plot(cos(thetas), sin(thetas), 'k--', LineWidth=1.5);
legend([scref, scdmd, sctd, scon, uc], {'Real', 'DMD', 'TDMD', 'FODMD', 'Unit circle'}, Location='southwest')
title('Residual map and the spectrum results - After')
xlabel('Re')
ylabel('Im')
axis equal
xlim([0.09, 1])
ylim([0., 0.6])

figname = 'mapafter';
saveas(gcf, [case_name, '_', figname, '.png']);

%% 2. phase - modulus - ref ver.
id = steps;
Pdmd = P_dmd{id}; vrdmd = vr_dmd{id}; Ddmd = evals_dmd{id};
Ptd = P_td{id}; vrtd = vr_td{id}; Dtd = evals_td{id};
Pon = P_on{id}; vron = vr_on{id}; Don = evals_on{id};
figure()
hold on
scref = scatter(angle(Dk), abs(Dk), 64, 'blue', 'o', 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
scdmd = scatter(angle(Ddmd), abs(Ddmd), 'square', 'MarkerFaceAlpha', 0., 'LineWidth', 1.);
sctd = scatter(angle(Dtd), abs(Dtd), 'green', 'diamond', 'MarkerFaceAlpha', 0., 'LineWidth', 1.);
scon= scatter(angle(Don), abs(Don), 'red', '+', 'LineWidth', 1.5);
legend([scref, scdmd, sctd, scon], {'Real', 'DMD', 'TDMD', 'FODMD'}, Location='best')
title('Phase-modulus diagram of the final spectrums')
xlabel('Phase')
ylabel('Modulus')
yscale log
axis padded
box on

%% 3. residual - real spec.
figure()
hold on
pcolor(re_map, im_map, res_map)
colormap sky
hold on
shading interp
cb = colorbar(); cb.Label.String = 'Residual';
contour(re_map, im_map, res_map, linspace(min_res, max_res, 8), '-k', LineWidth=0.5)
scref = scatter(real(Dall), imag(Dall), 'red', 'o', 'filled', 'MarkerFaceAlpha', 1., 'MarkerEdgeAlpha', 1.);
uc = plot(cos(thetas), sin(thetas), 'k--', LineWidth=1.5);
legend([scref, uc], {'Real spectrums', 'Unit circle'}, Location='southwest')
title('Residual map and the real spectrums')
xlabel('Re')
ylabel('Im')
axis equal

figname = 'mapreal';
saveas(gcf, [case_name, '_', figname, '.png']);

%% 4. span dist
timeline = init: steps;
for i = timeline
    dist_dmd(i) = span_distance(P_dmd{i}, Vk);
    dist_td(i) = span_distance(P_td{i}, Vk);
    dist_on(i) = span_distance(P_on{i}, Vk);
end
figure()
hold on
plot(timeline-init, dist_dmd(timeline), '-*')
plot(timeline-init, dist_td(timeline), '-*')
plot(timeline-init, dist_on(timeline), '-o', LineWidth=1.5)
legend('DMD', 'TDMD', 'Online-debiased', Location='best')
title('Eigenspace error by iterations')
xlabel('Iteration')
ylabel('Residual')
xlim tight

figname = 'iterdist';
saveas(gcf, [case_name, '_', figname, '.png']);

% 5. subspace projection err (same as 4.)
timeline = init: steps;
for i = timeline
    dist_dmd(i) = subspace_distance(Vk, P_dmd{i});
    dist_td(i) = subspace_distance(Vk, P_td{i});
    dist_on(i) = subspace_distance(Vk, P_on{i});
end
figure()
hold on
plot(timeline-init, dist_dmd(timeline), '-*')
plot(timeline-init, dist_td(timeline), '-*')
plot(timeline-init, dist_on(timeline), '-o', LineWidth=1.5)
legend('DMD', 'TDMD', 'Online-debiased', Location='best')
title('Subspace projection error of the real eigenmatrix by iterations')
xlabel('Iteration')
ylabel('Residual')
xlim tight





