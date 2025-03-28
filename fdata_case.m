clear all
close all
rng(2024)


%% select data
case_name = 'cylinder';
% case_name = 'dam';    % recommended
% case_name = 'tube';   % recommended
% case_name = 0;
tr = false;
noise_scale = 0;

switch case_name
    case 'artificial'
        n = 600;    % dimension
        steps = 201;    % data size
        %% generate matrix and data
        % [A_org, evals, evecs] = rand_mat_real(n);
        % [A_org, evals, evecs] = rand_mat_sym(n);
        % [A_org, evals, evecs] = case1(n);
        [A_org, evals, evecs] = case2(n);
        distri = ones(n, 1);
        x0 = evecs * distri;    % control the beginning components
        data = zeros([n, steps]);
        data(:, 1) = x0;
        for i = 2:steps
            data(:, i) = A_org * data(:, i-1);
        end
        % data = data + 1e-8 * eye(size(data));
        % data = data + 1e-6 * randn(size(data));     % Robust!

        nx = n; ny = 1;

    case 'cavity'
        reyn = 13;     % 13, or 16, 19, 20, 30
        load_time = 200;
        file_name = ['Cavity', num2str(reyn), 'k.mat'];
        load(file_name);
        nx = VelocityField.N+1; ny = nx;    tr = true;
        data = VelocityField.Psi(:, 1: load_time);

    case 'cylinder'
        file_name = 'CYLINDER_ALL.mat';
        load(file_name);
        data = VORTALL;

    case 'dam'
        file_name = 'dam_bc_case0000.mat';
        load(file_name);
        data = vort;
        nx = 64; ny = 64;   tr = true;

    case 'tube'
        file_name = 'tube_bc_case0000.mat';
        load(file_name);
        data = vort;
        nx = 64; ny = 64;   tr = true;

    case 'neuron'
        file_name = 'ecog_window.mat';
        load(file_name);
        data = X(:, 1:200);
        nx = 59; ny = 1;

    otherwise
        disp('================================================================')
        disp('Failed to find data!')
        disp('================================================================')
        return
end


disp(['Mean absolute value of clean data:', num2str(mean(abs(data), 'all'))])
noise_scale = 0.05 * mean(abs(data), 'all');
raw_data = data;
n = size(raw_data, 1);
steps = round(size(data, 2) / 3 * 2);
veri_steps = size(data, 2) - steps;

clean_data = raw_data(:, 1:steps);
data = clean_data + noise_scale * randn(n ,steps);
veri_data = raw_data(:, steps+1: steps+veri_steps);
X = data(:, 1:end-1); Y = data(:, 2:end);

%% algorithm hyper-parameters
% n = 400;    % dimension
% steps = 110;    % data size
% veri_steps = steps*2; % preserve some data for verification
r = 10;     % lower-rank
init = steps-r;  % initial step for the alg
w = init;

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
legend([scdmd, sctd, scon, uc], {'DMD', 'TDMD', 'Online-debiased', 'Unit circle'}, Location='southwest')
title('Residual map and the spectrum results - Before')
xlabel('Re')
ylabel('Im')
axis equal
pause(0.5)
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
% scref = scatter(real(Dk), imag(Dk), 64, 'blue', 'o', 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
scdmd = scatter(real(Ddmd), imag(Ddmd), 'square', 'MarkerFaceAlpha', 0., 'LineWidth', 1.);
sctd = scatter(real(Dtd), imag(Dtd), 'green', 'diamond', 'MarkerFaceAlpha', 0., 'LineWidth', 1.);
scon= scatter(real(Don), imag(Don), 'red', '+', 'LineWidth', 1.5);

uc = plot(cos(thetas), sin(thetas), 'k--', LineWidth=1.5);
legend([scdmd, sctd, scon, uc], {'DMD', 'TDMD', 'Online-debiased', 'Unit circle'}, Location='southwest')
title('Residual map and the spectrum results - After')
xlabel('Re')
ylabel('Im')
axis equal
pause(0.5)
figname = 'mapafter';
saveas(gcf, [case_name, '_', figname, '.png']);

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
legend('DMD', 'TDMD', 'Online-debiased', Location='best')
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
legend('DMD', 'TDMD', 'Online-debiased', Location='northwest')
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
legend('DMD', 'TDMD', 'Online-debiased', Location='best')
title('Residual by iterations')
xlabel('Iteration')
ylabel('Residual')
xlim tight

figname = 'iterres';
saveas(gcf, [case_name, '_', figname, '.png']);

%% 5. mode - field
id = steps;
Pdmd = P_dmd{id}; vrdmd = vr_dmd{id}; Ddmd = evals_dmd{id};
Ptd = P_td{id}; vrtd = vr_td{id}; Dtd = evals_td{id};
Pon = P_on{id}; vron = vr_on{id}; Don = evals_on{id};

Von = Pon*vron; Vtd = Ptd*vrtd;
V = Von; D = Don;
figure(Position=[100, 100, 1400, 700])
hold on
sgtitle('Top 4 debiased modes')
for j = 5:8
    i = j-4;
    piece = reshape(V(:, j), nx, ny);
    if tr
        piece = piece.';
    end
    subplot(4, 2, 2*i-1)
    hold on
    title(['Real part of mode ', '# ', num2str(i), ', Eigenvalue = ', num2str(D(i), '%.2f')])
    colormap hsv
    pcolor(real(piece));
    shading interp
    cb = colorbar(); cb.Label.String = 'Re';
    axis tight

    subplot(4, 2, 2*i)
    hold on
    title(['Imag part of mode ', '# ', num2str(i)])
    colormap hsv
    pcolor(imag(piece));
    shading interp
    cb = colorbar(); cb.Label.String = 'Im';
    axis tight
end

figname = '4db_modes';
saveas(gcf, [case_name, '_', figname, '.png']);

% ================================
V = Vtd; D = Dtd;
figure(Position=[100, 100, 1400, 700])
hold on
sgtitle('Top 4 DMD modes')
for i = 1:4
    piece = reshape(V(:, i), nx, ny);
    if tr
        piece = piece.';
    end
    subplot(4, 2, 2*i-1)
    hold on
    title(['Real part of mode ', '# ', num2str(i), ', Eigenvalue = ', num2str(D(i), '%.2f')])
    colormap hsv
    pcolor(real(piece));
    shading interp
    cb = colorbar(); cb.Label.String = 'Re';
    axis tight

    subplot(4, 2, 2*i)
    hold on
    title(['Imag part of mode ', '# ', num2str(i)])
    colormap hsv
    pcolor(imag(piece));
    shading interp
    cb = colorbar(); cb.Label.String = 'Im';
    axis tight
end

figname = '4td_modes';
saveas(gcf, [case_name, '_', figname, '.png']);

%% view
% figure(Position=[100, 100, 1500, 400])
% hold on
% for i = 1: veri_steps
%     subplot(131)
%     hold on
%     title('Data')
%     piece = reshape(veri_data(:, i), [nx, ny]);
%     if tr
%         piece = piece.';
%     end
%     colormap hsv
%     pcolor(real(piece));
%     shading interp
%     colorbar()
%     axis tight
% 
%     subplot(132)
%     hold on
%     title('DMD')
%     piece = reshape(pred_td(:, i), [nx, ny]);
%     if tr
%         piece = piece.';
%     end
%     colormap hsv
%     pcolor(real(piece));
%     shading interp
%     colorbar()
%     axis tight
% 
%     subplot(133)
%     hold on
%     title('DB')
%     piece = reshape(pred_on(:, i), [nx, ny]);
%     if tr
%         piece = piece.';
%     end
%     colormap hsv
%     pcolor(real(piece));
%     shading interp
%     colorbar()
%     axis tight
% 
%     pause(0.5)
% end
    



% %% Reference based comparison ================================================================
% %% 1. residual map - ref ver.
% id = init;
% Pdmd = P_dmd{id}; vrdmd = vr_dmd{id}; Ddmd = evals_dmd{id};
% Ptd = P_td{id}; vrtd = vr_td{id}; Dtd = evals_td{id};
% Pon = P_on{id}; vron = vr_on{id}; Don = evals_on{id};
% 
% figure()
% hold on
% pcolor(re_map, im_map, res_map)
% colormap sky
% hold on
% shading interp
% cb = colorbar(); cb.Label.String = 'Residual';
% contour(re_map, im_map, res_map, linspace(min_res, max_res, 8), '-k', LineWidth=0.5)
% scref = scatter(real(Dk), imag(Dk), 64, 'blue', 'o', 'filled', 'MarkerFaceAlpha', 1., 'MarkerEdgeAlpha', 1.);
% scdmd = scatter(real(Ddmd), imag(Ddmd), 'square', 'MarkerFaceAlpha', 0., 'LineWidth', 1.);
% sctd = scatter(real(Dtd), imag(Dtd), 'green', 'diamond', 'MarkerFaceAlpha', 0., 'LineWidth', 1.);
% scon= scatter(real(Don), imag(Don), 'red', '+', 'LineWidth', 1.5);
% 
% uc = plot(cos(thetas), sin(thetas), 'k--', LineWidth=1.5);
% legend([scref, scdmd, sctd, scon, uc], {'Real', 'DMD', 'TDMD', 'Online-debiased', 'Unit circle'}, Location='southwest')
% title('Residual map and the spectrum results - Before')
% xlabel('Re')
% ylabel('Im')
% axis equal
% xlim([0.12, 1])
% ylim([0., 0.88])
% 
% figname = 'mapbefore';
% saveas(gcf, [case_name, '_', figname, '.png']);
% 
% % ================================
% 
% id = steps;
% Pdmd = P_dmd{id}; vrdmd = vr_dmd{id}; Ddmd = evals_dmd{id};
% Ptd = P_td{id}; vrtd = vr_td{id}; Dtd = evals_td{id};
% Pon = P_on{id}; vron = vr_on{id}; Don = evals_on{id};
% 
% figure()
% hold on
% pcolor(re_map, im_map, res_map)
% colormap sky
% hold on
% shading interp
% cb = colorbar(); cb.Label.String = 'Residual';
% contour(re_map, im_map, res_map, linspace(min_res, max_res, 8), '-k', LineWidth=0.5)
% scref = scatter(real(Dk), imag(Dk), 64, 'blue', 'o', 'filled', 'MarkerFaceAlpha', 1., 'MarkerEdgeAlpha', 1.);
% scdmd = scatter(real(Ddmd), imag(Ddmd), 'square', 'MarkerFaceAlpha', 0., 'LineWidth', 1.);
% sctd = scatter(real(Dtd), imag(Dtd), 'green', 'diamond', 'MarkerFaceAlpha', 0., 'LineWidth', 1.);
% scon= scatter(real(Don), imag(Don), 'red', '+', 'LineWidth', 1.5);
% 
% uc = plot(cos(thetas), sin(thetas), 'k--', LineWidth=1.5);
% legend([scref, scdmd, sctd, scon, uc], {'Real', 'DMD', 'TDMD', 'Online-debiased', 'Unit circle'}, Location='southwest')
% title('Residual map and the spectrum results - After')
% xlabel('Re')
% ylabel('Im')
% axis equal
% xlim([0.12, 1])
% ylim([0., 0.88])
% 
% figname = 'mapafter';
% saveas(gcf, [case_name, '_', figname, '.png']);
% 
% %% 2. phase - modulus - ref ver.
% id = steps;
% Pdmd = P_dmd{id}; vrdmd = vr_dmd{id}; Ddmd = evals_dmd{id};
% Ptd = P_td{id}; vrtd = vr_td{id}; Dtd = evals_td{id};
% Pon = P_on{id}; vron = vr_on{id}; Don = evals_on{id};
% figure()
% hold on
% scref = scatter(angle(Dk), abs(Dk), 64, 'blue', 'o', 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
% scdmd = scatter(angle(Ddmd), abs(Ddmd), 'square', 'MarkerFaceAlpha', 0., 'LineWidth', 1.);
% sctd = scatter(angle(Dtd), abs(Dtd), 'green', 'diamond', 'MarkerFaceAlpha', 0., 'LineWidth', 1.);
% scon= scatter(angle(Don), abs(Don), 'red', '+', 'LineWidth', 1.5);
% legend([scref, scdmd, sctd, scon], {'Real', 'DMD', 'TDMD', 'Online-debiased'}, Location='best')
% title('Phase-modulus diagram of the final spectrums')
% xlabel('Phase')
% ylabel('Modulus')
% yscale log
% axis padded
% box on
% 
% %% 3. residual - real spec.
% figure()
% hold on
% pcolor(re_map, im_map, res_map)
% colormap sky
% hold on
% shading interp
% cb = colorbar(); cb.Label.String = 'Residual';
% contour(re_map, im_map, res_map, linspace(min_res, max_res, 8), '-k', LineWidth=0.5)
% scref = scatter(real(Dall), imag(Dall), 'red', 'o', 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
% uc = plot(cos(thetas), sin(thetas), 'k--', LineWidth=1.5);
% legend([scref, uc], {'Real spectrums', 'Unit circle'}, Location='southwest')
% title('Residual map and the real spectrums')
% xlabel('Re')
% ylabel('Im')
% axis equal
% 
% figname = 'mapreal';
% saveas(gcf, [case_name, '_', figname, '.png']);
% 
% %% 4. span dist
% timeline = init: steps;
% for i = timeline
%     dist_dmd(i) = span_distance(P_dmd{i}, Vk);
%     dist_td(i) = span_distance(P_td{i}, Vk);
%     dist_on(i) = span_distance(P_on{i}, Vk);
% end
% figure()
% hold on
% plot(timeline-init, dist_dmd(timeline), '-*')
% plot(timeline-init, dist_td(timeline), '-*')
% plot(timeline-init, dist_on(timeline), '-o', LineWidth=1.5)
% legend('DMD', 'TDMD', 'Online-debiased', Location='best')
% title('Space-distance by iterations')
% xlabel('Iteration')
% ylabel('Residual')
% xlim tight
% 
% figname = 'iterdist';
% saveas(gcf, [case_name, '_', figname, '.png']);
% 
% %% 5. subspace projection err (same as 4.)
% % timeline = init: steps;
% % for i = timeline
% %     dist_dmd(i) = span_distance(Vk, P_dmd{i});
% %     dist_td(i) = span_distance(Vk, P_td{i});
% %     dist_on(i) = span_distance(Vk, P_on{i});
% % end
% % figure()
% % hold on
% % plot(timeline-init, dist_dmd(timeline), '-*')
% % plot(timeline-init, dist_td(timeline), '-*')
% % plot(timeline-init, dist_on(timeline), '-o', LineWidth=1.5)
% % legend('DMD', 'TDMD', 'Online-debiased', Location='best')
% % title('Subspace projection error of the real eigenmatrix by iterations')
% % xlabel('Iteration')
% % ylabel('Residual')
% % xlim tight





