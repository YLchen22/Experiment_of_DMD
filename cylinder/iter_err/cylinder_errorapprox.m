% compare the system reconstruction of the original DMD and the new idea
% that reconstruct the system, and reconstruct the error of the
% reconstructed system, and add them together as debaising
clc
clear
close all
load cylinder.mat;  %%% data, nx, ny
%%% data = 89351 * 151 matrix, we use last 51 time steps as validation set\
valid_data = data;

% ez DMD
r = 10;
n = 60;
X = data(:, 1:n); Y = data(:, 2:n+1);
[mode, eigenvalue] = dmd_decom(X, Y, r);

% ez recon DMD
step_min = 0; step_max = 150;
zero_state = data(:, 1);
recon = dmd_recon(mode, eigenvalue, zero_state, step_min, step_max);
recon1 = recon;
err1 = data - recon1;

% iter DMD and adding error term
% ez DMD uses 1:61 directly. here uses 1:31 to start and n columns per iter

per_col = 10;    %%% set this thing to decide added data per iter
iter = 30 / per_col + 1;

sum_recon = 0; err = data; %%% don't touch these two
for i = 1:iter
    n = 30 + per_col*(i-1);
    %%% this code is runned by magic, don't touch
    X = err(:, 1:n); Y = err(:, 2:n+1);
    zero_state = X(:, 1);

    [mode, eigenvalue] = dmd_decom(X, Y, r);
    step_min = 0; step_max = 150;
    recon = dmd_recon(mode, eigenvalue, zero_state, step_min, step_max);
    sum_recon = sum_recon + recon;
    modes{i} = mode; eigs{i} = eigenvalue;
    recons{i} = recon; sum_recons{i} = sum_recon;  %%% save number i-th err term

    err = data - sum_recon;
    errs{i} = err;
end

% recon with err terms being added
% step_min = 0; step_max = 150;
% zero_state = data(:, 1);
% sum_recons = cell(1, iter); errs = cell(1, iter); recon = 0;
% for i = 1:iter
%     recon = recon + recons{i}; err = data - recon;
%     sum_recons{i} = recon; errs{i} = err;
% end

%%% visualize the MSE of 2 methods
mse1 = mean(err1.^2);

figure()
hold on
plot(mse1, 'Color', 'black', 'LineWidth', 1.)
legendtext = {'MSE of DMD'};
%%% make lines and legends
for i = 1:iter
    mse2 = mean(errs{i}.^2);
    if i == iter
        plot(mse2, 'Color', 'blue', 'LineWidth', 1.)
    else
        plot(mse2)
    end
    str = ['MSE with iter-', num2str(i)];
    legendtext = [legendtext, str];
end
xlabel('Time step')
ylabel('MSE of reconstruction')
legend(legendtext); legend('Location', 'best')

figure()
hold on
legendtext = {};
%%% make lines and legends
for i = 2:iter
    mse2 = mean(errs{i}.^2);
    if i == iter
        plot(mse2, 'Color', 'blue', 'LineWidth', 1.)
    else
        plot(mse2)
    end
    str = ['MSE with iter-', num2str(i)];
    legendtext = [legendtext, str];
end
xlabel('Time step')
ylabel('MSE of reconstruction')
legend(legendtext); legend('Location', 'best')

%%% visualize the prediction and error
fig = figure('Position', [100, 100, 1500, 600]);
%%% visualization step number 
steps = 51;
file_name = ['Cylinder_errorapprox', num2str(r), 'rank.gif'];

figure(fig)
for i = 1:steps
    clf
    set(fig, 'WindowStyle', 'normal');
    subplot(231)
    piece = reshape(DMD_prediction(:, i), nx, ny);
    plotCylinder_m(piece)
    title(['Time = ', num2str(100+i), ', DMD prediction'])

    subplot(232)
    piece = reshape(DMD_debiased(:, i), nx, ny);
    plotCylinder_m(piece)
    title('Debiased DMD prediction')

    subplot(233)
    piece = reshape(target_bias(:, i), nx, ny);
    plotCylinder_m(piece)
    title('Target bias of the half data DMD')

    subplot(234)
    piece = reshape(DMD_prediction_error(:, i), nx, ny);
    plotCylinder_m(piece)
    title('Error of the DMD prediction')

    subplot(235)
    piece = reshape(DMD_debiased_error(:, i), nx, ny);
    plotCylinder_m(piece)
    title('Error of the debiased DMD prediction')

    subplot(236)
    piece = reshape(approx_bias(:, i), nx, ny);
    plotCylinder_m(piece)
    title('Approximate bias')
    
    pause(0.01)

    %%% save gif for visualization
    F=getframe(gcf);
    I=frame2im(F);
    [I,map]=rgb2ind(I,256);
    pic_num = i;
    if pic_num == 1
        imwrite(I,map,file_name,'gif','Loopcount',inf,'DelayTime',0.15);
    else
        imwrite(I,map,file_name,'gif','WriteMode','append','DelayTime',0.15);
    end
    pic_num = pic_num + 1;
end















