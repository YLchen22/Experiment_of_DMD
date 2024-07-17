% compare the system reconstruction of the original DMD and the new idea
% that reconstruct the system, and reconstruct the error of the
% reconstructed system, and add them together as debaising
clc
clear
close all
load cylinder.mat;  %%% data, nx, ny
%%% data = 89351 * 151 matrix, we use last 51 time steps as validation set
valid_data = data(:, 101:151);
r = 10;

% DMD
X = data(:, 1:99); Y = data(:, 2:100);
step_min = 1; step_max = 51;
[recon, ~] = dmd_predict(X, Y, r, step_min, step_max);

recon1 = recon;
bias1 = valid_data - recon1;

% DMD with half data
X = data(:, 1:49); Y = data(:, 2:50);
step_min = 1; step_max = 101;
[recon, ~] = dmd_predict(X, Y, r, step_min, step_max);

recon_half = recon(:, 1:50);   % corresponding to time 51:100
recon2 = recon(:, 51:101);     % corresponding to time 101:151
bias2 = valid_data - recon2;

% DMD to approximate error propagation
error = data(:, 51:100) - recon_half;
X = error(:, 1:49); Y = error(:, 2:50);
step_min = 1; step_max = 51;
[recon, ~] = dmd_predict(X, Y, r, step_min, step_max);

recon3 = recon;
bias3 = (valid_data - recon2) - recon3;

DMD_prediction = recon1;
DMD_prediction_error = bias1;

DMD_debiased = recon2 + recon3;
DMD_half_error = valid_data - recon2;
DMD_debiased_error = valid_data - DMD_debiased;

target_bias = bias2;
approx_bias = recon3;

%%% visualize the MSE of 2 methods
MSE_DMD = mean(DMD_prediction_error.^2);
MSE_DMDhalf = mean(DMD_half_error.^2);
MSE_debiased = mean(DMD_debiased_error.^2);
figure()
hold on
plot(MSE_DMD)
plot(MSE_DMDhalf)
plot(MSE_debiased)
legend('MSE of the original DMD prediction',...
    'MSE of the original DMD prediction with half data',...
    'MSE of the debiased DMD prediction')

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















