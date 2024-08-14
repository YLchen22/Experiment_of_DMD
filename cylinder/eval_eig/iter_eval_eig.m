% compare the system reconstruction of the original DMD and the new idea
% that reconstruct the system, and reconstruct the error of the
% reconstructed system, and add them together as debaising
clc
clear
close all
load cylinder.mat;  %%% data, nx, ny
%%% data = 89351 * 151 matrix, we use last 51 time steps as validation set

% DMD
r = 10;
n = 40;
X = data(:, 1:n); Y = data(:, 2:n+1);
step_min = 0; step_max = 150;

[mode, eigenvalue] = dmd_decom(X, Y, r);
mode1 = mode; eigenvalue1 = eigenvalue;

zero_state = data(:, 1);
recon = dmd_recon(mode, eigenvalue, zero_state, step_min, step_max);
recon1 = recon;
err1 = data - recon1;

%%% debias with the direct pinv(X) being allowed
n_start = n;
n = 60;
X = data(:, 1:n_start); Y = data(:, 2:n_start+1);
[mode, eigenvalue] = dmd_decom(X, Y, r);

mode_b = mode; eigenvalue_b = eigenvalue;

%%% necessary preparation
steps = n_start+1;
for i = steps:n+1
    X = data(:, i-n_start: i-1);
    Y = data(:, i-n_start+1: i);
    mode_inv = pinv(mode); lambda = diag(eigenvalue); pinv_x = pinv(X);
    for j = 1: r
        %%% select one mode to do
        mode_i = mode(:, j);
        %%% this term eqls to (\bar A - \hat A)*v_i
        diff_A_vi = Y * (pinv_x*mode_i) - lambda(j) * mode_i;

        mode_b(:, j) = mode * pinv(eigenvalue(j) * eye(r) - lambda)...
            * (mode_inv * diff_A_vi);
        eigenvalue_b(j) = mode_i.' * diff_A_vi;
    end
    mode_new = mode + mode_b; eigenvalue_new = eigenvalue + eigenvalue_b;
    mode = mode_new; eigenvlue = eigenvalue_new;
    i
end

recon = dmd_recon(mode, eigenvalue, zero_state, step_min, step_max);
recon2 = recon;
err2 = data - recon2;

%%% visualize the MSE of 2 methods
MSE_DMD = mean(err1.^2);
MSE_db = mean(err2.^2);
figure()
hold on
plot(MSE_DMD, 'DisplayName', 'DMD')
plot(MSE_db, 'DisplayName', 'DMD with modified eig')
% plot(MSE_trun_db, 'DisplayName', 'MSE of the debiased one with truncated svd')
xlabel('Time step')
ylabel('MSE of reconstruction')
legend('Location', 'best')

figure()
hold on
plot(MSE_DMD - MSE_db)
xlabel('Time step')
ylabel('Difference of MSE')

eig_diff = eigenvalue_new - eigenvalue1;
re = real(eig_diff); im = imag(eig_diff);
figure()
scatter(re, im)
xlabel('R(\Delta \lambda)')
ylabel('I(\Delta \lambda')

%%%
r = 5;
X = data(:, 1:n); Y = data(:, 2:n+1);
step_min = 0; step_max = 150;

[mode, eigenvalue] = dmd_decom(X, Y, r);
mode1 = mode; eigenvalue1 = eigenvalue;

zero_state = data(:, 1);
recon = dmd_recon(mode, eigenvalue, zero_state, step_min, step_max);
recon1 = recon;
err1 = data - recon1;

MSE_DMD = mean(err1.^2);
MSE_db = mean(err2.^2);
figure()
hold on
plot(MSE_DMD, 'DisplayName', 'DMD')
plot(MSE_db, 'DisplayName', 'DMD with modified eig')
% plot(MSE_trun_db, 'DisplayName', 'MSE of the debiased one with truncated svd')
xlabel('Time step')
ylabel('MSE of reconstruction')
legend('Location', 'best')


%%% visualize the prediction and error
% fig = figure('Position', [100, 100, 1500, 600]);
% %%% visualization step number 
% steps = 51;
% file_name = ['Cylinder_errorapprox', num2str(r), 'rank.gif'];
% 
% figure(fig)
% for i = 1:steps
%     clf
%     set(fig, 'WindowStyle', 'normal');
%     subplot(231)
%     piece = reshape(DMD_prediction(:, i), nx, ny);
%     plotCylinder_m(piece)
%     title(['Time = ', num2str(100+i), ', DMD prediction'])
% 
%     subplot(232)
%     piece = reshape(DMD_debiased(:, i), nx, ny);
%     plotCylinder_m(piece)
%     title('Debiased DMD prediction')
% 
%     subplot(233)
%     piece = reshape(target_bias(:, i), nx, ny);
%     plotCylinder_m(piece)
%     title('Target bias of the half data DMD')
% 
%     subplot(234)
%     piece = reshape(DMD_prediction_error(:, i), nx, ny);
%     plotCylinder_m(piece)
%     title('Error of the DMD prediction')
% 
%     subplot(235)
%     piece = reshape(DMD_debiased_error(:, i), nx, ny);
%     plotCylinder_m(piece)
%     title('Error of the debiased DMD prediction')
% 
%     subplot(236)
%     piece = reshape(approx_bias(:, i), nx, ny);
%     plotCylinder_m(piece)
%     title('Approximate bias')
%     
%     pause(0.01)
% 
%     %%% save gif for visualization
%     F=getframe(gcf);
%     I=frame2im(F);
%     [I,map]=rgb2ind(I,256);
%     pic_num = i;
%     if pic_num == 1
%         imwrite(I,map,file_name,'gif','Loopcount',inf,'DelayTime',0.15);
%     else
%         imwrite(I,map,file_name,'gif','WriteMode','append','DelayTime',0.15);
%     end
%     pic_num = pic_num + 1;
% end















