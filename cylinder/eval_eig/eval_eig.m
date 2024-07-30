% compare the system reconstruction of the original DMD and the new idea
% that reconstruct the system, and reconstruct the error of the
% reconstructed system, and add them together as debaising
clc
clear
close all
load cylinder.mat;  %%% data, nx, ny
%%% data = 89351 * 151 matrix, we use last 51 time steps as validation set

% DMD
r = 5; 
X = data(:, 1:30); Y = data(:, 2:31);
step_min = 0; step_max = 150;
tic;
[mode, eigenvalue] = dmd_decom(X, Y, r);
time_dmd = toc;

zero_state = data(:, 1);
recon = dmd_recon(mode, eigenvalue, zero_state, step_min, step_max);
recon1 = recon;
err1 = data - recon1;

%%% debias with the direct pinv(X) being allowed
mode_b = mode; eigenvalue_b = eigenvalue;
%%% necessary preparation
tic;
mode_inv = pinv(mode); lambda = diag(eigenvalue); pinv_x = pinv(X);
for i = 1: r    
    %%% select one mode to do
    mode_i = mode(:, i);
    %%% this term eqls to (\bar A - \hat A)*v_i
    diff_A_vi = Y * (pinv_x*mode_i) - lambda(i) * mode_i;

    mode_b(:, i) = mode * pinv(eigenvalue(i) * eye(r) - lambda)...
        * (mode_inv * diff_A_vi);
    eigenvalue_b(i) = mode_i.' * diff_A_vi;
end
mode_new = mode + mode_b; eigenvalue_new = eigenvalue + eigenvalue_b;
time_full = toc;

recon = dmd_recon(mode_new, eigenvalue, zero_state, step_min, step_max);
recon2 = recon;
err2 = data - recon2;

%%% debias with only the truncated pinv(X)
mode_b_tr = mode; eigenvalue_b_tr = eigenvalue;
%%% necessary preparation
tic;
[U, S, V] = svds(X, r);
mode_inv = pinv(mode); lambda = diag(eigenvalue); pinv_x = V / S * U';
for i = 1: r    
    %%% select one mode to do
    mode_i = mode(:, i);
    %%% this term eqls to (\bar A - \hat A)*v_i
    diff_A_vi = Y * (pinv_x*mode_i) - lambda(i) * mode_i;

    mode_b_tr(:, i) = mode * pinv(eigenvalue(i) * eye(r) - lambda)...
        * (mode_inv * diff_A_vi);
    eigenvalue_b_tr(i) = mode_i.' * diff_A_vi;
end
mode_new_tr = mode + mode_b_tr; eigenvalue_new_tr = eigenvalue + eigenvalue_b_tr;
time_trun = toc;

recon = dmd_recon(mode_new_tr, eigenvalue, zero_state, step_min, step_max);
recon3 = recon;
err3 = data - recon3;

disp('standard dmd run time:')
disp(time_dmd)
disp('econ svd run time:')
disp(time_full)
disp('truncated svd run time:')
disp(time_trun)

%%% visualize the MSE of 2 methods
MSE_DMD = mean(err1.^2);
MSE_full_db = mean(err2.^2);
MSE_trun_db = mean(err3.^2);
figure()
hold on
plot(MSE_DMD, 'DisplayName', 'MSE of the original DMD prediction')
plot(MSE_full_db, 'DisplayName', 'MSE of the debiased one with econ svd')
% plot(MSE_trun_db, 'DisplayName', 'MSE of the debiased one with truncated svd')
xline(101, '--', 'DisplayName', 'Time step mark')
xlabel('Time step')
ylabel('MSE of reconstruction')
legend()


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















