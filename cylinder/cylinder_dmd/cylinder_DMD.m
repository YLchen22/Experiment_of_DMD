%% POD_Cylinder Wake
clc
clear
close all
load cylinder.mat;  %%% data, nx, ny
X = data; nt = size(X, 2);  %%% row=space=nx*ny, column=time
X = data(:, 1:end-1); Y = data(:, 2:end);

% DMD
%%% rank matters. 
%%% 10 is enough for acceptable result, but for example 5 is too small
r = 10;
[U, S, V] = svds(X, r);
tilde_A = U' * Y * V / S;   % is a r*r matrix, DMD

[tilde_Phi, tilde_Lambda] = eig(tilde_A); tilde_Lambda = diag(tilde_Lambda);
re_tilde = real(tilde_Lambda); im_tilde = imag(tilde_Lambda);

%%% try to use difference to approximately debias
mode = U * tilde_Phi;
%%% DMD system reconstruction
n = nt;     % reconstruction time steps
ampl = diag(pinv(mode) * X(:, 1));      % r*r
evol = tilde_Lambda .^ (0: n-1);        % a trick to define vandermonde matrix
recon_dmd = real(mode * ampl * evol);

bias = data - recon_dmd;

fig = figure('Position', [100, 100, 600, 800]);

%%% visualization step number 
steps = 50;
file_name = ['CylinderDMD_', num2str(steps), 'steps.gif'];

for i = 1:steps
    figure(fig)
    clf
    subplot(311)
    piece = reshape(data(:, i), nx, ny);
    plotCylinder_m(piece)
    title('Data')

    subplot(312)
    piece = reshape(X(:, i), nx, ny);
    plotCylinder_m(piece)
    title('DMD reconstruction')

    subplot(313)
    piece = reshape(bias(:, i), nx, ny);
    plotCylinder_m(piece)
    title('Data - DMD')
    
    pause(0.01)

    %%% save gif for visualization
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
end














