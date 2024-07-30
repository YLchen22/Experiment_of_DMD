%%% See the propagation of error with different data amount
clc
clear
close all
load cylinder.mat;  %%% data, nx, ny
X = data; nt = size(X, 2);  %%% row=space=nx*ny, column=time

% DMD
X = data(:, 1:100); Y = data(:, 2:101);     %%% here use 100+1
r = 5;
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
recon1 = recon_dmd; bias1 = bias;

% DMD with less data
%%% here use 30+1. 30 is a boundary, because this is almost the time of a
%%% cycle. if you use data less than 30, prediction will be seriously wrong
X = data(:, 1:30); Y = data(:, 2:31);
r = 5;
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
recon2 = recon_dmd; bias2 = bias;

fig = figure('Position', [100, 100, 1200, 600], 'WindowStyle', 'normal');

%%% visualization step number 
steps = nt;
file_name = 'CylinderDMD_errprop.gif';

for i = 1:steps
    figure(fig)
    clf
    subplot(221)
    piece = reshape(recon1(:, i), nx, ny);
    plotCylinder_m(piece)
    title(['DMD reconstruction, step=', num2str(i)])

    subplot(222)
    piece = reshape(recon2(:, i), nx, ny);
    plotCylinder_m(piece)
    title('Less data DMD reconstruction')

    subplot(223)
    piece = reshape(bias1(:, i), nx, ny);
    plotCylinder_m(piece)
    title('Data - DMD')
    
    subplot(224)
    piece = reshape(bias1(:, i), nx, ny);
    plotCylinder_m(piece)
    title('Data - LessDMD')
    
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














