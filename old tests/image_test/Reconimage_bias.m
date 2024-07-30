%%% reconstruct image using data
clear
close all

%%% display original image
im = imread('image.jpg');
figure()
imshow(im); title('original image')

grayim = rgb2gray(im);
grayim = double(grayim) / 255;
figure()
imshow(grayim); title('original image (grayscale)')

%%% add noise
std_noise = 0.1;
A = grayim + randn(size(grayim)) * std_noise;
figure()
imshow(A); title('noised image')
%%% observe data
rng(2024)

n_obs = 500; r = 200;
X = randn([length(A), n_obs]); Y = A * X;
[U, S, V] = svd(X, 'econ');
U1 = U(:,1:r); S1 = S(1:r,1:r); V1 = V(:,1:r);
A1 = Y * V / S * U';
A2 = (U * U') * A1 * (U * U');
figure()
plot(diag(S)); title('singular values of X')
figure()
subplot(121)
imshow(A1); title('original')
subplot(122)
imshow(A2); title('hat matrix')
sgtitle([num2str(n_obs), ' observations and ', num2str(r), ' rank'])

[phi_real, lambda_real] = eig(A); lambda_real = diag(lambda_real);
re_real = real(lambda_real); im_real = imag(lambda_real);

[phi_1, lambda_1] = eig(A1); lambda_1 = diag(lambda_1);
re_1 = real(lambda_1); im_1 = imag(lambda_1);

[phi_hat, lambda_hat] = eig(A2); lambda_hat = diag(lambda_hat);
re_hat = real(lambda_hat); im_hat = imag(lambda_hat);

figure()
hold on
scatter(re_real, im_real, 'o', 'DisplayName', 'real results')
scatter(re_1, im_1, 'x', 'LineWidth', 2, 'DisplayName', 'reconstruction results')
scatter(re_hat, im_hat, 'x', 'LineWidth', 2, 'DisplayName', 'A hat results')
theta = linspace(0, 2*pi, 100); x = cos(theta); y = sin(theta);
plot(x, y, 'Color', 'black', 'HandleVisibility', 'off')
sgtitle('Eigenvalues'); legend();

%%% this part (calculating bias term) runs a long time
phi_db = phi_hat; lambda_db = lambda_hat;
for i = 1 : length(lambda_hat)
    new_Phi = phi_hat(:, i);
    lambda_term(i) = new_Phi' * (A1 - A2) * new_Phi;
    lambda_db(i) = lambda_hat(i) + lambda_term(i);
    phi_term(:, i) = pinv(lambda_hat(i) * eye(size(A)) - A2) * (A1 - A2) * new_Phi;
    phi_db(:, i) = new_Phi + phi_term(i);
end
re_db = real(lambda_db); im_db = imag(lambda_db);

figure()
hold on
scatter(re_real, im_real, 'o', 'DisplayName', 'real results')
scatter(re_1, im_1, 'x', 'LineWidth', 2, 'DisplayName', 'reconstruction results')
scatter(re_db, im_db, 'x', 'LineWidth', 2, 'DisplayName', 'debiased A hat results')
theta = linspace(0, 2*pi, 100); x = cos(theta); y = sin(theta);
plot(x, y, 'Color', 'black', 'HandleVisibility', 'off')
sgtitle('Eigenvalues'); legend();

err_hat = lambda_1 - lambda_hat;
err_db = lambda_1 - lambda_db;

figure()
hold on
plot(abs(err_hat))
plot(abs(err_db))
legend('abs(error hat)', 'abs(error db)')

figure()
hold on
plot(abs(lambda_term))
legend('abs(debias term)')

matdisp({A1-A2})

figure()
hold on
scatter(re_real(1:50), im_real(1:50), 'o', 'DisplayName', 'real results')
scatter(re_1(1:50), im_1(1:50), 'x', 'LineWidth', 2, 'DisplayName', 'reconstruction results')
scatter(re_hat(1:50), im_hat(1:50), 'x', 'LineWidth', 2, 'DisplayName', 'A hat results')
sgtitle('Eigenvalues'); legend();




