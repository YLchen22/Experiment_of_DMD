%%% reconstruct image using data
clear all

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
n_obs = 300;
X = randn([length(A), n_obs]); Y = A * X;




A1 = Y * pinv(X);
figure()
imshow(A1); title(['reconstructed matrix with ', num2str(n_obs), ' observations'])










