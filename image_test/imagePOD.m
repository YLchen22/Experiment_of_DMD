%%% directly used an image and do svd
clear all

%%% display original image
im = imread('image.jpg');
figure()
imshow(im); title('original image')

grayim = rgb2gray(im);
grayim = double(grayim) / 255;
figure()
imshow(grayim); title('original image (grayscale)')

%%% POD image
A = grayim; n = size(A, 2);
A_mean = mean(A, 2);
A = A - A_mean;
cov = A*A' / n;
[eigvec, eigval] = eig(cov);
eigval = diag(eigval); [~,idx] = sort(eigval, 'descend');
eigval = eigval(idx); eigvec = eigvec(:, idx);

r=10;
vec = eigvec(:, 1:r);
new_im = vec * vec' * A + A_mean;
figure()
imshow(new_im); title([num2str(r) 'rank reconstruction'])

r=100;
vec = eigvec(:, 1:r);
new_im = vec * vec' * A + A_mean;
figure()
imshow(new_im); title([num2str(r) 'rank reconstruction'])

r=200;
vec = eigvec(:, 1:r);
new_im = vec * vec' * A + A_mean;
figure()
imshow(new_im); title([num2str(r) 'rank reconstruction'])

r=300;
vec = eigvec(:, 1:r);
new_im = vec * vec' * A + A_mean;
figure()
imshow(new_im); title([num2str(r) 'rank reconstruction'])













