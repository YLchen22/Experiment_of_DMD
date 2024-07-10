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

%%% SVD image
[U, S, V] = svd(grayim, 'econ');
singular = diag(S);

figure()
plot(singular)

r=10;
U1 = U(:,1:r); S1 = S(1:r,1:r); V1 = V(:,1:r);
new_im = U1*S1*V1';
figure()
imshow(new_im); title([num2str(r) 'rank reconstruction'])

r=100;
U2 = U(:,1:r); S2 = S(1:r,1:r); V2 = V(:,1:r);
new_im = U2*S2*V2';
figure()
imshow(new_im); title([num2str(r) 'rank reconstruction'])

r=200;
U5 = U(:,1:r); S5 = S(1:r,1:r); V5 = V(:,1:r);
new_im = U5*S5*V5';
figure()
imshow(new_im); title([num2str(r) 'rank reconstruction'])














