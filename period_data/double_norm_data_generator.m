function [X, Y, data] = double_norm_data_generator(dim, step, period)
% DATA_GENERATOR
% this func generates periodic double gaussian snapshots
% x is [-1, 1]*[-1, 1] grid in vector
% u(x, t) = gaussian_1(x,t) + gaussian_2(x,t)
% dim must be squared number

unit = dim ^ 0.5;
% field of x
x = linspace(-1, 1, unit);
y = linspace(-1, 1, unit);
[X, Y] = ndgrid(x, y);

data = zeros(dim, step);
for t = 1: step
    p1 = abs( t-1 - period/2 ) / period * 4 - 1; p2 = -p1;
    mu1 = [1, p1]; mu2 = [-1, p2];
    snapshot = mvnpdf([X(:), Y(:)], mu1) + mvnpdf([X(:), Y(:)], mu2);
    snapshot = snapshot * 10;
    data(:, t) = snapshot;
end

end

