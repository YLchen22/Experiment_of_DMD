function [X, Y, data] = circle_norm_data_generator(dim, step, period)
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
    mux = sin((t-1) / period * 2 * pi); muy = cos((t-1) / period * 2 * pi);
    mu = [mux, muy];
    snapshot = mvnpdf([X(:), Y(:)], mu);
    snapshot = snapshot * 5 + 1;
    data(:, t) = snapshot;
end

end

