function data = sin_data_generator(dim, step, period)

% DATA_GENERATOR
% this func generates periodic sin snapshots
% u(x, t) = sin(x + (t-1)*2pi / p)

% field of x, 0: 2pi
field_x = linspace(0, 2*pi, dim)';

% obtain all x,0
field_xt = field_x * ones(1, step);

% obtain all x,t
t = 1: step;
field_xt = field_xt + (t-1)*2*pi / period;

% u = sin(x,...t)
data = sin(field_xt);

end

