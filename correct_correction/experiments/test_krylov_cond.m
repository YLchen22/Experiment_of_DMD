clear
rng(2024);

n = 400;
dim = n;
[A, evals, evecs] = rand_mat_real(n);

% initial state: a certain combination of the eigen-vectors
x0 = evecs * (1:n)';

% lower rank for dmd, and the steps for the simulation snapshots
r = 10; steps = 100;
or_evals = evals(1:r); or_evecs = evecs(:, 1:r);

% simulate snapshots observation
d = zeros([n, steps]);
for i = 1:steps
    d(:, i) = A^(i-1) * x0;
end

for i = 1:steps
    dp = d(:, 1:i);
    c(i) = cond(dp);
end

figure()
plot(c)
yscale log



function [A, evals, evecs] = rand_mat_real(n)
    evals = logspace(0, -2, n) .* (1 + 0.1*randn(1, n));
    evecs = rand_col(n, n);
    A = evecs * diag(evals) / evecs;
end

function vecs = rand_col(dim, num)
% generate a matrix of column vecters, randomly genereted with directions
% input:
% dim: dimension of the vectors, row of vecs
% num: number of the random vectors. column of vecs
% output:
% vecs: dim*num matrix. each column is a random vector
vecs = unifrnd(0, 2*pi, [dim, num]);
vecs = cos(vecs);
len = vecnorm(vecs, 2, 1);
vecs = vecs ./ len;
end

