clear all
rng(2024)

dim = 400;
step = 50; % how many snapshots i want to generate
period = 200; % length (step number) of the period

[Xdata, Ydata, data] = circle_norm_data_generator(dim, step, period);
dmin = min(data, [], 'all'); dmax = max(data, [], 'all');

figure()
for t = 1:step

    state = data(:, t);
    Zdata = reshape(state, dim^0.5, dim^0.5);

    contourf(Xdata, Ydata, Zdata);
    colormap(jet);
    colorbar;
    clim([dmin, dmax]);

    pause(0.02)
end

X = data(:, 1:end-1); Y = data(:, 2:end);
A = Y * pinv(X);

norm(Y - A*X, 'fro')

% use first r-rank eigens to reconstruct the system
r = step;
[V, D] = main_eig(A, r);
A = V * diag(D) * pinv(V);
figure(1)
plot(abs(D))

approx_data = data;
state = data(:, 1);

figure(2)
for t = 1: 200
    approx_data(:, t) = state;

    Zdata = reshape(state, dim^0.5, dim^0.5);

    contourf(Xdata, Ydata, real(Zdata));
    colormap(jet);
    colorbar;
    clim([dmin, dmax]);

    pause(0.02)
    state = A * state;
end

% norm(data - approx_data, 'fro')



function [Vs, Ds] = sort_by_val(V, D)
    %%% simply sort the [eigenvectors, eigenvalues] by the module of
    %%% eigenvalues
    eigenvalues = abs(diag(D));
    [~, indices] = sort(eigenvalues, 'descend');
    Vs = V(:, indices);
    D = diag(D);
    Ds = diag(D(indices));
end


function [Vs, Ds] = main_eig(A, r)
    %%% simply sort the [eigenvectors, eigenvalues] by the module of
    %%% eigenvalues
    [V, D] = eig(A);
    [Vs, Ds] = sort_by_val(V, D);
    Vs = Vs(:, 1:r);
    Ds = diag(Ds(1:r, 1:r));
end


