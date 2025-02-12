clear
rng(2024);

% matrix size
n = 20;
% randomly generate eigenvlaues, eigenvectors recover matrix
[A, evals, evecs] = rand_mat_real(n);
[u, s, v] = svd(A);

% initial state: a certain combination of the eigen-vectors
x0 = evecs * (1:n)';

% lower rank for dmd, and the steps for the simulation snapshots
r = 10; steps = 11;

% simulate snapshots observation
d = zeros([n, steps]);
for i = 1:steps
    d(:, i) = A^(i-1) * x0;
end
X = d(:, 1:end-1); Y = d(:, 2:end);

A_pred = Y * pinv(X);
[pvecs, pvals] = eig(A_pred);
[pvecs, pvals] = sort_by_val(pvecs, pvals);


function [Vs, Ds] = sort_by_val(V, D)
    %%% simply sort the [eigenvectors, eigenvalues] by the module of
    %%% eigenvalues
    eigenvalues = abs(diag(D));
    [~, indices] = sort(eigenvalues, 'descend');
    Vs = V(:, indices);
    D = diag(D);
    Ds = diag(D(indices));
end


% random matrix with complex eigenvalues and vectors
% (same as the test case of diff-svd)
function A = rand_mat(n)
    mat = normrnd(0, 1, [n, n]);
    [u, ~] = qr(mat);
    
    mat = normrnd(0, 1, [n, n]);
    [v, ~] = qr(mat);
    
    s = logspace(0, -1, n);
    
    A = u * diag(s) * v';
end


% random matrix with real eigenvalues and vectors
function [A, evals, evecs] = rand_mat_real(n)
    evals = logspace(0, -1, n);
    evecs = rand_vec(n, n);
    A = evecs * diag(evals) / evecs;
end


function prod = vec_similarity(a, b)
    a = a / norm(a);
    b = b / norm(b);
    prod = norm(a'*b) / norm(a) / norm(b);
end
