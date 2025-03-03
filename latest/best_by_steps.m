function [record_evals, record_vr, record_P, record_B] = best_by_steps(data, init, r, Q2V, evals)

% input full data matrix, initial number, and low-rank number.
% perform DMD or Direct Regression to compare with the online iteration
% return: records of four main indicators as cells.
%   record_evals: the predicted eigenvalues.
%   record_evecs: the predicted eigenvectors.
%   record_P: the predicted subspace.

steps = size(data, 2);

for i = init: steps

    X = data(:, 1:i-1);
    Y = data(:, 2:i);

    

    record_evals{i} = evals;
    record_evecs{i} = evecs;
    record_P{i} = P;
end

end


function [V, D, Q] = simple_dmd(X, Y, r)
    [u, s, v] = svds(X, r);
    R = u' * Y * v * pinv(s);
    [V, D] = main_eig(R, r);
    Q = u;
    V = Q * V;
end
