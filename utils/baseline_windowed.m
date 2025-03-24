function [record_evals, record_vr, record_P] = baseline_windowed(data, init, r, w, option)

% input full data matrix, initial number, and low-rank number.
% perform DMD or Direct Regression to compare with the online iteration
% return: records of four main indicators as cells.
%   record_evals: the predicted eigenvalues.
%   record_evecs: the predicted eigenvectors.
%   record_P: the predicted subspace.

if nargin < 5
    option = 'dmd';
end

steps = size(data, 2);

for i = init: steps
    window = min([i, w]);
    X = data(:, 1: i-1); Y = data(:, 2: i);
    D_window = data(:, i-window+1: i);

    if strcmp(option, 'dmd')
        [evecs, evals, P] = simple_dmd(X, Y, r);
        vr = P'*evecs;
    else
        [evecs, evals] = main_eig(Y*pinv(X), r);
        [P, ~] = qr(evecs, 'econ');
        [vr, evals] = main_eig(P' * Y * pinv(P'*X), r);
    end

    record_evals{i} = evals;
    record_vr{i} = vr;
    record_P{i} = P;
end

end


function [V, D, Q] = simple_dmd(X, Y, r)
    [u, s, v] = svds(X, r);
    % R = u' * Y * v * pinv(s);
    [V, D] = main_eig(u' * Y * pinv(u'*X), r);
    Q = u;
    V = Q * V;
end