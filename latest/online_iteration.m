function [record_evals, record_vr, record_P, record_B] = online_iteration(data, init, r, option)

% input full data matrix, initial number, and low-rank number.
% perform CHEAP or EXPENSIVE version of online debiasing DMD (selectable)
% return: records of four main indicators as cells.
%   record_evals: the predicted eigenvalues.
%   record_vr: the predicted low-ranked eigenvectors.
%   record_P: the predicted subspace.
%   record_B: the predicted parameter matrix of subspace.
%   ***record{t} is the result when t-columns of data are used.

%   Specifically, we have (predicted) evecs = P * vr.
%   and expect: P{i} = X{i} * B{i}

if nargin < 4
    option = 'cheap';
end

steps = size(data, 2);

for i = init: steps

    Dt = data(:, 1: i);
    X = Dt(:, 1: end-1);
    Y = Dt(:, 2: end);

    %% for the first iter, use DMD to start
    if i == init
        [u, s, v] = svds(X, r);
        Ar = u' * Y * v / s;
        P = u;
        B = pinv(X) * P;
    else
        % after the first, assume P = XB here
        Ar = P' * Y * B;    % should equal to P' * Y * pinv(X) * P
    end

    [vr, evals] = main_eig(Ar, r);
    evecs = P * vr;

    %% record current iteration
    record_evals{i} = evals;
    record_vr{i} = vr;
    record_P{i} = P;
    record_B{i} = B;    % P{t} == Dt(:, 1:t-1) * B{t}
    
    %% update!
    B_new = [B; zeros(1, r)];

    for j = 1:r
        % these ar verified
        temp1 = B * vr(:, j) / evals(j);
        term1 = [0; temp1];
        
        lambda_minus = diag(evals(j) - evals);
        temp2 = B * vr * pinv(lambda_minus) / vr * P' * Y * B * vr(:, j);
        term2 = [temp2; 0];

        B_new(:, j) = term1 + term2;
    end

    vecs_new = Dt * B_new;
    [P_new, R] = qr(vecs_new, 'econ');

    if strcmp(option, 'cheap')
        B_new = B_new / R;
    else
        B_new = pinv(Dt) * P_new;   % theoretically, P_new = Dt * B
    end

    %% pass to the next iter: 
    P = P_new;
    B = B_new;
end

end

