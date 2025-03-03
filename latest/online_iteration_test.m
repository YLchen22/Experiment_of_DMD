function [record_evals, record_vr, record_P, record_B] = online_iteration_test(data, init, r, w)

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
    w = 100;
end

k = 0; r = r+k;
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
        B = v / s;
    end

    % after the first, assume P = XB here
    Ar = P' * Y * B;    % should equal to P' * Y * pinv(X) * P
    [vr, evals] = main_eig(Ar, r);

    %% record current iteration
    record_evals{i} = evals;
    record_vr{i} = vr;
    record_P{i} = P;
    record_B{i} = B;    % P{t} == Dt(:, 1:t-1) * B{t}
    
    %% update!
    BU_new = [B; zeros(1, r)];

    for j = 1:r
        % these ar verified
        temp1 = B * vr(:, j) / evals(j);
        term1 = [0; temp1];
        
        lambda_minus = diag(evals(j) - evals);
        temp2 = B * vr * pinv(lambda_minus) / vr * evals(j) * vr(:, j);
        term2 = [temp2; 0];

        BU_new(:, j) = term1 + term2;
    end

    vecs_new = Dt * BU_new;
    % vecs_new = vecs_new + eye(size(vecs_new)) * 1e-4;     
    % regularization is not necessary, but seems better if added?
    [P_new, R] = qr(vecs_new, 'econ');
    
    %% this seems to be a mystery but useful trick...
    %% 前面额外多乘一项 P_new' * P_new 是用来消除误差的，我也不知道误差从哪来，而且为什么乘完之后就正常了。。。

    % redun = 0;
    % Dt_proj = P_new' * P_new * P_new' * Dt(:, end-r-redun+1: end); % use r+redun columns
    % B_proj = pinv(Dt_proj);   % pinv(r * r+redun)
    % B_new = [zeros(i-r-redun, r); B_proj];

    window = min([i, w]);
    Dt_proj = P_new' * P_new * P_new' * Dt(:, end-window+1: end); % use w columns
    B_proj = pinv(Dt_proj);     % pinv(r * i)
    B_new = [zeros(i-window, r); B_proj];

    %% pass to the next iter: 
    P = P_new;
    B = B_new;
end

end

