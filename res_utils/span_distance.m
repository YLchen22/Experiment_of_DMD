function delta = span_distance(vecs1, vecs2)
%% given two set of eigenvectors and compare their span distance
    
    r = size(vecs1, 2);
    [Q, ~] = qr(vecs1);
    Q1 = Q(:, 1:r); Qc1 = Q(:, r+1:end);

    r = size(vecs2, 2);
    [Q, ~] = qr(vecs2);
    Q2 = Q(:, 1:r); Qc2 = Q(:, r+1:end);

    % double way distance
    % d1 = norm(Qc1' * Q2);
    % d2 = norm(Q1' * Qc2);
    d1 = norm(Qc1' * Q2, 'fro');
    d2 = norm(Q1' * Qc2, 'fro');
    
    delta = max([d1, d2]);

end

