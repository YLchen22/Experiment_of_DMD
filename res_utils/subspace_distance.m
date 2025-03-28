function error = span_distance(ref_vec, Q)
%% given a set of real eigenvectors and a subspace (orthogonal matrix)
% compare the projection error of ref_vec to Q space

    error = norm(ref_vec - Q * Q' * ref_vec, 'fro');
end

