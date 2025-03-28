function [res, best_vec] = mini_residual(lambda, featx, featy)
%MINI_RESIDUAL: Solve a best mode by minimizing the residual
%   lambda: pseudo-eigenvalue
%   featx, featy: data projected to the feature space
%   RETURN:
%   res, best_vec: corresponding residual and the best mode

%   weight matrix W is not included (considered as identity)

    [Q,R] = qr(featx,"econ");
    RY = featy / R;
    RYTY = RY' * RY;
    RXTY = Q' * RY;
    C = RYTY - lambda * RXTY' - lambda' * RXTY + abs(lambda)^2 * eye(size(RYTY));

    [rvec, res2] = eigs(C, 1, 'smallestabs');
    best_vec = featx * rvec;
    res = sqrt(abs(res2));

end

