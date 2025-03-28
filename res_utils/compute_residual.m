function res = compute_residual(lambda, feat_vec, featx, featy)
%MODE_RESIDUAL calculate the mode residual
%   lambda, feat_vec: eigenvalue and eigenvector of the feature
%   featx, featy: data projected to the feature space
%   RETURN
%   res: corresponding residual

%   weight matrix W is not included (considered as identity)
    featy = featy'; featx = featx';

    A = featy'*featy - lambda*(featx'*featy)' - lambda'*(featx'*featy) + abs(lambda)^2 * (featx'*featx);
    B = featx'*featx;

    res = (feat_vec' * A * feat_vec) / (feat_vec' * B * feat_vec);
    res = sqrt(real(res));
end

