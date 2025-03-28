function res = sum_residual(lambda, feat_vec, featx, featy)
%MODE_RESIDUAL calculate the mode residual
%   lambda, feat_vec: MULTIPLE eigenvalues and eigenvectors of the feature
%   featx, featy: data projected to the feature space
%   RETURN
%   res: corresponding summed residual

%   weight matrix W is not included (considered as identity)

    r = length(lambda);
    res = 0;

    featy = featy'; featx = featx';

    for i = 1:r
        lam = lambda(i);
        g = feat_vec(:, i);
        
        Xg = featx*g; Yg = featy*g;
        A1 = Yg'*Yg;
        A2 = Yg'*Xg;
        A3 = Xg'*Yg;
        A4B = Xg'*Xg;
        
        A = A1 + lam*A2 + lam'*A3 + abs(lam)^2 * A4B;
        B = A4B;
        res_ = A/B;
        res = res + sqrt(real(res_));
    end
    
end

