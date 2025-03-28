function pred = make_prediction(P, vr, evals, steps, x0)
%MAKE_PREDICTION: make prediction from a set of solved modes
% initial state is included. from x0 to x_t-1, t=steps in total.
%   INPUT:
%   P: orthogonal basis (corresponding to low-rank subspace)
%   vr: low-rank eigenvectors
%   evals: eigenvalues
%   stpes, x0: steps and initial state

    ampl = diag(vr \ P' * x0);      % r*r
    evol = evals .^ (0: steps-1);       % a trick to define vandermonde matrix
    pred = real(P * vr * ampl * evol);

end

