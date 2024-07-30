function recon = dmd_recon(mode, eigenvalue, zero_state, step_min, step_max)
%%% Predictor function after DMD, push forward system prediction from a
%%% given zero-state
%%% input: [modes, eigenvalues, zero-state, step_min, step_max]
%%% output: reconstruction

ampl = diag(pinv(mode) * zero_state);      % r*r
evol = eigenvalue .^ (step_min: step_max);        % a trick to define vandermonde matrix
recon = real(mode * ampl * evol);
end


