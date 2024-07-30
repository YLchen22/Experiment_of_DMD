function [recon, mode, ampl, evol, tilde_Phi, tilde_Lambda] = ...
    dmd_predict(X, Y, rank, step_min, step_max)
%%% one-step DMD predictor function, push forward system prediction from
%%% the last snapshot given
%%% input: [X_data, Y_data, lower rank, step_min, step_max]
%%% output: [reconstruction, mode, amplitute, evolution, til_Phi, til_Lambda]

[U, S, V] = svds(X, rank);
tilde_A = U' * Y * V / S;   % is a r*r matrix, DMD

[tilde_Phi, tilde_Lambda] = eig(tilde_A); tilde_Lambda = diag(tilde_Lambda);

mode = U * tilde_Phi;
%%% DMD system prediction
ampl = diag(pinv(mode) * Y(:, end));      % r*r
evol = tilde_Lambda .^ (step_min: step_max);        % a trick to define vandermonde matrix
recon = real(mode * ampl * evol);
end

