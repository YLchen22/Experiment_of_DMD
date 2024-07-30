function [mode, eigenvalue] = dmd_decom(X, Y, rank)
%%% one-step DMD function, return first r modes and eigenvalues
%%% input: [X_data, Y_data, lower rank]
%%% output: [modes, eigenvalues]

[U, S, V] = svds(X, rank);
tilde_A = U' * Y * V / S;   % is a r*r matrix, DMD

[tilde_Phi, tilde_Lambda] = eig(tilde_A); tilde_Lambda = diag(tilde_Lambda);

mode = U * tilde_Phi;
eigenvalue = tilde_Lambda;
end

