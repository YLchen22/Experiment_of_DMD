function view_pred(file_name, t, from_beginning, gif_name)
%VIEW_PRED : use solved modes to predict system dynamics from the first
%snapshot or the last snapshot of data
%   file_name: a saved .mat file in the main function
%   t: steps wanted to process
%   from_beginning: true to start from the first snapshot, false the last
%   gif_name: If you want to save a gif, input the wanted gif name

    if nargin < 3
        from_beginning = true;
    end

    if nargin < 4
        gif_name = false;
    end

    % includes 'evals', 'vr', 'P', 'init', 'steps', 'nx', 'ny', 'tr', 'x0', 'xn'
    load(file_name)
    if from_beginning
        x = x0;
    else
        x = xn;
    end

    ampl = diag(vr{end} \ P{end}' * x);      % r*r
    evol = evals{end} .^ (0: t-1);            % a trick to define vandermonde matrix
    prediction = real(P{end} * vr{end} * ampl * evol);

    view_data(prediction, nx, ny, tr, gif_name)

end






