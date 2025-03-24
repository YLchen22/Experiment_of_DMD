function view_mode(file_name, r)
%VIEW_MODE : load and visualize solved modes on a given field
%   file_name: a saved .mat file in the main function
%   r: the rank of mode wanted to view

    dt = 0.01;

    % includes 'evals', 'vr', 'P', 'init', 'steps', 'nx', 'ny', 'tr'
    load(file_name)   

    figure('Position', [100, 100, 1200, 600])
    % use a better color setting
    colormap hsv
    % load CCcool.mat
    % colormap(CC);
    % load batlow.mat
    % colormap(batlow)

    for i = init: steps
        evecs = P{i} * vr{i};
        all_evecs(:, i-init+1) = evecs(:, r);
    end
    rmax = max(real(all_evecs), [], 'all');
    rmin = min(real(all_evecs), [], 'all');
    imax = max(imag(all_evecs), [], 'all');
    imin = min(imag(all_evecs), [], 'all');
    
    for i = init: steps

        % view mode
        piece = reshape(all_evecs(:, i-init+1), [nx, ny]);

        if nx == 1
            piece = [piece; piece];
        elseif ny == 1
            piece = [piece, piece];
        end

        if tr
            piece = piece.';
        end
        
        subplot(121)
        pcolor(real(piece));
        clim([rmin, rmax]);
        shading interp
        colorbar()
        title('Real part')
        
        subplot(122)
        pcolor(imag(piece));
        clim([imin, imax]);
        shading interp
        colorbar()
        title('Imag part')

        sgtitle(['Mode #', num2str(r), ...
            ', with ', num2str(i), ' snapshots being used', ...
            '. ', 'Predicted eigenvalue = ', num2str(evals{i}(r))])

        pause(dt)
    end

end






