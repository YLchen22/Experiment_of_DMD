function view_modes_by_rank(file_name, dt)
%VIEW_DATA : load and visualize solved modes on a given field
%   data_name: the used dataset

    if nargin < 2
        dt = 1.;
    end

    % includes 'evals', 'vr', 'P', 'init', 'steps', 'nx', 'ny', 'tr'
    load(file_name)   

    figure('Position', [100, 100, 1200, 600])
    % use a better color setting
    colormap hsv
    % load CCcool.mat
    % colormap(CC);
    % load batlow.mat
    % colormap(batlow)

    evecs = P{end} * vr{end};
    rank = size(evecs, 2);
    
    fig = figure(1);
    for r = 1:rank
        % view mode
        piece = reshape(evecs(:, r), [nx, ny]);

        if tr
            piece = piece.';
        end
        
        subplot(121)
        pcolor(real(piece));
        shading interp
        colorbar()
        title('Real part')
        
        subplot(122)
        pcolor(imag(piece));
        shading interp
        colorbar()
        title('Imag part')

        sgtitle(['Mode # ', num2str(r)])

        pause(dt)
    end

    % frame = getframe(fig);  % 获取当前图像帧
    % im = frame2im(frame); % 转换为图像
    % [A, map] = rgb2ind(im, 256); % 转换为索引图像
    % 
    % % 追加到 GIF
    % if i == init
    %     imwrite(A, map, [file_name, '_modes.gif'], 'gif', 'LoopCount', inf, 'DelayTime', 0.3);
    % else
    %     imwrite(A, map, [file_name, '_modes.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', 0.3);
    % end



end






