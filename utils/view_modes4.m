function view_modes4(file_name, gif_name)
%VIEW_MODES4 : load and visualize THE TOP 4 solved modes on a given field
%   file_name: a saved .mat file in the main function
%   gif_name: If you want to save a gif, input the wanted gif name

    if nargin < 2
        gif_name = false;
    end

    dt = 0.01;

    % includes 'evals', 'vr', 'P', 'init', 'steps', 'nx', 'ny', 'tr'
    load(file_name)   

    fig = figure('Position', [100, 100, 1200, 600]);
    % use a better color setting
    colormap hsv
    % load CCcool.mat
    % colormap(CC);
    % load batlow.mat
    % colormap(batlow)

    for i = init: steps
        for r = 1:4
            evecs = P{i} * vr{i};
            all_evecs{r}(:, i-init+1) = evecs(:, r);
        end
    end

    for r = 1:4
        rmax{r} = max(real(all_evecs{r}), [], 'all');
        rmin{r} = min(real(all_evecs{r}), [], 'all');
        imax{r} = max(imag(all_evecs{r}), [], 'all');
        imin{r} = min(imag(all_evecs{r}), [], 'all');
    end
    

    for i = init: steps
        for r = 1:4
            % view mode
            piece = reshape(all_evecs{r}(:, i-init+1), [nx, ny]);
    
            if tr
                piece = piece.';
            end
            
            subplot(4, 2, 2*r-1)
            pcolor(real(piece));
            clim([rmin{r}, rmax{r}]);
            shading interp
            colorbar()
            title('Real part')
            
            subplot(4, 2, 2*r)
            pcolor(imag(piece));
            clim([imin{r}, imax{r}]);
            shading interp
            colorbar()
            title('Imag part')
    
            sgtitle(['Mode # 1-4', ...
                ', with ', num2str(i), ' snapshots being used.'])
    
            pause(dt)
        end

        if gif_name
            frame = getframe(fig);  % 获取当前图像帧
            im = frame2im(frame); % 转换为图像
            [A, map] = rgb2ind(im, 256); % 转换为索引图像
            
            % 追加到 GIF
            if i == init
                imwrite(A, map, [gif_name, '.gif'], 'gif', 'LoopCount', inf, 'DelayTime', 0.5);
            else
                imwrite(A, map, [gif_name, '.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
            end
        end

    end

end






