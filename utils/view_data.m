function view_data(data, nx, ny, tr, gif_name, dt)
%VIEW_DATA : visualize snapshots data as matrix on a given field
%   data: snapshots data matrix
%   nx, ny: x and y field
%   tr: set as TRUE if data needs transpose
%   gif_name: If you want to save a gif, input the wanted gif name
%   dt: time to pause between two snapshots

    if nargin < 4
        tr = false;
    end

    if nargin < 5
        gif_name = False;
    end

    if nargin < 6
        dt = 0.01;
    end

    dmin = min(data, [], "all");
    dmax = max(data, [], "all");

    fig = figure();
    % use a better color setting
    colormap hsv
    % load CCcool.mat
    % colormap(CC);
    % load batlow.mat
    % colormap(batlow)

    % view snapshots
    for i = 1: size(data, 2)
        piece = reshape(data(:, i), [nx, ny]);

        if nx == 1
            piece = [piece; piece];
        elseif ny == 1
            piece = [piece, piece];
        end

        if tr
            piece = piece.';
        end

        pcolor(real(piece));
        clim([dmin, dmax]);
        hold on
        shading interp
        colorbar()

        if (nx ~= 1) && (ny ~= 1)
            if dmin * dmax > 0
                contour(piece, [linspace(dmin, dmax, 10)],'--k','LineWidth',1)
            else
                contour(piece, [linspace(dmin, dmin/20, 5)],'--k','LineWidth',1)
                contour(piece, [linspace(dmax/20, dmax, 5)],'--k','LineWidth',1)
            end
        end
        
        title(['Step = ', num2str(i)])
        pause(dt)

        if gif_name
            frame = getframe(fig);  % 获取当前图像帧
            im = frame2im(frame); % 转换为图像
            [A, map] = rgb2ind(im, 256); % 转换为索引图像
            
            % 追加到 GIF
            if i == 1
                imwrite(A, map, [gif_name, '.gif'], 'gif', 'LoopCount', inf, 'DelayTime', 0.1);
            else
                imwrite(A, map, [gif_name, '.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
            end
        end
    end

end






