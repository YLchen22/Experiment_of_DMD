function plot_snapshot(snapshot, nx, ny, tr)
%PLOT_SNAPSHOT 此处显示有关此函数的摘要
%   此处显示详细说明

    if nargin < 4
        tr = false;
    end

    piece = reshape(snapshot, [nx, ny]);

    if nx == 1
        piece = [piece; piece];
    elseif ny == 1
        piece = [piece, piece];
    end

    if tr
        piece = piece.';
    end

    figure()
    colormap hsv
    pcolor(real(piece));
    hold on
    shading interp
    colorbar()

end

