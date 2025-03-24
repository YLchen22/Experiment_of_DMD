clc; clear; close all;

% 读取已有 GIF 文件
old_filename = 'modes_of_cylinder.gif';  % 你的 GIF 文件
new_filename = 'updated_animation.gif';  % 新的 GIF 文件

% 设定新的帧间隔时间（单位：秒）
new_delay_time = 0.5;

% 获取 GIF 信息
gif_info = imfinfo(old_filename);
num_frames = numel(gif_info);  % 获取帧数

% 重新保存 GIF，修改帧间隔
for k = 1:num_frames
    % 逐帧读取 GIF
    [im, map] = imread(old_filename, 'gif', 'Frames', k);
    
    % 写入新的 GIF 文件并修改帧间隔
    if k == 1
        imwrite(im, map, new_filename, 'gif', 'LoopCount', inf, 'DelayTime', new_delay_time);
    else
        imwrite(im, map, new_filename, 'gif', 'WriteMode', 'append', 'DelayTime', new_delay_time);
    end
end

disp('新的 GIF 生成完成，帧间隔已调整。');
