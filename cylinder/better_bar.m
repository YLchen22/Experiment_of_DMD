function better_bar(data)
% Input a vector containing only true or false, plot a bar figure to
% visualize it that True is better and False is worse at any index

% 将 true 和 false 转换为相应的高度值
heights = double(data) * 2 - 1;

% 分别找到 true 和 false 的索引
true_indices = find(data);
false_indices = find(~data);

% 创建一个新的图形窗口
figure()

% 绘制 true 的柱形图，设置颜色为绿色
hold on;
bar(true_indices, heights(true_indices), 'b');

% 绘制 false 的柱形图，设置颜色为红色
bar(false_indices, heights(false_indices), 'r');

% 设置 x 轴标签和 y 轴标签
xlabel('Index');
ylabel('Value');

% 设置 y 轴刻度范围（可选）
ylim([-1.5, 1.5]);

% 添加图例
legend('Better', 'worse', 'Location', 'best');

% 关闭 hold 状态
hold off;

end

