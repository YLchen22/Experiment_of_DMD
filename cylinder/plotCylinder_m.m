function f1=plotCylinder_m(VORT)
%2022.2.19
%流场后处理
%输入变量VORT，为空间变量转化为列向量的流场快照矩阵
%输入变量ny，nx，为空间尺度大小
pcolor(VORT);
shading interp;
load CCcool.mat
colormap(CC);
caxis([-min(abs(min(min(VORT))),max(max(VORT))),min(abs(min(min(VORT))),max(max(VORT)))]);

%重置坐标轴
set(gca,'XTick',[1 50 100 150 200 250 300 350 400 449],'XTickLabel',{'-1','0','1','2','3','4','5','6','7','8'})
set(gca,'YTick',[1 50 100 150 199],'YTickLabel',{'2','1','0','-1','-2'});
set(gca,'Layer','top');
ylabel('y','FontName','Arial', 'FontSize',14);
xlabel('x','FontName','Arial', 'FontSize',14);
axis equal
hold on

% 增加云图控制线
contour(VORT,[linspace(-max(max(VORT)),-max(max(VORT))/35,6)],'-k','LineWidth',1)
contour(VORT,[linspace(max(max(VORT))/35,max(max(VORT)),6)],'--k','LineWidth',1)
theta = (1:100)/100'*2*pi;
x = 49+25*sin(theta);
y = 99+25*cos(theta);
fill(x,y,[.3 .3 .3])  % place cylinder
plot(x,y,'k','LineWidth',1.2) % cylinder boundary
set(gcf,'PaperPositionMode','auto') %

colorbar()