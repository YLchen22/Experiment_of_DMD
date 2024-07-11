%% POD_Cylinder Wake
clc
clear
close all
load CYLINDER_ALL.mat;
X=VORTALL';
Y=[X;X];
[U0x,An,phiU,Ds]=POD_SVD_M(Y);

%% 简单看一下瞬态计算结果，顺便保存个gif图片
pic_num = 1;
for i=1:100
    clf
    plotCylinder_m(reshape(VORTALL(:,i),nx,ny),nx,ny);
    %hold on
    pause(0.05)
    F=getframe(gcf);
    I=frame2im(F);
    [I,map]=rgb2ind(I,256);
    if pic_num == 1
        imwrite(I,map,'test4.gif','gif','Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(I,map,'test4.gif','gif','WriteMode','append','DelayTime',0.2);
    end
    pic_num = pic_num + 1;
end

%% 查看一下平均流场的信息，并保存图片
plotCylinder_m(reshape(U0x,nx,ny),nx,ny);
print(gcf, '-dpng', '-r600', './U0x.png');

%% 动态展示前十阶POD分解模态结果
for k=1:10
    clf
    plotCylinder_m(reshape(An(1,k).*phiU(:,k),nx,ny),nx,ny);
    drawnow
    pause(0.1)
end

%% 前二十模态能量大小(SVD具有自动排序功能，不需要再对特征值大小进行排序)
figure(1)
plot(1:20,Ds(1:20)/sum(Ds),'o--');
ylabel('Energy order','FontName','Arial', 'FontSize',14);
xlabel('k','FontName','Arial', 'FontSize',14);
grid on
name1=['Energy order.png'];
print(gcf, '-dpng', '-r600',name1);
figure(2)
semilogy(1:20,cumsum(Ds(1:20))/sum(Ds),'o--');
ylabel('Total Energy','FontName','Arial', 'FontSize',14);
xlabel('k','FontName','Arial', 'FontSize',14);
grid on
name2=['Energy cum.png'];
print(gcf, '-dpng', '-r600',name2);

%% 前6阶模态叠加还原流场
Sigma=zeros(size(phiU));
for i=1:6
    V{i}=An(:,i).*phiU(:,i)';
    Sigma=Sigma+V{i}';
end
for j=1:100
    clf
    plotCylinder_m(reshape(Sigma(:,j)'+U0x,nx,ny),nx,ny);
    pause(0.05)
end
