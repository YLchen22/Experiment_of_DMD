%自己试一试DMD
clear
close all

%定义流场时间和空间信息
x=-8:0.2:8;
y=-5:0.2:5;
t=0:0.05:6;
Fs=1/(t(2)-t(1));
[X,Y,T]=meshgrid(x,y,t);
[Ny,Nx,Nt]=size(X);

%自定义流场,U是沿x方向的速度分量，V是y方向的
U0=-1*Y.^2+5;
V0=0*Y;%UV0不随时间变化，设为定常流场
U1=-5*sin(Y).*cos(2*pi*0.3*Y).*(exp(T/5));
V1=5*sin(X-1*pi*T).*cos(2*pi*0.3*Y).*(exp(T/5));
U3=0.01*rand(Ny,Nx,Nt);
V3=0.01*rand(Ny,Nx,Nt);

U_Sum=U0+U1+U3;
V_Sum=V0+V1+V3;

%1计算UV向量
U_xt=Uxyt_to_Uxt(U_Sum);%把2维问题转化为1维问题
V_xt=Uxyt_to_Uxt(V_Sum);%把2维问题转化为1维问题
%然后把UV向量合并
UV_xt=[U_xt;V_xt];
%之前的这些都属于数据准备和整理部分


%计算X和Y
UX=UV_xt(:,1:end-1);
UY=UV_xt(:,2:end);
%计算DMD
[Dd,b,Phi,Time_DMD,Energy]=DMD_CLASS(UX,UY);

%后处理
figure(1)
%图1，绘制频率和衰减图
wa=log(Dd)*Fs;
scatter(real(wa(end:-1:1)),imag(wa(end:-1:1))/2/pi,30,-log(Energy(end:-1:1)),'filled')
xlabel('衰减率');ylabel('频率')
ylim([-6,6]);xlim([-1,1])
hold on
plot([0,0],ylim,'b--')
plot(xlim,[0,0],'b--')
hold off
box on
set(gcf,'position',[488   342   400   350])

figure(2)
%图3，绘制频率-能量排序
Freq=imag(wa)/2/pi;
k=find(Freq>=0);
stem(Freq(k),log10(Energy(k)),'BaseValue',-4,'MarkerFaceColor','auto');
ylim([-4,8]);
xlim([-1,11])
set(gca,'YTickLabel',{'1e-4','1e-2','0','1e2','1e4','1e6','1e8'},'YTick',[-4:2:8])
set(gca,'XTick',[0:2:10])
xlabel('频率hz');ylabel('能量')

figure(3)
%绘制前10阶模态能量占比
Cumsum_Energy=cumsum(Energy);
subplot(2,1,1)
bar( 1:10 , Cumsum_Energy(1:10)/Cumsum_Energy(end) ,'BarWidth',1)
ylim([0,1]);
subplot(2,1,2)
plot( 1:10 , Energy(1:10)/Cumsum_Energy(end) ,'-o')
ylim([0,1]);



figure(4)
X=X(:,:,1);
Y=Y(:,:,1);

k=1;
%绘制模态
subplot(2,2,1)%总
[Uxy0,Vxy0]=UV2UxyVxy(UV_xt(:,k),Ny,Nx);
hold on
pcolor(X,Y,curl(X,Y,Uxy0,Vxy0));shading interp
quiver(X(1:5:end,1:5:end),Y(1:5:end,1:5:end),Uxy0(1:5:end,1:5:end),Vxy0(1:5:end,1:5:end),'color','k')
hold off
axis equal off
title('总')

subplot(2,2,2)%1
UVxt_DMD_k=real(Phi(:,1) * Time_DMD(1,:));
[Uxyk,Vxyk]=UV2UxyVxy(UVxt_DMD_k(:,k),Ny,Nx);
hold on
pcolor(X,Y,curl(X,Y,Uxyk,Vxyk));shading interp
quiver(X(1:5:end,1:5:end),Y(1:5:end,1:5:end),Uxyk(1:5:end,1:5:end),Vxyk(1:5:end,1:5:end),'color','k')
hold off
axis equal off
title('1阶')

subplot(2,2,3)%2
UVxt_DMD_k=real(Phi(:,2) * Time_DMD(2,:));
[Uxyk,Vxyk]=UV2UxyVxy(UVxt_DMD_k(:,k),Ny,Nx);
hold on
pcolor(X,Y,curl(X,Y,Uxyk,Vxyk));shading interp
quiver(X(1:5:end,1:5:end),Y(1:5:end,1:5:end),Uxyk(1:5:end,1:5:end),Vxyk(1:5:end,1:5:end),'color','k')
hold off
axis equal off
title('2阶')

subplot(2,2,4)%3&4
UVxt_DMD_k=real(Phi(:,3) * Time_DMD(3,:));
[Uxy3,Vxy3]=UV2UxyVxy(UVxt_DMD_k(:,k),Ny,Nx);
UVxt_DMD_k=real(Phi(:,4) * Time_DMD(4,:));
[Uxy4,Vxy4]=UV2UxyVxy(UVxt_DMD_k(:,k),Ny,Nx);
Uxyk=Uxy3+Uxy4;
Vxyk=Vxy3+Vxy4;
hold on
pcolor(X,Y,curl(X,Y,Uxyk,Vxyk));shading interp
quiver(X(1:5:end,1:5:end),Y(1:5:end,1:5:end),Uxyk(1:5:end,1:5:end),Vxyk(1:5:end,1:5:end),'color','k')
hold off
axis equal off
title('3+4阶')




function [Dd,b,Phi,Time_DMD,Energy]=DMD_CLASS(X,Y)
%DMD分解函数
%输入：
%X，Y，DMD分解的数据矩阵
%
%输出：
%Dd，特征根
%b，初始状态
%Phi，DMD分解的模态
%Time_DMD，DMD还原信号所用到的时间项
%Energy，每个模态对应的能量，从大到小排序

N=size(X,2);
%1 SVD分解
[U,S,V] = svd(X,'econ');
%删除奇异值约等于0的模态，防止计算发散
Sd=diag(S);
r=sum(Sd>1e-6);
U=U(:,1:r);
S=S(1:r,1:r);
V=V(:,1:r);
%2 求解矩阵A
A=U'*(Y*V/S);
%3 求矩阵A的特征值和特征向量
[Om,D]=eig(A);
%4 求DMD模态
Phi=Y*V/S*Om;

Dd=diag(D);%求特征值转为向量形式

%5 求解初始状态b
b=Phi\X(:,1);

Q=Dd.^(0:N-1);%计算出范德蒙矩阵
Time_DMD=b.*Q;%旧版本matlab可以用右边形式替换Time_DMD=(b*ones(1,N)).*Q;


% %对信号初始振幅排序，速度快
% [~,Ib] = sort(abs(b),'descend');
%求出所有模态的信号，并计算能量
Energy=zeros(size(Phi,2),1);
for k=1:size(Phi,2)
Uxt_DMD_k=real(Phi(:,k) * Time_DMD(k,:));
E_k=sum(sum(Uxt_DMD_k.^2));
Energy(k)=E_k;
end
[Energy,Ie] = sort(Energy,'descend');%对每个模态的能量进行排序
%按顺序输出特征值、初始状态b、模态Phi
Dd=Dd(Ie);
b=b(Ie);
Phi=Phi(:,Ie);
Time_DMD=Time_DMD(Ie,:);
end

function Uxt=Uxyt_to_Uxt(Uxyt)
%把3维矩阵的xyt压缩为xt
[Ny,Nx,Nt]=size(Uxyt);%这里是matlab的meshgrid定义导致的，x和y是反着的，注意。
Nxy=Ny*Nx;
Uxt=reshape(Uxyt,Nxy,Nt);
end

function [Uxy,Vxy]=UV2UxyVxy(UVx,Ny,Nx)
Ux=UVx(1:Ny*Nx);
Vx=UVx(Ny*Nx+1:end);
Uxy=reshape(Ux,Ny,Nx);
Vxy=reshape(Vx,Ny,Nx);
end
