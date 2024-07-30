%一维信号DMD分解示意程序
clear
close all

x=0:0.01:5;m=length(x);%空间
t=0:0.05:5;N=length(t);%时间
Fs=1/(t(2)-t(1));%采样频率
[X,T]=meshgrid(x,t);

%创建信号,Utx的格式为时间*空间的矩阵
U_tx=1.2*exp(-0.5*T).*sin(2*pi*(X+2*T))+0.8*exp(0.3*T).*sin(2*pi*(3*X+4*T))+0.1*rand(N,m)+1.1;

%计算X和Y
U_xt=U_tx';
X=U_xt(:,1:end-1);
Y=U_xt(:,2:end);
%计算DMD
[Dd,b,Phi,Time_DMD,Energy]=DMD_CLASS(X,Y);

%之后全部是后处理
figure(1)
%图1，输出特征根分布
scatter(real(Dd(end:-1:1)),imag(Dd(end:-1:1)),30,-log(Energy(end:-1:1)),'filled')
axis equal
set(gcf,'position',[488   342   400   350])

figure(2)
%图2，绘制频率和衰减图
wa=log(Dd)*Fs;
scatter(real(wa(end:-1:1)),imag(wa(end:-1:1))/2/pi,30,-log(Energy(end:-1:1)),'filled')
xlabel('衰减率σ');ylabel('频率w')
ylim([-6,6]);xlim([-1,1])
hold on
plot([0,0],ylim,'b--')
plot(xlim,[0,0],'b--')
hold off
box on
set(gcf,'position',[488   342   400   350])

figure(3)
%图3，绘制频率-能量排序
Freq=imag(wa)/2/pi;
k=find(Freq>=0);
stem(Freq(k),log10(Energy(k)),'BaseValue',-6,'MarkerFaceColor','auto');
ylim([-6,6]);
xlim([-1,11])
set(gca,'YTickLabel',{'1e-6','1e-4','1e-2','0','1e2','1e4','1e6'},'YTick',[-6:2:6])
set(gca,'XTick',[0:2:10])
xlabel('频率hz');ylabel('能量')

figure(4)
%绘制模态
subplot(3,2,1)%总
U_x1=U_xt(:,1);
plot(x,U_x1);
xlim([0,5])
subplot(3,2,2)%一阶模态
Uxt_DMD_k=real(Phi(:,1) * Time_DMD(1,:));
plot(x,Uxt_DMD_k(:,1));
ylim([0,3]);xlim([0,5])
subplot(3,2,3)%二阶模态
Uxt_DMD_k=real(Phi(:,2) * Time_DMD(2,:));
plot(x,Uxt_DMD_k(:,1));
xlim([0,5])
subplot(3,2,4)%三阶模态
Uxt_DMD_k=real(Phi(:,3) * Time_DMD(3,:));
plot(x,Uxt_DMD_k(:,1));
xlim([0,5])
subplot(3,2,5)%四阶模态
Uxt_DMD_k=real(Phi(:,4) * Time_DMD(4,:));
plot(x,Uxt_DMD_k(:,1));
xlim([0,5])
subplot(3,2,6)%五阶模态
Uxt_DMD_k=real(Phi(:,5) * Time_DMD(5,:));
plot(x,Uxt_DMD_k(:,1));
xlim([0,5])

figure(5)
%绘制前10阶模态能量占比
Cumsum_Energy=cumsum(Energy);
subplot(2,1,1)
bar( 1:10 , Cumsum_Energy(1:10)/Cumsum_Energy(end) ,'BarWidth',1)
ylim([0,1]);
subplot(2,1,2)
plot( 1:10 , Energy(1:10)/Cumsum_Energy(end) ,'-o')
ylim([0,1]);

figure(6)
%还原信号与原信号对比(利用前5阶模态还原)
Uxt_DMD_k=real(Phi(:,1:5) * Time_DMD(1:5,:));
plot(x,U_xt(:,1),x,Uxt_DMD_k(:,1))

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
Dd=diag(D);%求特征值转为向量形式
%4 求DMD模态
Phi=Y*V/S*Om;
%5 求解初始状态b
b=Phi\X(:,1);
%6 对模态进行排序
Q=Dd.^(0:N-1);%计算出范德蒙矩阵
Time_DMD=b.*Q;%旧版本matlab可以用右边形式替换Time_DMD=(b*ones(1,N)).*Q;

% %对信号初始振幅排序，速度快，输出能量可以近似用振幅平方代替
% [b,Ib] = sort(abs(b),'descend');
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