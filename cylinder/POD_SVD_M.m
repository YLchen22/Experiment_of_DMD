function [U0x,An,phiU,Ds]=POD_SVD_M(Utx)
%2022.2.18
%SVD-POD算法
%输入Utx，其中时间离散长度N=size(Utx,1)，空间离散长度m=size(Utx,2)
%输出U0x，0阶模态，可以看做定常平均值
%输出An：时间变量，对应模态的幅值随时间的变化，可以用来做时间序列分析
%输出phiU：POD的模态
%输出Ds：特征值Ds反映了每一个模态对应的能量，可以用来排序

m=size(Utx,2);%空间长度
N=size(Utx,1);%时间长度

%除去稳态值
U0x=mean(Utx,1);
Utx=Utx-U0x.*ones(N,m);
[U,S,phiU]=svd(Utx,'econ');
An=U*S;
Ds=diag(S).^2/N;
end
