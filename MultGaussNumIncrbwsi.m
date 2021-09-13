close all
clear all
clc

r12=0.5;
r13=0.5;
r23=0.5;
muX=[1 1 1]
SigX=[1 r12 r13; r12 1 r23; r13 r23 1];
b=[0 0]';
A=[4 -2 1; 2 5 -1];
%n=8192;
%X=mvnrnd(muX,SigX,n);
ns=[100 400 500 1000 10000 20000 100000];
mm=@(n) min(ceil(n^(2/(7+tanh((1500-n)/500)))),48);

k=3;
%%
lambda=1;
tic
for q=1:length(ns)
    n=ns(q)
    M=mm(n)
u=sobolpoints(n,k);
X=norminv(u)*chol(SigX)+muX;
Y=(A*X'+b)';
%M=mm;
     [Wass2,WassSep,SinkSep,Mean2Sep]=bwsi(X,Y,M,lambda);
     WW(q,:)=mean(sqrt(WassSep));
     Sk(q,:)=mean(sqrt(SinkSep));
end
tbures=toc
%xlswrite('W2Wasseropt20000.xlsx',S2)
WW
Sk
%%
figure
W2analyt=[6.47 6.52 2.86];
hold on
% Wass 1
aa=max(WW(:,:));
axis([.5 size(WW,1)+1 0 max(aa)+0.6]);
plot([1:size(WW,1)],WW,'d','LineWidth',2)
set(gca,'xtick',[1:size(WW,1)])

set(gca,'xticklabel',num2str(ns'),'XTickLabelRotation',-45)
a=axis();
hold on
plot([a(1) a(2)],W2analyt(1)*[1 1])
plot([a(1) a(2)],W2analyt(2)*[1 1])
plot([a(1) a(2)],W2analyt(3)*[1 1])
legend({'Numerical X_1','Numerical X_2','Numerical X_3','Analytical X_1','Analytical X_2','Analytical X_3'},'NumColumns',2)
xlabel('Sample Size')
ylabel('\xi^{W2} using the Wasserstein-Bures estimation')
set(gca,'FontSize',24)

% analytic:  6.47 6.52 2.86
