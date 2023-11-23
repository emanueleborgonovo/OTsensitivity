%% Incresing Sample Sizes 0.50 Corr
close all
clear all
clc
tic
ns=[50 100 200 300 400 500 1000 10000 20000 50000 100000 500000];% 1000 10000 20000 50000 100000 500000];% 1000 10000 20000 100000 500000];%
%m=[2 4 5 6 8 10 15 20 40 50];
m=@(n) min(ceil(n^(2/(7+tanh((1500-n)/500)))),48);
%ns=ceil(logspace(1,4,15)); %(1 is the first power of 10; 4 the last power,
%15 intervals)
dint=333; % d=7; % d=1;
dim=3*dint;
w=kron([4,-2,1],ones(1,dint))';
%
%model=@(x)x*w;
% Distribution
rho=0.5;
Sig=rho*ones(dim,dim);
for i=1:dim
        Sig(i,i)=1;
end
mu=ones(1,dim);
%xall=mvnrnd(mu,Sig,ns(end));
%xall=mu+randn(ns(end),dim)*chol(Sig);
% takes ages, construct symmetric root by hand
% for this all-off diagonal constant 
v=ones(dim,1)/sqrt(dim);
Z=(v*v')*sqrt(1+(dim-1)*rho);
for i=1:dim-1
    v=[ones(i,1);
       -i;
       zeros(dim-i-1,1)];
    v=v/sqrt(i*(i+1));
    Z=Z+sqrt(1-rho)*(v*v');
end
% without the sqrt in line 29 and 35, Z=Sig; 
%xall=mu+randn(ns(end),dim)*Z;
%yall=xall*w;
I=[1 dint+1 2*dint+1];

xall=zeros(ns(end),length(I));
yall=zeros(ns(end),1);
if(~isempty(gcp('nocreate'))), parpool('threads'); end
parfor i=1:ns(end)
    %j=min(100,ns(end)-i+1);
    x0=mu+randn(1,dim)*Z;
    %xall(i+(1:j)-1,:)=x0(I);
    %yall(i+(1:j)-1)=x0*w;
    xall(i,:)=x0(I);
    yall(i)=x0*w;
end
V=var(yall);
%
tic
for r=1:20
%    WW=[];
%dd=[];
%SSi=[];
r
ii=randn(dim);
for i=1:length(ns)
    n=ns(i);
    mm=m(n);
    ii=randperm(ns(end),n-(n==ns(end))); % jackknife for last sample size
    y=yall(ii);
    xx=xall(ii,:); %,I);
    W=wassersi(xx,y,mm);
    W2(r,i,:)=W.W22/(2*V);    
end
end
toc
%%
%clf

%W1analyt=[433.144	430.024	431.583];
%W8analyt=[634.337	627.964	631.144];
%W16analyt=[773.582	765.077	769.319];
%danalyt=[0.306	0.303	0.305];
%Sanalyt=[0.501	0.495	0.498];

%figure
% Wass 1
%subplot(3,2,1)
%boxplot(W1(:,:,1),'plotstyle','compact')
%hold on
%boxplot(W1(:,:,2),'plotstyle','compact')
%hold on
%boxplot(W1(:,:,3),'plotstyle','compact')
%aa=max(W1(:,:,1));
%axis([.5 11.5 0 max(aa)]);
%set(gca,'xtick',[1:11])

%set(gca,'xticklabel',num2str(ns'),'XTickLabelRotation',-45)
%a=axis();
%plot([a(1) a(2)],W1analyt(1)*[1 1])
%plot([a(1) a(2)],W1analyt(2)*[1 1])
%plot([a(1) a(2)],W1analyt(3)*[1 1])

%title('\xi_i^{W1}, Wasserstein 1','FontSize',16)

%set(gca,'FontSize',16)
%legend('Analytical X1','Analytical X2','Analytical X3')
%%
%subplot(3,2,2)
W2analyt=[0.293 0.289 0.291];
%W2analyt=[16 4 1];
figure
%boxplot(W2(:,:,1),'plotstyle','compact')
%hold on
%boxplot(W2(:,:,2),'plotstyle','compact')
%hold on
%boxplot(W2(:,:,3),'plotstyle','compact')
boxplot(reshape(permute(W2,[1,3,2]),20,12*3),'PlotStyle','Compact','Color','bgr')
hold on
aa=max(W2(:,:,1));
axis([.5 3*12+.5 0 max(aa)]);
set(gca,'xtick',3*[1:length(ns)]-2.4); % -1 to center but because of space trick below

%set(gca,'xticklabel',num2str(ns'),'XTickLabelRotation',-45)
set(gca,'xticklabel',[repmat('     ',12,1),num2str(ns')],'XTickLabelRotation',-45)
a=axis();
plot([a(1) a(2)],W2analyt(1)*[1 1],'b')
plot([a(1) a(2)],W2analyt(2)*[1 1],'g')
plot([a(1) a(2)],W2analyt(3)*[1 1],'r')
for i=1:12
    plot((3*i+.5)*[1,1],[a(3) a(4)],':k');
end
title('$\iota(Y,X_i)$, Wasserstein 2 squared','FontSize',16,'interpreter','latex')
xlabel('Sample size','FontSize',24)
ylabel('$\iota(Y,X_i)$','FontSize',24,'interpreter','latex')
set(gca,'FontSize',16)
legend('i=1','i=2','i=3')
%legend('Analytical \iota(Y,X_i)$','Analytical X2','Analytical X3')

