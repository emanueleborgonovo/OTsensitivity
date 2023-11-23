close all
clear all
%% Increasing Sample Sizes 0.50 Corr
clc

ns=[50 100 200 300 400 500 1000 10000 20000 50000 100000 500000];% 1000 10000 20000 50000 100000 500000];% 1000 10000 20000 100000 500000];%
%m=[2 4 5 6 8 10 15 20 40 50];
m=@(n) min(ceil(n^(2/(7+tanh((1500-n)/500)))),48);
%ns=ceil(logspace(1,4,15)); %(1 is the first power of 10; 4 the last power,
%15 intervals)+
%aux=0;
figure
%for dint=[1 3 33 333 3333];
%    aux=aux+1;
%    aux
dint=3333;
    % d=7; % d=1;
dim=3*dint;
w=kron([4,-2,1],ones(1,dint))';
%
%model=@(x)x*w;
% Distribution
rho=0.5;
Sig=rho*ones(dim)+(1-rho)*diag(ones(dim,1));

mu=ones(1,dim);
disp('sample creation')
tic
xall=mvnrnd(mu,Sig,ns(end));
toc
yall=xall*w;
V=var(yall);

I=[1 dint+1 2*dint+1];
xall=xall(:,I); % throw away all unused input dims

%%
R=20;
disp('replicates')
W2=zeros(R,length(ns),length(I));
for r=1:R
    r
%    WW=[];
% dd=[];
% SSi=[];
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
%%


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
%W2analyt=[0.571833333	0.008366667	0.143666667
%0.357194444	0.089034722	0.192048611
%0.296740995	0.257194587	0.276443682
%0.293233831	0.289054726	0.291144279
%0.293	0.293	0.293];
%%
W2analyt=[.293 .293 .293 ]; % GUESSING
boxplot(reshape(permute(W2,[1,3,2]),20,length(ns)*3),'PlotStyle','Compact','Color','bgr')
hold on
aa=max(W2(:,:,1));
axis([.5 3*length(ns)+.5 0 max(aa)]);
set(gca,'xtick',3*[1:length(ns)]-2.4); % -1 to center but because of space trick below

%set(gca,'xticklabel',num2str(ns'),'XTickLabelRotation',-45)
set(gca,'xticklabel',[repmat('     ',length(ns),1),num2str(ns')],'XTickLabelRotation',-45)
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