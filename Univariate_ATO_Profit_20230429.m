%% Incresing Sample Sizes 0.50 Corr
close all
%clear all
clc
load('ATOdataInvyProfitTime2.mat')
tic
ns=length(y);% 500000 700000];% 50000 100000 500000];% 1000 10000 20000 50000 100000 500000];% 1000 10000 20000 100000 500000];%
%m=[2 4 5 6 8 10 15 20 40 50];
m=[20];% 60 60];
%ns=ceil(logspace(1,4,15)); %(1 is the first power of 10; 4 the last power,
%15 intervals)
%
%model=@(x)x*w;
% Distribution
%ishigami
%ishigamifunction=model
xall=x;
yall=y;
V=var(yall);
%    WW=[];
%dd=[];
%SSi=[];
for i=1:length(ns)
    i;
    n=ns(i);
    mm=m(i);
    y=yall(1:n);
    xx=xall(1:n,:);
    W=wassersi(xx,y,mm);
    W2(i,:)=W.W22/(2*V);
   [Wass2,WassSep,SinkSep,Mean2Sep]=bwsi(x,y,mm)
     Adv(i,:)=mean(Mean2Sep)/(2*var(y));
     WB(i,:)=mean(WassSep)/(2*var(y));
end
%%
toc
%%

figure1 = figure('WindowState','maximized');

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

AA=[W2' WB' Adv'];
bar(AA)
xticks=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22];
set(axes1,'FontSize',14,'XTickLabel',...
    {'X_2','X_4','X_6','X_{8}','X_{10}','X_{12}','X_{14}','X_{16}','X_{18}','X_{20}','X_{22}'});

title('Global Sensitivity Measures for the ATO Univariate Output','FontSize',30,'interpreter','latex')
legend('$\iota(Y,X_i)$','$\iota^{WB}(Y,X_i)$','$Adv(Y,X_i)$','FontSize',30,'interpreter','latex')
xlabel('Input','FontSize',30,'interpreter','latex')
ylabel('$\widehat{\iota}(Y,X_i)$,  $\widehat{\iota^{WB}}(Y,X_i)$,  $\widehat{Adv}(Y,X_i)$','FontSize',30,'interpreter','latex')
