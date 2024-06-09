%% Incresing Sample Sizes 0.50 Corr
close all
%clear all
clc
tic
ns=[50 100 200 300 400 500 1000 10000 20000 50000 100000 200000];% 500000 700000];% 50000 100000 500000];% 1000 10000 20000 50000 100000 500000];% 1000 10000 20000 100000 500000];%
%m=[2 4 5 6 8 10 15 20 40 50];
m=[8 15 11 12 13 15 20 25 40 50 60 60];% 60 60];
%ns=ceil(logspace(1,4,15)); %(1 is the first power of 10; 4 the last power,
%15 intervals)
%
%model=@(x)x*w;
% Distribution
%ishigami
%ishigamifunction=model
dim=3;
p=sobolset(dim);
u=net(p,ns(end)+1);
xall=u*2*pi-pi;
yall=ishigamifunction(xall);
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
    W2(i,:)=W.W22/V;
   [Wass2,WassSep,SinkSep,Mean2Sep]=bwsi(xx,y,mm);
     Adv(i,:)=mean(Mean2Sep)/V;
     WB(i,:)=mean(WassSep)/V;
end
%%
toc

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
W2analyt=[0.423 0.423];
%W2analyt=[16 4 1];

figure
subplot(1,2,1)
xf=[1:length(ns)];
plot(xf,W2(:,1),'-o','Color','k','MarkerSize',10,'LineWidth',2)
hold on
plot(xf,W2(:,2),'-o','Color','r','MarkerSize',10,'LineWidth',2)
plot(xf,W2(:,3),'-o','Color','b','MarkerSize',10,'LineWidth',2)
hold on
plot(xf,WB(:,1),'-s','Color','k','MarkerSize',10,'LineWidth',2)
hold on
plot(xf,WB(:,2),'-s','Color','r','MarkerSize',10,'LineWidth',2)
hold on
plot(xf,WB(:,3),'-s','Color','b','MarkerSize',10,'LineWidth',2)
hold on
plot(xf,Adv(:,1),'-d','Color','k','MarkerSize',10,'LineWidth',2)
hold on
plot(xf,Adv(:,2),'-d','Color','r','MarkerSize',10,'LineWidth',2)
hold on
plot(xf,Adv(:,3),'-d','Color','b','MarkerSize',10,'LineWidth',2)


aa=max([max(W2(:,1)),max(W2(:,2)),max(W2(:,3))]);
axis([.5 length(ns)+0.5 0 max(aa)]);
set(gca,'xtick',[1:length(ns)])

set(gca,'xticklabel',num2str(ns'),'XTickLabelRotation',-45)
a=axis();
grid on
%plot([a(1) a(2)],W2analyt(1)*[1 1])
%plot([a(1) a(2)],W2analyt(2)*[1 1])
%plot([a(1) a(2)],W2analyt(3)*[1 1])
set(gca,'FontSize',16)
title('$\iota(Y,X_i)$, $\iota^{WB}(Y,X_i)$, $Adv(Y,X_i)$ for the Ishigami Function','FontSize',24,'interpreter','latex')
xlabel('Sample Size','FontSize',24,'interpreter','latex')
ylabel('$\iota(Y,X_i)$, $WB(Y,X_i)$, $Adv(Y,X_i)$','FontSize',24,'interpreter','latex')
legend({'$\iota(Y,X_1)$','$\iota(Y,X_2)$','$\iota(Y,X_3)$','$WB(Y,X_1)$','$WB(Y,X_2)$','$WB(Y,X_3)$','$Adv(Y,X_1)$','$Adv(Y,X_2)$','$Adv(Y,X_3)$'}, ...
       'FontSize',16,'interpreter','latex','NumColumns',3)
%%

%figure
subplot(1,2,2)
AA=[W2(end,:);WB(end,:);Adv(end,:)];
b = bar(AA); %barh(AA);
for i=1:3
xtips = b(i).XEndPoints;
ytips = b(i).YEndPoints;
%txtlabel = string(b(i).YData); % sprintf('%0.2g\n',..)
for j=1:3
text(xtips(j),ytips(j),sprintf('%0.2g',b(i).YData(j)),'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',16)
end
end
set(gca,'xticklabel',{'\iota(Y,X_i)', '\iota^{WB}(Y,X_i)', 'Adv(Y,X_i)'})
set(gca,'FontSize',24)
%% All in the same subplot
if(0)
figure
subplot(3,1,1)
plot(xf,W2(:,1),'-o','Color','k','MarkerSize',10,'LineWidth',2)
hold on
plot(xf,W2(:,2),'-o','Color','r','MarkerSize',10,'LineWidth',2)
plot(xf,W2(:,3),'-o','Color','b','MarkerSize',10,'LineWidth',2)
grid on
aa=max([max(W2(:,1)),max(W2(:,2)),max(W2(:,3))]);
axis([.5 length(ns)+0.5 0 max(aa)]);
set(gca,'xtick',[1:length(ns)])

set(gca,'xticklabel',num2str(ns'),'XTickLabelRotation',-45)
a=axis();
grid on

title('$\iota(Y,X_i)$, Wasserstein 2 squared','FontSize',24,'interpreter','latex')
xlabel('Sample Size','FontSize',24,'interpreter','latex')
ylabel('$\iota(Y,X_i)$,','FontSize',24,'interpreter','latex')
set(gca,'FontSize',16)
legend('$\iota(Y,X_1)$','$\iota(Y,X_2)$','$\iota(Y,X_3)$','FontSize',24,'interpreter','latex')

%
subplot(3,1,2)
plot(xf,WB(:,1),'-s','Color','k','MarkerSize',10,'LineWidth',2)
hold on
plot(xf,WB(:,2),'-s','Color','r','MarkerSize',10,'LineWidth',2)
hold on
plot(xf,WB(:,3),'-s','Color','b','MarkerSize',10,'LineWidth',2)
grid on

aa=max([max(W2(:,1)),max(W2(:,2)),max(W2(:,3))]);
axis([.5 length(ns)+0.5 0 max(aa)]);
set(gca,'xtick',[1:length(ns)])

set(gca,'xticklabel',num2str(ns'),'XTickLabelRotation',-45)
a=axis();
grid on

title('$\iota(Y,X_i)$, Wasserstein 2 squared','FontSize',24,'interpreter','latex')
xlabel('Sample Size','FontSize',24,'interpreter','latex')
ylabel('$WB(Y,X_i)$,','FontSize',24,'interpreter','latex')
set(gca,'FontSize',16)
legend('$WB(Y,X_1)$','$WB(Y,X_2)$','$WB(Y,X_3)$','$Adv(Y,X_1)$','$Adv(Y,X_2)$','$Adv(Y,X_3)$','FontSize',24,'interpreter','latex')

%
subplot(3,1,3)
plot(xf,Adv(:,1),'-o','Color','k','MarkerSize',10,'LineWidth',2)
hold on
plot(xf,Adv(:,2),'-o','Color','r','MarkerSize',10,'LineWidth',2)
plot(xf,Adv(:,3),'-o','Color','b','MarkerSize',10,'LineWidth',2)
grid on
title('$Adv(Y,X_i)=\xi^V{Y,X_i}$, ','FontSize',24,'interpreter','latex')
xlabel('Sample Size','FontSize',24,'interpreter','latex')
ylabel('$\iota(Y,X_i)$','FontSize',24,'interpreter','latex')
set(gca,'FontSize',16)
legend('$Adv(Y,X_1)$','$Adv(Y,X_2)$','$Adv(Y,X_3)$','FontSize',24,'interpreter','latex')


%plot([a(1) a(2)],W2analyt(1)*[1 1])
%plot([a(1) a(2)],W2analyt(2)*[1 1])
%plot([a(1) a(2)],W2analyt(3)*[1 1])

title('$\iota(Y,X_i)$, Wasserstein 2 squared','FontSize',24,'interpreter','latex')
xlabel('Sample Size','FontSize',24,'interpreter','latex')
ylabel('$\iota(Y,X_i)$, $WB(Y,X_i)$, $Adv(Y,X_i)$,','FontSize',24,'interpreter','latex')
set(gca,'FontSize',16)
legend('$\iota(Y,X_1)$','$\iota(Y,X_2)$','$\iota(Y,X_3)$','$WB(Y,X_1)$','$WB(Y,X_2)$','$WB(Y,X_3)$','$Adv(Y,X_1)$','$Adv(Y,X_2)$','$Adv(Y,X_3)$','FontSize',24,'interpreter','latex')

%%
end