%% Main ATO unconditional distribution
close all
clear all
clc
%% FInite Change Indices (Code for Figures 3, 4, 5, 6)
x0=[2 1 3];
xplus=[3 2 4];
xminus=[1 0.5 2];
%% Uncertainty Quantification (Figure 7)
load('ATOdataInvyProfitTime2.mat')
N=size(y,1);
yEOQ=y;
Ey=mean(yEOQ)
Vy=var(yEOQ)
% Unconditional Output Distribution
figure
[f,xi] = ksdensity(yEOQ);
plot(xi,f,'b','LineWidth',2)
xlabel('$y$','Interpreter','latex','FontSize',28)
ylabel('$\hat{f}_Y(y),\hat{f}_{Y|X_1}(y)$','FontSize',28,'Interpreter','latex')
title('Empirical Density of the EOQ','FontSize',20)
text(0.70,0.75,'E[EOQ]=0.77','FontSize',24,'Color','b')
text(0.70,0.65,'V[EOQ]=0.10','FontSize',24,'Color','b')
xline(prctile(yEOQ,5),'b')
xline(prctile(yEOQ,95),'b')




%% Pearson's given data intuition and conditiona densities in Figure 9 
figure
subplot(2,4,1)
scatter(x(:,20),yEOQ)
%for i=1:10
%    xline(prctile(x(:,20),i*10),'LineWidth',2)
%end
xlabel('$X_{20}$','FontSize',24,'interpreter','latex')
ylabel('$Y$','FontSize',24,'interpreter','latex')
title('Graph (a): $X_{20}$','FontSize',24,'interpreter','latex')

subplot(2,4,5)
for i=1:10
hold on
aa=find(x(:,20)==i);
[f,xi] = ksdensity(yEOQ(aa));
plot(xi,f,'b','LineWidth',2)
hold on 
[f,xi] = ksdensity(yEOQ);
plot(xi,f,'k','LineWidth',3)
end
title('Conditional Densities','FontSize',24,'interpreter','latex')
xlabel('$Y$','FontSize',24,'interpreter','latex')
ylabel('$f_{Y|X_{20}\in \mathcal{X}_{20}^m}$','FontSize',24,'interpreter','latex')

subplot(2,4,2)
scatter(x(:,18),yEOQ)
%for i=1:10
%    xline(prctile(x(:,18),i*10),'LineWidth',2)
%end
xlabel('$X_{18}$','FontSize',24,'interpreter','latex')
ylabel('$Y$','FontSize',24,'interpreter','latex')
title('Graph (b): $X_{18}$','FontSize',24,'interpreter','latex')

subplot(2,4,6)
for i=1:10
hold on
aa=find(x(:,18)==i);
[f,xi] = ksdensity(yEOQ(aa));
plot(xi,f,'b','LineWidth',1)
hold on 
[f,xi] = ksdensity(yEOQ);
plot(xi,f,'k','LineWidth',3)
axis([-200 300 0 0.014])
end
title('Conditional Densities','FontSize',24,'interpreter','latex')
xlabel('$Y$','FontSize',24,'interpreter','latex')
ylabel('$f_{Y|X_{18}\in \mathcal{X}_{18}^m}$','FontSize',24,'interpreter','latex')

subplot(2,4,3)
scatter(x(:,1),yEOQ)
for i=1:10
    xline(prctile(x(:,1),i*10),'LineWidth',1)
end
xlabel('$X_1$','FontSize',24,'interpreter','latex')
ylabel('$Y$','FontSize',24,'interpreter','latex')
title('Graph (c): $X_1$','FontSize',24,'interpreter','latex')

subplot(2,4,7)
for i=1:10
hold on
aa=find(and(x(:,1)>prctile(x(:,1),(i-1)),x(:,1)<prctile(x(:,1),(i))));
[f,xi] = ksdensity(yEOQ(aa));
axis([-200 300 0 0.014])
plot(xi,f,'b','LineWidth',1)
hold on 
[f,xi] = ksdensity(yEOQ);
plot(xi,f,'k','LineWidth',3)
end
title('Conditional Densities','FontSize',24,'interpreter','latex')
xlabel('$Y$','FontSize',24,'interpreter','latex')
ylabel('$f_{Y|X_1\in \mathcal{X}_1^m}$','FontSize',24,'interpreter','latex')

subplot(2,4,4)
scatter(x(:,3),yEOQ)
for i=1:10
    xline(prctile(x(:,3),i*10),'LineWidth',1)
end
xlabel('$X_{3}$','FontSize',24,'interpreter','latex')
ylabel('$Y$','FontSize',24,'interpreter','latex')
title('Graph (d): $X_3$','FontSize',24,'interpreter','latex')

subplot(2,4,8)
for i=1:10
hold on
aa=find(and(x(:,3)>prctile(x(:,3),(i-1)),x(:,3)<prctile(x(:,3),(i))));
[f,xi] = ksdensity(yEOQ(aa));
axis([-200 300 0 0.014])
plot(xi,f,'b','LineWidth',1)
hold on 
[f,xi] = ksdensity(yEOQ);
plot(xi,f,'k','LineWidth',3)
end
title('Conditional Densities','FontSize',24,'interpreter','latex')
xlabel('$Y$','FontSize',24,'interpreter','latex')
ylabel('$f_{Y|X_{3}\in \mathcal{X}_3^m}$','FontSize',24,'interpreter','latex')



%% Plot for estimator convergence (Figure 10)
% Varying both partition cardinality M and sample size N

m=[5 10 15 25];
N=[200 1000 10000 16000];
for i=1:length(m)
  [diaux,etaiaux]=deltamim(x(1:N(i),:),yEOQ(1:N(i)),m(i));
dm(i,:)=diaux;
etam(i,:)=etaiaux;
end
figure
subplot(2,1,1)
title('Asymptotic Convergence','FontSize',18)
plot(N,etam,'-d','LineWidth',2)
xlabel('Sample Size N','FontSize',24)
ylabel('$\hat{\eta}_i,i=1,2,3$','FontSize',24,'Interpreter','latex')
subplot(2,1,2)
plot(m,dm,'-d','LineWidth',2)
xlabel('Sample Size N','FontSize',24)
ylabel('$\hat{\delta}_i,i=1,2,3$','FontSize',24,'Interpreter','latex')
legend('$X_1$','$X_2$','$X_3$','FontSize',24,'interpreter','latex')


% Varying the partition cardinality M at fixed N
m=[5 10 20 30 50 100];
for i=1:length(m)
  [diaux,etaiaux]=deltamim(x,yEOQ,m(i));
dm(i,:)=diaux;
etam(i,:)=etaiaux;
end
figure
subplot(2,1,1)
title('Asymptotic Convergence','FontSize',18)
plot(m,etam,'-d','LineWidth',2)
xlabel('Partition Size','FontSize',24)
ylabel('$\hat{\eta}_i,i=1,2,3$','FontSize',24,'Interpreter','latex')
subplot(2,1,2)
plot(m,dm,'-d','LineWidth',2)
xlabel('Partition Size','FontSize',24)
ylabel('$\hat{\delta}_i,i=1,2,3$','FontSize',24,'Interpreter','latex')
legend('$X_1$','$X_2$','$X_3$','FontSize',24,'interpreter','latex')



%% Alternative Figure
figure

subplot(1,4,1)
for i=1:10
hold on
aa=find(x(:,20)==i);
[f,xi] = ksdensity(yEOQ(aa));
plot(xi,f,'b','LineWidth',2)
hold on 
[f,xi] = ksdensity(yEOQ);
plot(xi,f,'k','LineWidth',3)
end
title('Graph (a): $X_{20}$','FontSize',24,'interpreter','latex')
xlabel('$Y$','FontSize',24,'interpreter','latex')
ylabel('$f_{Y|X_{20}\in \mathcal{X}_{20}^m}$','FontSize',24,'interpreter','latex')

subplot(2,4,2)
scatter(x(:,18),yEOQ)
%for i=1:10
%    xline(prctile(x(:,18),i*10),'LineWidth',2)
%end
xlabel('$X_{18}$','FontSize',24,'interpreter','latex')
ylabel('$Y$','FontSize',24,'interpreter','latex')
title('Graph (b): $X_{18}$','FontSize',24,'interpreter','latex')

subplot(1,4,2)
for i=1:10
hold on
aa=find(x(:,18)==i);
[f,xi] = ksdensity(yEOQ(aa));
plot(xi,f,'b','LineWidth',1)
hold on 
[f,xi] = ksdensity(yEOQ);
plot(xi,f,'k','LineWidth',3)
axis([-200 300 0 0.014])
end
title('Graph (b): $X_{18}$','FontSize',24,'interpreter','latex')
xlabel('$Y$','FontSize',24,'interpreter','latex')
ylabel('$f_{Y|X_{18}\in \mathcal{X}_{18}^m}$','FontSize',24,'interpreter','latex')


subplot(1,4,3)
for i=1:10
hold on
aa=find(and(x(:,1)>prctile(x(:,1),(i-1)),x(:,1)<prctile(x(:,1),(i))));
[f,xi] = ksdensity(yEOQ(aa));
axis([-200 300 0 0.014])
plot(xi,f,'b','LineWidth',1)
hold on 
[f,xi] = ksdensity(yEOQ);
plot(xi,f,'k','LineWidth',3)
end
title('Graph (c): $X_1$','FontSize',24,'interpreter','latex')
xlabel('$Y$','FontSize',24,'interpreter','latex')
ylabel('$f_{Y|X_1\in \mathcal{X}_1^m}$','FontSize',24,'interpreter','latex')


subplot(1,4,4)
for i=1:10
hold on
aa=find(and(x(:,3)>prctile(x(:,3),(i-1)),x(:,3)<prctile(x(:,3),(i))));
[f,xi] = ksdensity(yEOQ(aa));
axis([-200 300 0 0.014])
plot(xi,f,'b','LineWidth',1)
hold on 
[f,xi] = ksdensity(yEOQ);
plot(xi,f,'k','LineWidth',3)
end
title('Graph (d): $X_3$','FontSize',24,'interpreter','latex')
xlabel('$Y$','FontSize',24,'interpreter','latex')
ylabel('$f_{Y|X_{3}\in \mathcal{X}_3^m}$','FontSize',24,'interpreter','latex')


