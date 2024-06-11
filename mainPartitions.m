%% Main tutORial
close all
%clear all
clc
%% FInite Change Indices (Code for Figures 3, 4, 5, 6)
x0=[2 1 3];
xplus=[3 2 4];
xminus=[1 0.5 2];
%% Uncertainty Quantification (Figure 7)
k=3;
N=2^18;
p=haltonset(k);
u=net(p,N);
x=u.*(xplus-xminus)+xminus;
yEOQ=x(:,1).*x(:,2).*x(:,3);
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
% Conditional Output Distribution
x1=2*ones(N,1);
xx=x;
xx(:,1)=x1;
yEOQ1=xx(:,1).*xx(:,2).*xx(:,3);
Ey1=mean(yEOQ1)
Vy1=var(yEOQ1)
[f1,xi1] = ksdensity(yEOQ1);
text(1.2,0.75,'E[EOQ|X_1=2]=0.78','FontSize',24,'Color','r')
text(1.2,0.65,'V[EOQ|X_1=2]=0.08','FontSize',24,'Color','r')
xline(prctile(yEOQ1,5),'r')
xline(prctile(yEOQ1,95),'r')
hold on
plot(xi1,f1,'-.r','LineWidth',2)
%% Global Sensitivity Measures (Figure 8)



%% Pearson's given data intuition and conditiona densities in Figure 9 
figure
subplot(2,2,1)
scatter(x(:,1),yEOQ)
for i=1:10
    xline(prctile(x(:,1),i*10),'LineWidth',2)
end
xlabel('$X_i$','FontSize',24,'interpreter','latex')
ylabel('$Y$','FontSize',24,'interpreter','latex')
title('Graph (a): H=10','FontSize',24,'interpreter','latex')

subplot(2,2,3)
for i=1:10
hold on
aa=find(and(x(:,1)>prctile(x(:,1),(i-1)*10),x(:,1)<prctile(x(:,1),(i)*10)));
[f,xi] = ksdensity(yEOQ(aa));
plot(xi,f,'b','LineWidth',2)
hold on 
[f,xi] = ksdensity(yEOQ);
plot(xi,f,'k','LineWidth',3)
end
title('Conditional Densities for Graph (a)','FontSize',24,'interpreter','latex')
xlabel('$Y$','FontSize',24,'interpreter','latex')
ylabel('$f_{Y|X_1\in \mathcal{X}_i^m}, m=1,2,...,10$','FontSize',24,'interpreter','latex')

subplot(2,2,2)
scatter(x(:,1),yEOQ)
for i=1:100
    xline(prctile(x(:,1),i),'LineWidth',1)
end
xlabel('$X_i$','FontSize',24,'interpreter','latex')
ylabel('$Y$','FontSize',24,'interpreter','latex')
title('Graph (b): H=100','FontSize',24,'interpreter','latex')

subplot(2,2,4)
for i=1:100
hold on
aa=find(and(x(:,1)>prctile(x(:,1),(i-1)),x(:,1)<prctile(x(:,1),(i))));
[f,xi] = ksdensity(yEOQ(aa));
plot(xi,f,'b','LineWidth',1)
hold on 
[f,xi] = ksdensity(yEOQ);
plot(xi,f,'k','LineWidth',3)
end
title('Conditional Densities for Graph (b)','FontSize',24,'interpreter','latex')
xlabel('$Y$','FontSize',24,'interpreter','latex')
ylabel('$f_{Y|X_i\in \mathcal{X}_i^m}, m=1,2,...,100$','FontSize',24,'interpreter','latex')

%% Plot for estimator convergence (Figure 10)
% Varying both partition cardinality M and sample size N

m=[5 10 20 30 50 100];
N=[200 1000 10000 20000 50000 200000];
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






