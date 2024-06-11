close all
clear all
clc
tic
r12=0.5;
r13=0.5;
r23=0.5;
muX=[1 1 1]
SigX=[1 r12 r13; r12 1 r23; r13 r23 1];
b=[0 0]';
A=[4 -2 1; 2 5 -1];
%n=8192;
%X=mvnrnd(muX,SigX,n);
ns=[50 100 250 500 1000];% 1000 10000]% 10000]; %100000];
mm=[2  5   7   8   10];
filename='ResultsGaussNormalPaper1SimplexSmallGood.xlsx';
k=3;
%%
lambda=1;
tic
for q=1:length(ns)
    n=ns(q)
    M=mm(q)
u=sobolpoints(n,k);
X=norminv(u)*chol(SigX)+muX;
Y=(A*X'+b)';
%M=mm;%[Wass2,WassSep,SinkSep,Mean2Sep]=bwsi(X,Y,M,lambda);
     [result]=bigwassersteintest(X,Y,M,2^0+2^1+2^2+2^6);
     results{q}=result;
%placing results in vectors
VY2=2*sum(var(Y));
W2analyt=[0.492 0.507 0.117];

   Wsimplexresn(q,:)=mean(results{q}.Wsimplex)/VY2;
    Wsinkhornn(q,:)=mean(results{q}.Wsinkhorn)/VY2;
    Wswapn(q,:)=mean(results{q}.Wswap)/VY2;
    %Wswapupscalen(q,:)=mean(results{q}.Wswapupscale)/VY2;
    %Whungariann(q,:)=mean(results{q}.Whungarian)/VY2;
    %Whungarianupscalen(q,:)=mean(results{q}.Whungarianupscale)/VY2;
    Wburesn(q,:)=mean(results{q}.Wbures)/VY2;
    %Wgradientn(q,:)=mean(results{q}.Wgradient)/VY2;
    %Wdivergencen(q,:)=mean(results{q}.Wdivergence)/VY2;
    %Wsimplexcomplementn(q,:)=mean(results{q}.Wsimplexcomplement)/VY2;
    W2analytvn(q,:)=W2analyt;
%Results reporting at each cycle to save them in case of memory failure
    
xlswrite(filename,Wsimplexresn,'Wsimplexresn')
xlswrite(filename,Wsinkhornn,'Wsinkhornn')
xlswrite(filename,Wswapn,'Wswapn')
%xlswrite(filename,Wswapupscalen,'Wswapupscalen')
%xlswrite(filename,Whungariann,'Whungariann')
%xlswrite(filename,Whungarianupscalen,'Whungarianupscalen')
xlswrite(filename, Wburesn,' Wburesn')
%xlswrite(filename,Wgradientn,'Wgradientn')
%xlswrite(filename,Wdivergencen,'Wdivergencen')
%xlswrite(filename,Wsimplexcomplementn,'Wsimplexcomplementn')
xlswrite(filename,W2analytvn,'W2analytic')

%Times    
    Tsimplexresn(q,:)=results{q}.Tsimplex;
   Tsinkhornn(q,:)=results{q}.Tsinkhorn;
   Tswapn(q,:)=results{q}.Tswap;
  %  Tswapupscalen(q,:)=mean(sqrt(results{q}.Tswapupscale));
   % Thungariann(q,:)=mean(sqrt(results{q}.Thungarian));
    %Thungarianupscalen(q,:)=mean(sqrt(results{q}.Thungarianupscale));
    Tburesn(q,:)=mean(results{q}.Tbures);
    %Tgradientn(q,:)=mean(sqrt(results{q}.Tgradient));
    %Tdivergencen(q,:)=mean(sqrt(results{q}.Tdivergence));
    %Tsimplexcomplementn(q,:)=mean(sqrt(results{q}.Tsimplexcomplement));
    %Results reporting at each cycle to save them in case of memory failure
    T=[Tsimplexresn';
        Tsinkhornn';Tswapn'%;Tswapupscalen';Thungariann';Thungarianupscalen'
        ;Tburesn'%;Tgradientn';Tdivergencen';Tsimplexcomplementn'
        ];
xlswrite(filename,T,'T')

end
toc
%% Subset of plots
figure

subplot(4,1,1)
plot(ns,Wsimplexresn,'--d', 'LineWidth',1.5)
hold on
set(gca,'ColorOrderIndex',1)
plot(ns,W2analytvn, 'LineWidth',1.5)
title('Simplex','FontSize',14)
xlabel('N','FontSize',14)
ylabel('\iota(Y,X_i)','FontSize',14)


subplot(4,1,2)
plot(ns,Wsinkhornn,'--d', 'LineWidth',1.5)
hold on
set(gca,'ColorOrderIndex',1)
plot(ns,W2analytvn, 'LineWidth',1.5)
title('Sinkhorn','FontSize',14)
xlabel('N','FontSize',14)
ylabel('\iota(Y,X_i)','FontSize',14)

%legend('\xi_1 Estimate','\xi_2 Estimate','\xi_3 Estimate','\xi_1 Analytical','\xi_2 Analytical','\xi_3 Analytical')
subplot(4,1,3)
plot(ns,Wswapn,'--d', 'LineWidth',1.5)
hold on
set(gca,'ColorOrderIndex',1)
plot(ns,W2analytvn, 'LineWidth',1.5)
title('Swap','FontSize',14)
xlabel('N','FontSize',14)
ylabel('\iota(Y,X_i)','FontSize',14)

subplot(4,1,4)
plot(ns,Wburesn,'--d', 'LineWidth',1.5)
hold on
set(gca,'ColorOrderIndex',1)
plot(ns,W2analytvn, 'LineWidth',1.5)
title('Bures','FontSize',14)
xlabel('N','FontSize',14)
ylabel('\iota(Y,X_i)','FontSize',14)

legend('\iota(Y,X_1) Estimate','\iota(Y,X_2) Estimate','\iota(Y,X_3) Estimate','\iota(Y,X_1) Analytical','\iota(Y,X_2) Analytical','\iota(Y,X_3) Analytical','Fontsize',18)


