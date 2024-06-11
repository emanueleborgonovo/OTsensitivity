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
ns=[50000 50000 50000 50000 50000 50000 50000 50000 50000 50000 50000 50000];% 1000 10000]% 10000]; %100000];
mm=[5 10 20 30 40 50 60 70 80 100 150 200];
filename='ResultsGaussNormalPaperpartitionsEntropic.xlsx';
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
     [dummy,WassSep,SinkSep,Mean2Sep]=bwsi(X,Y,M,lambda);
     
%placing results in vectors
VY2=2*sum(var(Y));
W2analyt=[0.492 0.507 0.117];
W2analytentr=[0.554 0.575 0.199];

   % Wsimplexresn(q,:)=mean(results{q}.Wsimplex)/VY2;
%    Wsinkhornn(q,:)=mean(results{q}.Wsinkhorn)/VY2;
    %Wswapn(q,:)=mean(results{q}.Wswap)/VY2;
    %Wswapupscalen(q,:)=mean(results{q}.Wswapupscale)/VY2;
    %Whungariann(q,:)=mean(results{q}.Whungarian)/VY2;
    %Whungarianupscalen(q,:)=mean(results{q}.Whungarianupscale)/VY2;
    Wburesn(q,:)=mean(WassSep)/VY2;
    WCuturin(q,:)=mean(SinkSep)/VY2;
    %Wgradientn(q,:)=mean(results{q}.Wgradient)/VY2;
    %Wdivergencen(q,:)=mean(results{q}.Wdivergence)/VY2;
    %Wsimplexcomplementn(q,:)=mean(results{q}.Wsimplexcomplement)/VY2;
    W2analytvn(q,:)=W2analyt;
    W2analytentrn(q,:)=W2analytentr;
%Results reporting at each cycle to save them in case of memory failure
    
%xlswrite(filename,Wsimplexresn,'Wsimplexresn')
%xlswrite(filename,Wsinkhornn,'Wsinkhornn')
%xlswrite(filename,Wswapn,'Wswapn')
%xlswrite(filename,Wswapupscalen,'Wswapupscalen')
%xlswrite(filename,Whungariann,'Whungariann')
%xlswrite(filename,Whungarianupscalen,'Whungarianupscalen')
xlswrite(filename, Wburesn,' Wburesn')
xlswrite(filename,WCuturin,'WCuturin')
%xlswrite(filename,Wdivergencen,'Wdivergencen')
%xlswrite(filename,Wsimplexcomplementn,'Wsimplexcomplementn')
xlswrite(filename,W2analytvn,'W2analytic')

%Times    
 %   Tsimplexresn(q,:)=mean(sqrt(results{q}.Tsimplex));
 %  Tsinkhornn(q,:)=mean(sqrt(results{q}.Tsinkhorn));
  % Tswapn(q,:)=mean(sqrt(results{q}.Tswap));
  %  Tswapupscalen(q,:)=mean(sqrt(results{q}.Tswapupscale));
   % Thungariann(q,:)=mean(sqrt(results{q}.Thungarian));
    %Thungarianupscalen(q,:)=mean(sqrt(results{q}.Thungarianupscale));
%    Tburesn(q,:)=mean(sqrt(results{q}.Tbures));
    %Tgradientn(q,:)=mean(sqrt(results{q}.Tgradient));
    %Tdivergencen(q,:)=mean(sqrt(results{q}.Tdivergence));
    %Tsimplexcomplementn(q,:)=mean(sqrt(results{q}.Tsimplexcomplement));
    %Results reporting at each cycle to save them in case of memory failure
 %   T=[%Tsimplexresn';
        %Tsinkhornn';Tswapn'%;Tswapupscalen';Thungariann';Thungarianupscalen';
  %      Tburesn'%;Tgradientn';Tdivergencen';Tsimplexcomplementn'
       % ];
%xlswrite(filename,T,'T')

end
toc
%% Subset of plots
figure
%subplot(1,2,1)
plot(mm,Wburesn,'--d', 'LineWidth',1.5)
hold on
set(gca,'ColorOrderIndex',1)
plot(mm,W2analytvn, 'LineWidth',1.5)
col=get(gca,'ColorOrderIndex');
plot(mm,WCuturin,'--o', 'LineWidth',1.5)
set(gca,'ColorOrderIndex',col)
plot(mm,W2analytentrn, 'LineWidth',1.5)
hold off
title('Bures','FontSize',14)
xlabel('M','FontSize',14)
ylabel('Estimates of \iota(Y,X_i)','FontSize',14)



