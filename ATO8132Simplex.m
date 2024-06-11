close all
clear all
clc


%n=8192;
%X=mvnrnd(muX,SigX,n);
%%
load('ATO_OT_8192.mat')
filename='ATO8132.xlsx';
%%
[result]=bigwassersteintest(x(1:2000,:),Inventory(1:2000,:),7,2^1+2^2+2^6+2^11);
%placing results in vectors
%%
VY2=sum(var(Inventory))*2;
    Wsimplexresn=mean(result.Wrevsimplex)/VY2;
    Wsinkhornn=mean(result.Wsinkhorn)/VY2;
    Wswapn=mean(result.Wswap)/VY2;
    %Wswapupscalen(q,:)=mean(results{q}.Wswapupscale)/VY2;
    %Whungariann(q,:)=mean(results{q}.Whungarian)/VY2;
    %Whungarianupscalen(q,:)=mean(results{q}.Whungarianupscale)/VY2;
    Wburesn=mean(result.Wbures)/VY2;
    %Wgradientn(q,:)=mean(results{q}.Wgradient)/VY2;
    %Wdivergencen(q,:)=mean(results{q}.Wdivergence)/VY2;
    %Wsimplexcomplementn(q,:)=mean(results{q}.Wsimplexcomplement)/VY2;
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

%Times    
    Tsimplexresn=mean(result.Trevsimplex);
   Tsinkhornn=mean(result.Tsinkhorn);
   Tswapn=mean(result.Tswap);
  %  Tswapupscalen(q,:)=mean(sqrt(results{q}.Tswapupscale));
   % Thungariann(q,:)=mean(sqrt(results{q}.Thungarian));
    %Thungarianupscalen(q,:)=mean(sqrt(results{q}.Thungarianupscale));
    Tburesn=mean(result.Tbures);
    %Tgradientn(q,:)=mean(sqrt(results{q}.Tgradient));
    %Tdivergencen(q,:)=mean(sqrt(results{q}.Tdivergence));
    %Tsimplexcomplementn(q,:)=mean(sqrt(results{q}.Tsimplexcomplement));
    %Results reporting at each cycle to save them in case of memory failure
    T=[Tsimplexresn';
        Tsinkhornn';Tswapn';%Tswapupscalen';Thungariann';Thungarianupscalen'
        Tburesn'%;Tgradientn';Tdivergencen';Tsimplexcomplementn'
        ];
xlswrite(filename,T,'T')

%% Subset of plots
figure
A=[Wsimplexresn;Wsinkhornn;Wswapn;Wburesn];
bar(A)
%subplot(4,1,1)
%plot(ns,Wsimplexresn,'--d', 'LineWidth',1.5)
%hold on
%set(gca,'ColorOrderIndex',1)
%plot(ns,W2analytvn, 'LineWidth',1.5)
%title('Simplex','FontSize',14)
%xlabel('N','FontSize',14)
%ylabel('\iota(Y,X_i)','FontSize',14)

%legend('\iota(Y,X_1) Estimate','\iota(Y,X_2) Estimate','\iota(Y,X_3) Estimate','\iota(Y,X_1) Analytical','\iota(Y,X_2) Analytical','\iota(Y,X_3) Analytical','Fontsize',18)

