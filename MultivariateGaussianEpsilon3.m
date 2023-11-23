close all
clear all
clc

r12=0.5;
r13=0.5;
r23=0.5;
muX=[1 1 1];
SigX=[1 r12 r13; r12 1 r23; r13 r23 1];
b=[0 0]';
A=[4 -2 1; 2 5 -1];
%n=8192;
%ns=[50 100 200 300 400 500 1000 10000 20000 30000];% 500000];% 1000 10000 20000 100000 500000];%
%m=[2 4 5 6 8 10 15 20 40 50];

mm=@(n) min(ceil(n^(2/(7+tanh((1500-n)/500)))),48);


%% Sample generation
N=100000;%[50 100 200 300 400 500 1000 10000 20000 50000 100000];
XN=mvnrnd(muX,SigX,N);
save XN.txt -ascii XN;


%%

W2analyt=[6.47 6.52 2.86];
aux=0;

lambda=sqrt(0.1);
for q=1:length(ns)
    n=ns(q)
    M=mm(n)
%M=mm;

Y=(A*X'+b)';
[n,k]=size(X);
[nn,l]=size(Y);
ms=round(linspace(0,n,M+1));

if(n~=nn), error('size mismatch'); end
[~,ix]=sort(X);

W=zeros(M,k);
for m=1:M
    mc=ms(m+1)-ms(m);     
    for r=1:k
        aux=aux+1;
        tic
    % conditional output
     yc=Y(ix(ms(m)+1:ms(m+1),r),:);
     % compute Sinkhorn (p=2)  
     [W2,retargs]=sinkfast(lambda,Y,yc);
     WW2(m,r)=sqrt(max(0,W2));
     SkD(m,r)=sqrt(max(0,retargs.SinkhornDual));
     SkP(m,r)=sqrt(max(0,retargs.SinkhornPrimal));
     numiterations(q,m,r)=retargs.NumIterations;
     tocc=toc;
     telapsedfast(aux,:)=[tocc,M];
     samplesizes(aux)=n;
    end
end
    
WW2m(q,:)=mean(WW2);
SkDm(q,:)=mean(SkD);
SkPm(q,:)=mean(SkP);
end
Estimates=WW2m

tGaussSkFast=telapsedfast;

 %%
%xlswrite('MultivGaussResultsAll30000_t.xlsx',WW2m4,'Sk4Wass');
%xlswrite('MultivGaussResultsAll30000_t.xlsx',SkPm4,'Sk4Primal');
%xlswrite('MultivGaussResultsAll30000_t.xlsx',SkDm4,'Sk4Dual');
%xlswrite('MultivGaussResultsAll30000_t.xlsx',tMultivGaussSk4,'tMultivGaussSk4');
%xlswrite('MultivGaussResultsAll30000_t.xlsx',samplesizes4,'samplesizes4');

%% Sinkhorn4
ns=[5000];
aux=0;
lambda=1;
for q=1:length(ns)
    n=ns(q)
    M=mm(n)

%M=mm;
X=mvnrnd(muX,SigX,n);
Y=(A*X'+b)';
[n,k]=size(X);
[nn,l]=size(Y);
ms=round(linspace(0,n,M+1));

if(n~=nn), error('size mismatch'); end
[~,ix]=sort(X);

W=zeros(M,k);
for m=1:M
    mc=ms(m+1)-ms(m);     
    for r=1:k
        r
        tic
        aux=aux+1;
    % conditional output
     yc=Y(ix(ms(m)+1:ms(m+1),r),:);
     % compute Sinkhorn (p=2)  
     [W2,retargs]=sinkhorn4([],lambda,Y,yc);
     WW2(m,r)=sqrt(max(0,W2));
     SkD(m,r)=sqrt(max(0,retargs.SinkhornDual));
     SkP(m,r)=sqrt(max(0,retargs.SinkhornPrimal));
     numiterations(q,m,r)=retargs.NumIterations;
     tocc=toc;
     telapsed4(aux,:)=[tocc,m];
     samplesizes4(aux)=n;
    end
end 
WW2m4(q,:)=mean(WW2);
SkDm4(q,:)=mean(SkD);
SkPm4(q,:)=mean(SkP);
end
 tMultivGaussSk4=telapsed4;

Estimates=WW2m4

fhh=figure
plot(ns,WW2m4,'-o',ns,sqrt(SkPm4),'-.d','LineWidth',2)
legend('\xi^{W2}_1','\xi^{W2}_2','\xi^{W2}_3','\xi^{Sk}_1','\xi^{Sk}_2','\xi^{Sk}_3')
xlabel('Sample size, N')
ylabel('OT-based Importance Measures')
title('Sinkhorn with numerical stabilization')
set(gca,'FontSize',18)
  
%saveas(fhh,'MultivGaussSk4increasing30000','eps')
%saveas(fhh,'MultivGaussSk4increasing30000','fig')
%saveas(fhh,'MultivGaussSk4increasing30000','jpg')

