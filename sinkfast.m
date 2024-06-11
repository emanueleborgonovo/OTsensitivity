function [W2prime,retargs]=sinkfast(l,x,y,restartU)
% SINKFAST Sinkhorn iteration on EXP(-C/L) with C=||X-Y||^2.
% [W,R]=SINKFAST(L,X,Y) returns the Euclidean Wasserstein 
%   distance (squared) in W and further results in R.
%   Use fieldnames(R) for a list of available data.

% written by elmar.plischke@tu-clausthal.de
verbose=  false; % true; % 
n=size(x,1); % assume uniform weights
m=size(y,1);

maxerr=1e-3;maxiter=1000;
mx=mean(x);my=mean(y); % center 
x=x-mx;y=y-my;
Cneg=2*x*y'-sum(x.^2,2)-sum(y.^2,2)';
if(l<0), l=l*min(Cneg(:)); end % relative 
% entries in K will be wiped out by exp(- * /l)
lminr=-min(max(Cneg,[],2))/log(realmax); %row wipeout
lminc=-min(max(Cneg,[],1))/log(realmax); %column wipeout
rowcritical=lminr>lminc;
if(verbose)
 disp(['Critical lambdas ' num2str(lminr) ' (row) ' num2str(lminc) ' (col)']);
end
if isempty(l), l=max(lminr,lminc)*1.11; end

K=exp(Cneg/l);%clear('Cneg');
if(nargin<4) || isempty(restartU) || size(restartU,1)~=n
u=ones(n,1);
else
u=restartU;
end
err=inf;iter=0;

Ktu=K'*u;
if(~rowcritical)
 while(iter==0 || isfinite(err) && err>maxerr && iter <maxiter)
  v=(1/m)./Ktu;
  v(isinf(v))=1; % avoid 0*inf = nan
  u=(1/n)./(K*v);
  Ktu=K'*u;
  err=sum(abs(v.*Ktu-1/m));
  iter=iter+1;
 end
else
 while(iter==0 || isfinite(err) && err>maxerr && iter <maxiter)
  v=(1/m)./Ktu;
  u=(1/n)./(K*v);
  u(isinf(u))=1; % avoid 0*inf = nan
  Ktu=K'*u;
  err=sum(abs(v.*Ktu-1/m));
  iter=iter+1; 
 end
end
% sqrt <P,C>
% l*(mean(log(u))+mean(log(v)))
if(verbose)
 disp(sprintf('iteration %d, error %f\n',iter,err));
end

% potentials
f=l*log(u./(1/n));g=l*log(v./(1/m));
% u = exp(f/lambda).*alfa, v = exp(g/lambda).*beta compared to sinkhorn2
meanSqDiff=sum((mx-my).^2);
% cost (dual)
W22=meanSqDiff+(mean(f)+mean(g));

innerF=@(A,B)A(:)'*B(:); % Frobenius/trace inner product

P=v'.*K.*u; % = Gamma exp((f+g'-C)/l)
SinkhornDual=W22-l*sum(P(:)); 
% cost (primal) <P,C>
W2prime=meanSqDiff-innerF(P,Cneg);
SinkhornPrimal=W2prime+l*innerF(P,log(P*(n*m)+1e-20)-1); % need dpi/dalfa x dbeta !

% put all intermediates and results in a struct
 retargs = struct('W2Dual',W22,'W2Primal',W2prime,'Kernel',K,...
           'Cost',-Cneg,'Coupling',P,'ScalingU',u,'ScalingV',v,...
           'PotentialF',f,'PotentialG',g,'MeanSqShift',meanSqDiff,...
           'NumIterations',iter','LastError',err,...
           'Lambda',l,'LambdasCritical',[lminr,lminc],...
           'SinkhornPrimal',SinkhornPrimal,'SinkhornDual',SinkhornDual);
end