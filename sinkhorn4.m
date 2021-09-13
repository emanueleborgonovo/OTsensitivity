function [OT,retargs]=sinkhorn4(C,lambda,x,y,alfa,beta)
% SINKHORN4 Optimal Transport via symmetric Sinkhorn algorithm.
%
% [W,R]=SINKFAST(C,L,X,Y,P,Q) returns the Euclidean Wasserstein 
%   distance (squared) in W and further results in R.
%   Use fieldnames(R) for a list of available data
%   Inputs: C   cost function c(x,y) (default squared Euclidean) yielding N x M matrix
%	          L   regularization
%           X,Y     sample locations (sizes N and M)
%           P,Q     weight vectors (sizes N and M) (uniform defaults)
% partly taken from Jean Feydy 2020

 maxiter=1000;
 miniter=0;
 maxerr=1e-5;
 verbose=false;%true; %false; % 
 
% assume: 
% sum(p)==sum(q)==1
% size(c)==[n,n], all(all(c>=0)), size(p)==[n,1],size(q)==[n,1]
 n=size(x,1);
 m=size(y,1);
 if(nargin<5)
  alfa=ones(n,1)/n; % uniform masses
  beta=ones(m,1)/m;
 end
 % Euclidean norm squared
 if isempty(C)
   C=@(x,y)sum(x.^2,2)+sum(y.^2,2)'-2*x*y';
 end
 cxy=C(x,y);nc=norm(cxy);
 cxx=C(x,x);
 cyy=C(y,y);
 
 Delta=max(sqrt(cxy(:))); % diameter
 
 if(isempty(lambda))
  alfai=min([alfa;beta]); 
  if(alfai>0)
   lambda=-1/10/log(alfai);
  else
   lambda=max(cxy(:));
  end
  disp(['SINKHORN: lambda set to ' num2str(lambda) ]);
 end

 f = cxy*beta;
 g = (alfa'*cxy)';
 ff = cxx*alfa;
 gg = cyy*beta;
 
% v0=zeros(n,1);
% u0=zeros(n,1);

 converged=false;
 Tlold=0;
 Tl=inf();
 iter=0;
 while ~converged
    f_ =.5*(f-lambda*logsumexp((g'-cxy)/lambda,2,beta));
    g  =.5*(g-lambda*logsumexp((f-cxy)/lambda,1,alfa)');
    f  = f_; % fake parallel execution
    ff =.5*(ff-lambda*logsumexp((ff'-cxx)/lambda,2,alfa));
    gg =.5*(gg-lambda*logsumexp((gg-cyy)/lambda,1,beta)');
  iter=iter+1;
  Tlold=Tl;
  Tl=((f-ff)'*alfa+(g-gg)'*beta);
  err=abs(Tl-Tlold);
  if verbose   
   KK=exp((-cxy+f+g')/lambda);
   err2=mean(alfa'*KK-1);
   % err3=mean(KK*beta-1); % same as err2 in symmetric approach
  disp(sprintf('errors %g %g, iter %d',err,err2, iter)); 
  end
  % relative error dependent on cost range
  converged=(iter>miniter && err<maxerr*Delta) || isnan(Tl) || iter==maxiter;
 end
 if(isnan(Tl))
     OT=Tlold;
     T2=0;
 else
     OT=((f-ff)'*alfa+(g-gg)'*beta);
     T2=(f'*alfa+g'*beta); % without divergence correction
 end
 if(nargout>1)
  K=exp(-cxy/lambda);
  P=exp((f+g'-cxy)/lambda); % double-stochastic
  Q=beta'.*P.*alfa;         % prescibed marginals
  W2prime=(Q(:)'*cxy(:));
  SinkhornPrimal=W2prime+lambda*alfa'*(P.*(log(P+1e-20)-1))*beta; 
  SinkhornDual=T2-lambda*alfa'*P*beta; % sum(P(:)); 
   retargs = struct('SinkhornDivergence',OT,'W2Dual',T2,...
      'W2Primal',W2prime,'Kernel',K,'Lambda',lambda,...
      'PotentialF',f,'PotentialG',g,'PotentialFF',ff,'PotentialGG',gg,...
      'NumIterations',iter','LastError',err,...
      'Cost',cxy,'Coupling',Q,...
      'SinkhornPrimal',SinkhornPrimal,'SinkhornDual',SinkhornDual);
 end
end
function y=expz(x)
% EXPZ exp with expz(-infty)=0
% B. Schmitzer
 ii=isinf(x) & x<0;
 y=exp(x);
 y(ii)=0;
end


function testsink
%%
n=100;
a=randn(n,1);
b=randn(n,1)*2+1;
W22=mean((sort(a)-sort(b)).^2)

W22s=sinkhorn3([],[],a,b,ones(n,1)/n,ones(n,1)/n)

%%
end

function testsink2
%%
ishigami
n=1024;
x=trafo(sobolpoints(n,k));
y=model(x);
D=(y-y').^2/2;
p=ones(n,1)/n;
M=10;
xt=trafo(zeros(1,k));
for m=1:M
  m
  xu=trafo(ones(1,k)*m/M);
  % select between xt and xu
  ii=x<=xu & x>xt;
  q=ii./sum(ii);
  for i=1:k
    W(i,m)=sinkhorn2(D,.1,p,q(:,i)); %no light speed here
  end
  %
 % tic;
 % W(:,m)=sinkhorn(D*2,4,p,q);toc
  xt=xu;
end
plot(1:M,W)
mean(sqrt(W),2)
%%
end