function [W2prime,retargs]=sinkfastCost(l,C,restartU)
% SINKFAST Sinkhorn iteration on EXP(-C/L)  [ C=||X-Y||^2 ].
% [W,R]=SINKFAST(L,C) returns the Euclidean Wasserstein 
%   distance (squared) in W and further results in R.
%   Use fieldnames(R) for a list of available data.

% written by elmar.plischke@tu-clausthal.de
verbose=false; % true; % 
n=size(C,1);m=size(C,2); % assume uniform weights
if(l<0), l=-l*max(C(:)); end
% relative 
maxerr=1e-3;maxiter=1000;

%mx=mean(x);my=mean(y); % center 
%x=x-mx;y=y-my;
%Cneg=2*x*y'-sum(x.^2,2)-sum(y.^2,2)';

% entries in K will be wiped out by exp(- * /l)
lminr=max(min(C,[],2))/log(realmax); %row wipeout
lminc=max(min(C,[],1))/log(realmax); %column wipeout
rowcritical=lminr>lminc;
if(verbose)
 disp(['Critical lambdas ' num2str(lminr) ' (row) ' num2str(lminc) ' (col)']);
end
if isempty(l), l=max(lminr,lminc)*1.11; end

K=exp(-C/l);%clear('Cneg');
if(nargin<4) || isempty(restartU) || size(restartU,1)~=n
u=ones(n,1);
else
u=restartU;
end
err=inf;iter=0;

Ktu=K'*u;
while(iter==0 || isfinite(err) && err>maxerr && iter <maxiter)
 v=(1/m)./Ktu;
 if(~rowcritical), v(isinf(v))=1; end % avoid 0*inf = nan
 u=(1/n)./(K*v);
 if( rowcritical), u(isinf(u))=1; end % avoid 0*inf = nan
 Ktu=K'*u;
 err=sum(abs(v.*Ktu-1/m));
 iter=iter+1;
% if(isnan(err))
%  disp('time for breakpoint'); 
% end 
end
% sqrt <P,C>
% l*(mean(log(u))+mean(log(v)))
if(verbose)
 disp(sprintf('iteration %d, error %f\n',iter,err));
end

% potentials
f=l*log(u./(1/n));g=l*log(v./(1/m));
% u = exp(f/lambda).*alfa, v = exp(g/lambda).*beta compared to sinkhorn2
% cost (dual)
W22=(mean(f)+mean(g));

innerF=@(A,B)A(:)'*B(:); % Frobenius/trace inner product

P=v'.*K.*u; % = Gamma exp((f+g'-C)/l)
SinkhornDual=W22-l*sum(P(:)); 
% cost (primal) <P,C>
W2prime=innerF(P,C);
SinkhornPrimal=W2prime+l*innerF(P,log(P*(n*m)+1e-20)-1); % need dpi/dalfa x dbeta !

% put all intermediates and results in a struct
 retargs = struct('W2Dual',W22,'W2Primal',W2prime,'Kernel',K,...
           'Coupling',P,'ScalingU',u,'ScalingV',v,...
           'PotentialF',f,'PotentialG',g,...
           'NumIterations',iter','LastError',err,...
           'Lambda',l,'LambdasCritical',[lminr,lminc],...
           'SinkhornPrimal',SinkhornPrimal,'SinkhornDual',SinkhornDual);
end