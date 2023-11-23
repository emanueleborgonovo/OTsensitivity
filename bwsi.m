function [Wass2,WassSep,SinkSep,Mean2Sep]=bwsi(x,y,M,sigma,gfx)
% Bures-Wasserstein (squared) for Sensitivity

if(nargin<3) || isempty(M), M=8; end
if(nargin<4) || isempty(sigma), sigma=0; end
% unconditional
my=mean(y);Cy=cov(y);
%try
% Ry=sqrtm(Cy);
%catch
 Ry=mysqrtm(Cy);
%end
traceCy=trace(Cy); % precompute

[n,k]=size(x);
[nn,l]=size(y);

ms=round(linspace(0,n,M+1));
[~,ix]=sort(x);
W=zeros(M,k);
S=zeros(M,k);
Ms=zeros(M,k);
for m=1:M
 mc=ms(m+1)-ms(m);
 if(l>1)     
  for i=1:k
  % conditional output
   yc=y(ix(ms(m)+1:ms(m+1),i),:);

   mc=mean(yc);
   Cc=cov(yc);
     
   meanSqDists  = sum((my-mc).^2);
  % try
  %  traceTerms0 = traceCy+trace(Cc)-...
  %            2*trace(sqrtm(Ry*Cc*Ry));
  % catch
   traceTerms0 = traceCy+trace(Cc)-...
              2*tracesqrtm(Ry*Cc*Ry);
  % traceTerms1 = 2*(traceCy-tracesqrtm(Ry*Cc*Ry));
  % traceTerms0 = traceCy+trace(Cc)-...
  %             2*sum(svd(Ry*Rc));  % Rc not computed
  % end
   Ms(m,i)= meanSqDists;           
   W(m,i)= meanSqDists + traceTerms0;
   %S(m,i)= traceTerms1;
   if(sigma~=0)
    D = mysqrtm( 4*Ry*Cc*Ry + sigma^4*eye(size(Cy)));
    traceTerms = traceCy + trace(Cc) - trace(D);
    entropicTerms = sigma^2*(size(Cy,1)*(1-log(2*sigma^2))...
         ... was: log(det(D+sigma^2*eye(size(Cy)))))
         + trace(logm(D+sigma^2*eye(size(Cy)))));
    S(m,i)= meanSqDists + traceTerms + entropicTerms;
   end
  end
 else
  % scalar output
  for i=1:k
  % conditional output
   yc=y(ix(ms(m)+1:ms(m+1),i));

   mc=mean(yc);
   Cc=var(yc);
      
   Ms(m,i)= (my-mc).^2;           
   W(m,i)= (my-mc).^2+Cy+Cc-2*sqrt(Ry*Cc*Ry);
  end
 end
end
Wass2=mean(sqrt(W));
WassSep=W;
SinkSep=S;
Mean2Sep=Ms;
if(nargin>=5) && ~isempty(gfx)
 % figure
 plot(linspace(1/M,1,M),sqrt(W));
 title(gfx);ylabel('Bures Wasserstein Sensitivity')
 xlabel('Input quantile');
 legend(num2str((1:k)'))
end
end

function B=mysqrtm(A)
% MYSQRTM robust symmetric matrix square root
 [U,S,V]=svd(A,'econ');
 B=U*diag(sqrt(diag(S)))*V';
end

function b=tracesqrtm(A)
% TRACESQRTM sum of square roots of eigenvalues.
% D=max(eig(A),0); % remove spurious negative values
% b=sum(sqrt(D));
 b=sum(sqrt(svd(A)));
end