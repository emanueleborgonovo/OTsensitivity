function [sigma]=wasserstein_bubble(X,Y,precision)
% Compute the Wasserstein distance
% Turn Puccetti 2017 into Bubblesort

% work with X' and Y' for 10% performace increase
if(nargin<3)
precision=1e-2;
end
%inner=@(x,y)x(:)'*y(:); % frobenius
n=size(X,1);
sigma=1:n;
%d=zeros(n,1);
X=X';
Ys=Y';
converged=false;
for m=n:-1:2
%while(~converged)
  Yo=Ys; 
  updated=false;
 for i=1:m-1
    xi=X(:,i);
    yi=Ys(:,i);
    for j=i+1:n
%      yj=Ys(:,j);
%        if ( xi-X(:,j))'*(yi-yj)>0
        if ( xi-X(:,j))'*(yi-Ys(:,j))>0
%            sigma([i j])=sigma([j i]); % direct swap
            si=sigma(i);sigma(i)=sigma(j);sigma(j)=si;
%            Ys(:,[i,j])=[yj,yi]; yi=yj;
            yj=yi;
            yi=Ys(:,j);
            %Ys(:,[i j])=Ys(:,[j i]);
            Ys(:,j)=yj;Ys(:,i)=yi;
            updated=true;
        end
    end
 end
% no change? optimum
 if ~updated
     disp('Optimum found'); 
     break;
	 %converged=true;
 else
% accuracy reached ?
%  for i=1:n
%    d(i)=X(i,:)*(Y(sigma(i),:)-Y(sigma_old(i),:))';
%  end
 d=sum(X.*(Ys-Yo),1);
 s=mean(abs(d))
  if(s<precision), disp('Precision reached'); %converged=true; 
   break;
  end
 end
 %plot([X(:,1),Y(sigma,1)]',[X(:,2),Y(sigma,2)]'); drawnow
 
 %%
end % while
%%
end
function wasserstein_swap_test
%%
n=1000;d=2;
X=rand(n,d);Y=rand(n,d);
tic
ii=wasserstein_swap(X,Y);
toc
subplot(2,1,1); plot([X(:,1),Y(ii,1)]',[X(:,2),Y(ii,2)]');
title(num2str(sum(mean(X.*Y(ii,:)))));
tic
jj=wasserstein_swap(-X,Y);
toc
subplot(2,1,2); plot([X(:,1),Y(jj,1)]',[X(:,2),Y(jj,2)]');
title(num2str(sum(mean(X.*Y(jj,:)))))
%%
end