function [sigma]=wasserstein_swap(X,Y)
% Compute the Wasserstein distance
% Puccetti 2017

precision=1e-4;
%inner=@(x,y)x(:)'*y(:); % frobenius
[n,d]=size(X);
sigma=1:n;
d=zeros(n,1);
Ys=Y;
while(1)
  Yo=Ys; 
  updated=false;
 for i=1:n-1
    xi=X(i,:);        % row vectors
    yi=Ys(i,:);
    for j=i+1:n
      %yj=Ys(j,:);
      %  if ( xi-X(j,:))*(yi-yj)'>0
      if ( xi-X(j,:))*(yi-Ys(j,:))'>0
            sigma([i j])=sigma([j i]); % direct swap
            yi=Ys(j,:);
            Ys([i,j],:)=Ys([j,i],:);
            updated=true;
        end
    end
 end
% no change? optimum
 if ~updated
     disp('Optimum found'); 
     break; 
 else
% accuracy reached ?
%  for i=1:n
%    d(i)=X(i,:)*(Y(sigma(i),:)-Y(sigma_old(i),:))';
%  end
 d=sum(X.*(Ys-Yo),2);
 s=mean(abs(d));
  if(s<precision), disp('Precision reached'); break; end
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