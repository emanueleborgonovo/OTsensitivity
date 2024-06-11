function [W,D,E]=wassgrad2(x,C,M,l_or_opts)
% WASSGRAD2 Entropic Wasserstein Sensitivity with a gradient method.
% [W,D,E]=WASSGRAD(X,C,M,L)
% INPUTS
% X input matrix (N x K), 
% C cost matrix (N x N) or output vector (N x I, I~=N)
% M partition size
% L scaling factor for entropy, negative for relative to max(cost)
% OUTPUTS
% W primal separation matrix (M x K)

 [n,k]=size(x);
 if(size(C,2)~=n)
 % assume vector/multivariate output and no cost matrix
  C=-2*(C*C')+(sum(C.^2,2)+sum(C.^2,2)'); % Euclidean squared
 end
 if nargin<3 || isempty(M)
  M=ceil(n^(1/3));
 end
% if(nargin<4 || isempty(l))
%     l=max(C(:))*.002; %0.05 % *.0015 exploits floating point range.
% end
% 
 
 %% process options
% defaults
opts=struct('StepSize',10,...
	        'MaxIter',1000,...
            'MaxError',1e-4,...
            'Lambda',-.002,...
            'Verbose',false);

if(nargin>=4) && ~isempty(l_or_opts)
    if isstruct(l_or_opts)
        members=fieldnames(opts);
        for i=1:length(members)
            o=members{i};
            if isfield(l_or_opts,o), opts.(o)=l_or_opts.(o);end
        end
    else
        opts.Lambda=l_or_opts;
    end
end
l=opts.Lambda;
if(l<0) % interpret negative as relative 
    if(opts.Verbose)
        % exp(-746)=0
     if(-l<.00135), disp('Underflow likely.');end
    end
    l=max(C(:))*(-l); 
end 

if(opts.Verbose)
    fprintf('Lambda %g\n', l)
end
 ms=round(linspace(0,n,M+1));
 [~,ix]=sort(x);

 K=exp(-C/l); % Gibbs kernel
 W=zeros(M,k);

 maxIter=opts.MaxIter;
 tol=opts.MaxError;
 eta=opts.StepSize;       
 
 for m=1:M
     for i=1:k
       ii=ix(ms(m)+1:ms(m+1),i);nn=size(ii,1);
       % Accelerated gradient descent for OT
       % An Lai Gu

       % Cost matrix C, source and target weights mu and nu, approx.
       % parameter lambda, step size eta
       mu=1/n;nu=1/nn;

       phi=zeros(nn,1);
       z=zeros(nn,1);
       thta=1;
       grd=inf;
       t=0; % iteration count
       KK=K(:,ii);
       while~(t>maxIter || max(abs(grd))<tol)
           v=exp(phi/l);
           grd=v.*(KK'*(mu./(KK*v)))-nu;
          
           zplus=phi-eta*grd;zplus=zplus-mean(zplus);
           thta0=thta;
           thta=.5*(1+sqrt(1+4*thta^2));
           phi=zplus+(thta0-1)/thta*(zplus-z);
           z=zplus;
           t=t+1;
       end
       if(opts.Verbose)
        fprintf('Partition %2d factor %2d maxiter %4d maxerr %g \n', m,i, t,  max(abs(grd)) )
       end
       
      W(m,i)=-(sum(mu.*max(phi'-C(:,ii),[],2))-sum(nu.*phi));
      %
      v=exp(phi/l);
      u=mu./(KK*v);
      P=v'.*KK.*u; % plan
      %Q=(ones(n,1)*mu)*(ones(1,nn)*nu); % product plan
      f=l*log(u./mu);g=l*log(v./nu);
      % Sinkhorn primal
      %D(m,i)=sum(P.*(C(:,ii)+l*(log(P./(mu*nu))-1)),'all');
      % Sinkhorn dual
      D(m,i)=mean(f)+mean(g)-l*sum(P,'all');
     end
 end
end